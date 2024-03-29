#!/usr/bin/env python
# encoding: utf-8
"""
*Match the ATLAS and PanSTARRS exposures to the lvk event maps within the lvk database*

:Author:
    David Young

:Date Created:
    June 2, 2023


Before using this script in a gocart conda environment, you will need to:

```bash
conda install pymysql -c conda-forge
```

You will also need to add database credientials to the gocart.yaml file.

Usage:
    lvk_match_exposures_to_maps 

Options:

    -h, --help            show this help message
    -v, --version         show version
    -s, --settings        the settings file
"""
################# GLOBAL IMPORTS ####################
import sys
import os
from fundamentals import tools
from os.path import expanduser
import healpy as hp
import numpy as np
from fundamentals.mysql import writequery


def main(arguments=None):
    """
    *The main function used when ``lvk_map_pixels_to_db.py`` is run as a single script from the cl*
    """

    # SETUP THE COMMAND-LINE UTIL SETTINGS
    su = tools(
        arguments=arguments,
        docString=__doc__,
        logLevel="DEBUG",
        options_first=False,
        projectName="gocart",
        defaultSettingsFile=True
    )
    arguments, settings, log, dbConn = su.setup()

    # UNPACK REMAINING CL ARGUMENTS USING `EXEC` TO SETUP THE VARIABLE NAMES
    # AUTOMATICALLY
    a = {}
    for arg, val in list(arguments.items()):
        if arg[0] == "-":
            varname = arg.replace("-", "") + "Flag"
        else:
            varname = arg.replace("<", "").replace(">", "")
        a[varname] = val
        if arg == "--dbConn":
            dbConn = val
            a["dbConn"] = val
        log.debug('%s = %s' % (varname, val,))

    nside = 128
    maps = list_maps_still_to_be_covered(dbConn=dbConn, log=log)
    for index, mmap in enumerate(maps):
        atExps, psExps = get_exposures_in_maps_temporal_window(log=log, dbConn=dbConn, mmap=mmap, windowDays=7)
        match_exp_to_map_pixels(log=log, dbConn=dbConn, exps=atExps, mapId=mmap["mapId"], survey="atlas", nside=nside, pointingSide=5.46)
        match_exp_to_map_pixels(log=log, dbConn=dbConn, exps=psExps, mapId=mmap["mapId"], survey="ps", nside=nside, pointingSide=0.4)

    sqlQuery = f"""update exp_atlas set processed = 1 where processed = 0"""
    writequery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn
    )
    sqlQuery = f"""update exp_ps set processed = 1 where processed = 0"""
    writequery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn
    )

    return


def calulate_exposure_healpix_ids(
        series,
        pointingSide,
        nside):
    """*calculate the covered healpix IDs for exposures in dataframe*

    **Key Arguments:**

    - ``series`` -- the dataframe row/series to apply work on
    - ``pointingSide`` -- the length of the side of the square exposure, in degrees
    - ``nside`` -- size of the healpix pixels to calculate       
    """

    # DETERMINE THE CORNERS INDIVIDUAL EXPOSURE AS MAPPED TO THE
    decCorners = [series["decDeg"] - pointingSide / 2,
                  series["decDeg"] + pointingSide / 2]

    corners = []

    for i, d in enumerate(decCorners):
        if d > 90.:
            decCorners[i] = 180. - d
        elif d < -90.:
            decCorners[i] = -180 - d

    meanDec = np.mean(decCorners)

    for d in decCorners:
        raCorners = [series["raDeg"] - (pointingSide / 2) / np.cos(np.deg2rad(meanDec)),
                     series["raDeg"] + (pointingSide / 2) / np.cos(np.deg2rad(meanDec))]
        for i, r in enumerate(raCorners):
            if r > 360.:
                raCorners[i] = 720. - r
            if r < 0.:
                raCorners[i] = 360. + r
        for r in raCorners:
            corners.append(hp.ang2vec(r, d, lonlat=True))

     # FLIP CORNERS 3 & 4 SO HEALPY UNDERSTANDS POLYGON SHAPE
    corners = [corners[0], corners[1],
               corners[3], corners[2]]

    # RETURN ALL HEALPIXELS IN EXPOSURE AREA
    try:
        series["ipixs"] = hp.query_polygon(nside, np.array(
            corners), nest=True)
    except:
        pass

    return series


def list_maps_still_to_be_covered(
        dbConn,
        log):
    """*Generate a list of maps with area still to be covered*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger  
    """
    log.debug('starting the ``list_maps_still_to_be_covered`` function')

    from fundamentals.mysql import readquery
    sqlQuery = f"""
        select primaryId as mapId, map, mjd_obs from alerts where map is not null and significant = 1;
    """
    maps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )

    log.debug('completed the ``list_maps_still_to_be_covered`` function')
    return maps


def get_exposures_in_maps_temporal_window(
        log,
        dbConn,
        mmap,
        windowDays):
    """*get all ATLAS and PS1 exposures taken within a fixed temporal window after the map event*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger

    **Return:**

    - `atExps` -- atlas exposures as pandas dataframe
    - `psExps` -- panstarrs exposures as pandas dataframe    
    """
    log.debug('starting the ``get_exposures_in_maps_temporal_window`` function')

    from fundamentals.mysql import readquery
    import pandas as pd
    start = mmap["mjd_obs"]

    sqlQuery = f"""
        SELECT primaryId as expname, raDeg, decDeg FROM lvk.exp_atlas where mjd > {start} and mjd < {start}+{windowDays} and processed = 0 order by mjd asc;
    """
    atExps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )
    atExps = pd.DataFrame(atExps)

    sqlQuery = f"""
        SELECT e.primaryId as expname, raDeg, decDeg FROM lvk.exp_ps e, lvk.ps1_skycell_map m where e.skycell=m.skycell_id and mjd > {start} and mjd < {start}+{windowDays} and processed = 0 order by mjd asc;
    """
    psExps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )
    psExps = pd.DataFrame(psExps)

    log.debug('completed the ``get_exposures_in_maps_temporal_window`` function')
    return atExps, psExps


def match_exp_to_map_pixels(
        dbConn,
        log,
        exps,
        mapId,
        survey,
        nside,
        pointingSide):
    """*match the exposures to the event maps pixels and record in lvk database*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger
    - `exps` -- exposures dataframe
    - `mapId` -- primaryId in the alerts table
    - `survey` -- atlas or ps
    - `nside` -- size of healpix pixel being used
    - ``pointingSide`` -- the length of the side of the square exposure, in degrees   
    """
    log.debug('starting the ``match_exp_to_map_pixels`` function')

    if not len(exps.index):
        return

    exps = exps.apply(calulate_exposure_healpix_ids, axis=1, pointingSide=pointingSide, nside=nside)
    exps.dropna(axis='index', how='any', subset=['ipixs'], inplace=True)

    # ONLY DO THIS FOR SMALL DATAFRAMES - THIS IS AN ANTIPATTERN
    for index, row in exps.iterrows():
        if len(row["ipixs"]):
            expName = row["expname"]

            ipixs = (",").join(row["ipixs"].astype(str))
            sqlQuery = f"""update alert_pixels_128 set exp_{survey}_id = '{expName}' where ipix in ({ipixs}) and exp_{survey}_id is null and mapId = {mapId}"""
            writequery(
                log=log,
                sqlQuery=sqlQuery,
                dbConn=dbConn
            )

    log.debug('completed the ``match_exp_to_map_pixels`` function')
    return None


if __name__ == '__main__':
    main()
