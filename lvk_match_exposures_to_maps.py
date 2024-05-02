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

    create_db_tables(dbConn=dbConn, log=log)

    nside = 128
    maps = list_maps_still_to_be_covered(dbConn=dbConn, log=log)
    for index, mmap in enumerate(maps):

        atExps, psExps = get_exposures_in_maps_temporal_window(log=log, dbConn=dbConn, mmap=mmap, windowDays=14)

        print(f"{len(atExps.index)} ATLAS exposures")
        match_exp_to_map_pixels(log=log, dbConn=dbConn, exps=atExps, mapId=mmap["mapId"], survey="atlas", nside=nside, pointingSide=5.46)
        print(f"{len(psExps.index)} PS exposures")
        match_exp_to_map_pixels(log=log, dbConn=dbConn, exps=psExps, mapId=mmap["mapId"], survey="ps", nside=nside, pointingSide=0.4)

        if index > 1:
            # Cursor up one line and clear line
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")

        percent = (float(index) / float(len(mmap))) * 100.
        print(f'{index}/{len(maps)} ({percent:1.1f} done)')

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

    tmpDf = exps.copy()

    tmpDf["decCorner1"] = tmpDf["decDeg"] - pointingSide / 2
    tmpDf["decCorner2"] = tmpDf["decDeg"] + pointingSide / 2

    tmpDf.loc[(tmpDf['decCorner1'] > 90.), 'decCorner1'] = 180. - tmpDf.loc[(tmpDf['decCorner1'] > 90.)]
    tmpDf.loc[(tmpDf['decCorner1'] < -90.), 'decCorner1'] = -180. - tmpDf.loc[(tmpDf['decCorner1'] < -90.)]
    tmpDf.loc[(tmpDf['decCorner2'] > 90.), 'decCorner2'] = 180. - tmpDf.loc[(tmpDf['decCorner2'] > 90.)]
    tmpDf.loc[(tmpDf['decCorner2'] < -90.), 'decCorner2'] = -180. - tmpDf.loc[(tmpDf['decCorner2'] < -90.)]

    tmpDf['decMean'] = tmpDf[['decCorner1', 'decCorner2']].mean(axis=1)

    tmpDf['raCorner1'] = tmpDf["raDeg"] - (pointingSide / 2) / np.cos(np.deg2rad(tmpDf['decMean']))
    tmpDf['raCorner2'] = tmpDf["raDeg"] + (pointingSide / 2) / np.cos(np.deg2rad(tmpDf['decMean']))

    tmpDf.loc[(tmpDf['raCorner1'] > 360.), 'raCorner1'] = 720. - tmpDf.loc[(tmpDf['raCorner1'] > 360.)]
    tmpDf.loc[(tmpDf['raCorner1'] < 0.), 'raCorner1'] = 360. + tmpDf.loc[(tmpDf['raCorner1'] < 0.)]
    tmpDf.loc[(tmpDf['raCorner2'] > 360.), 'raCorner2'] = 720. - tmpDf.loc[(tmpDf['raCorner2'] > 360.)]
    tmpDf.loc[(tmpDf['raCorner2'] < 0.), 'raCorner2'] = 360. + tmpDf.loc[(tmpDf['raCorner2'] < 0.)]

    one = hp.ang2vec(tmpDf['raCorner1'].values, tmpDf['decCorner1'].values, lonlat=True)
    two = hp.ang2vec(tmpDf['raCorner2'].values, tmpDf['decCorner1'].values, lonlat=True)
    three = hp.ang2vec(tmpDf['raCorner1'].values, tmpDf['decCorner2'].values, lonlat=True)
    four = hp.ang2vec(tmpDf['raCorner2'].values, tmpDf['decCorner2'].values, lonlat=True)

    bigList = []
    # 1,2,4,3 IS NOT A BUG ... HEALPY NEEDS THIS ORDER
    bigList[:] = [[o, t, f, th] for o, t, th, f in zip(one, two, three, four)]
    tmpDf['corners'] = bigList

    ipix = []
    ipix[:] = [hp.query_polygon(nside, np.array(c), nest=True) for c in bigList]

    ipixStr = []
    ipixStr[:] = [(",").join(map(str, i)) for i in ipix]
    exps["ipixs"] = ipixStr

    exps.dropna(axis='index', how='any', subset=['ipixs'], inplace=True)

    # ONLY DO THIS FOR SMALL DATAFRAMES - THIS IS AN ANTIPATTERN
    print("QUERIES")
    sqlQueryList = []
    for index, row in exps.iterrows():
        if len(row["ipixs"]):
            expName = row["expname"]
            ipixs = row["ipixs"]
            sqlQuery = f"""update alert_pixels_128 set exp_{survey}_id = '{expName}' where ipix in ({ipixs}) and exp_{survey}_id is null and mapId = {mapId};"""
            sqlQueryList.append(sqlQuery)
    sqlQuery = ("\n".join(sqlQueryList))

    print("EXECUTING")
    writequery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn
    )
    print("DONE")

    log.debug('completed the ``match_exp_to_map_pixels`` function')
    return None


def create_db_tables(
        dbConn,
        log):
    """*create the required database tables*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger
    """
    log.debug('starting the ``create_db_tables`` function')

    from fundamentals.mysql import writequery
    sqlQuery = f"""CREATE TABLE IF NOT EXISTS `exp_atlas_alert_map_matches` (
      `expId` bigint(11) NOT NULL,
      `mapId` bigint(11) NOT NULL,
      `distmu_90` double DEFAULT NULL,
      `distsigma_90` double DEFAULT NULL,
      `distnorm_90` double DEFAULT NULL,
      `prob_90` double DEFAULT NULL,
      UNIQUE KEY `expid_mapid_uq` (`expId`,`mapId`)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;
    """
    writequery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
    )

    sqlQuery = f"""CREATE TABLE IF NOT EXISTS `exp_ps_alert_map_matches` (
      `expId` bigint(11) NOT NULL,
      `mapId` bigint(11) NOT NULL,
      `distmu_90` double DEFAULT NULL,
      `distsigma_90` double DEFAULT NULL,
      `distnorm_90` double DEFAULT NULL,
      `prob_90` double DEFAULT NULL,
      UNIQUE KEY `expid_mapid_uq` (`expId`,`mapId`)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;
    """
    writequery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
    )

    log.debug('completed the ``create_db_tables`` function')
    return None

# use the tab-trigger below for new function
# xt-def-function


if __name__ == '__main__':
    main()
