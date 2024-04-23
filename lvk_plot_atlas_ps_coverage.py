#!/usr/bin/env python
# encoding: utf-8
"""
*Plot the ATLAS and PanSTARRS coverage of each event map*

:Author:
    David Young

:Date Created:
    June 5, 2023

Before using this script in a gocart conda environment, you will need to:

```bash
conda install pymysql -c conda-forge
```

You will also need to add database credentials to the gocart.yaml file.

Usage:
    lvk_plot_atlas_ps_coverage [<daysAgo>]

Options:
    <daysAgo>             only replot events from last N days

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
import yaml
import matplotlib.path as mpath
import matplotlib.patches as patches
from gocart.convert import aitoff
from astropy.time import Time
from datetime import datetime
import pandas as pd
from tabulate import tabulate


def main(arguments=None):
    """
    *The main function used when ``lvk_plot_atlas_ps_coverage.py`` is run as a single script from the cl*
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

    import pandas as pd
    nside = 128
    pixelArea = float(hp.nside2pixarea(nside, degrees=True))
    maps = list_maps_to_be_plotted(dbConn=dbConn, log=log, daysAgo=a["daysAgo"])

    print(f"Generating {len(maps)} x 3 plots")
    count = len(maps)

    for index, mmap in enumerate(maps):

        exists = os.path.exists(mmap["map"])
        if not exists:
            print(f"The map '{mmap['map']}' does not exist on this file system")
            continue

        if index > 1:
            # Cursor up one line and clear line
            sys.stdout.write("\x1b[1A\x1b[2K")

        percent = (float(index) / float(count)) * 100.
        print(f'{index}/{count} ({percent:1.1f}% done)')

        mapMjd = mmap["mjd_obs"]

        # NOW WRITE OUT ALL EXPOSURES FOR ATLAS AND PS
        atlasExps, atlasStats = get_atlas_exposures_covering_map(log=log, dbConn=dbConn, mapId=mmap["mapId"], pixelArea=pixelArea, mjdUpper=mapMjd + 14)
        psExps, psStats = get_ps_skycells_covering_map(log=log, dbConn=dbConn, mapId=mmap["mapId"], pixelArea=pixelArea, mjdUpper=mapMjd + 14)

        outputFolder = os.path.dirname(mmap["map"])
        df = pd.DataFrame(atlasExps)
        df.to_csv(outputFolder + "/atlas_exposures.csv", index=False)
        df = pd.DataFrame(psExps)
        df.to_csv(outputFolder + "/ps_skycells.csv", index=False)

        continue

        coverageStats = []
        for rangeDays in [1, 3, 7]:
            atlasExps, atlasStats = get_atlas_exposures_covering_map(log=log, dbConn=dbConn, mapId=mmap["mapId"], mjdUpper=mapMjd + rangeDays, pixelArea=pixelArea)
            psExps, psStats = get_ps_skycells_covering_map(log=log, dbConn=dbConn, mapId=mmap["mapId"], mjdUpper=mapMjd + rangeDays, pixelArea=pixelArea)

            atlasStats["days since event"] = rangeDays
            psStats["days since event"] = rangeDays
            atlasStats["survey"] = "atlas"
            psStats["survey"] = "panstarrs"
            coverageStats.append(atlasStats)
            coverageStats.append(psStats)

            # GRAB META
            try:
                yamlFilePath = outputFolder + "/meta.yaml"
                with open(yamlFilePath, 'r') as stream:
                    meta = yaml.safe_load(stream)
            except:
                meta = {}

            atlasPatches = get_patches(log=log, exposures=atlasExps, pointingSide=5.46)
            psPatches = get_patches(log=log, exposures=psExps, pointingSide=0.4)

            converter = aitoff(
                log=log,
                mapPath=mmap["map"],
                outputFolder=outputFolder,
                settings=settings,
                plotName=f"atlas_coverage_{rangeDays}d.png",
                meta=meta,
                patches=atlasPatches,
                patchesColor="#d33682",
                patchesLabel=" ATLAS Exposure"
            )
            converter.convert()

            converter = aitoff(
                log=log,
                mapPath=mmap["map"],
                outputFolder=outputFolder,
                settings=settings,
                plotName=f"ps_coverage_{rangeDays}d.png",
                meta=meta,
                patches=psPatches,
                patchesColor="#859900",
                patchesLabel=" PanSTARRS Skycell"
            )
            converter.convert()

        coverageStats = pd.DataFrame(coverageStats)
        # SORT BY COLUMN NAME
        coverageStats.sort_values(['survey', 'days since event'],
                                  ascending=[True, True], inplace=True)

        if "ALERT" in meta:
            header = f"# {meta['ALERT']['superevent_id']}, {meta['ALERT']['alert_type']} Alert (issued {meta['ALERT']['time_created'].replace('Z','')} UTC)\n"
        else:
            header = f"# archive map\n"
        coverageStats = tabulate(coverageStats, headers='keys', tablefmt='psql', showindex=False)
        with open(outputFolder + "/map_coverage.txt", "w") as myFile:
            myFile.write(header)
            myFile.write(coverageStats)

    return


def list_maps_to_be_plotted(
        dbConn,
        log,
        daysAgo=False):
    """*Generate a list of maps needing maps to be plotted*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger  
    """
    log.debug('starting the ``list_maps_to_be_plotted`` function')

    extra = ""
    if daysAgo:
        utcnow = datetime.utcnow()
        mjdnow = Time([utcnow], scale='utc').mjd[0]
        extra = f"and mjd_obs > {mjdnow} - {daysAgo}"

    from fundamentals.mysql import readquery
    sqlQuery = f"""
        select primaryId as mapId, map, mjd_obs from alerts where map is not null and significant = 1 {extra};
    """

    maps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )

    log.debug('completed the ``list_maps_to_be_plotted`` function')
    return maps


def get_atlas_exposures_covering_map(
        log,
        dbConn,
        mapId,
        pixelArea,
        mjdUpper=700000000):
    """*Get all of the atlas exposures covering map*

    **Key Arguments:**

    - `log` -- logger
    - `dbConn` -- mysql database connection
    - `mapId` -- the primaryId of the map in database   
    - `pixelArea` -- healpix nside pixel area in square deg
    - `mjdUpper` -- return exposures taken before this mjd
    """
    log.debug('starting the ``get_atlas_exposures_covering_map`` function')

    from fundamentals.mysql import readquery
    sqlQuery = f"""
        SELECT DISTINCT
            mjd,
            expname,
            raDeg,
            decDeg,
            exp_time,
            filter,
            limiting_magnitude
        FROM
            exp_atlas e,
            alert_pixels_128 p
        WHERE
            p.mapId = {mapId}
                AND e.primaryId = p.exp_atlas_id
                and e.mjd < {mjdUpper}
        ORDER BY mjd;
    """
    atlasExps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )

    sqlQuery = f"""
        select count(*) as count, sum(p.prob)*100 as prob, count(*)*{pixelArea} as area from exp_atlas e, alert_pixels_128 p where p.mapId = {mapId} and e.primaryId = p.exp_atlas_id and e.mjd < {mjdUpper};
    """
    pixels = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )

    if pixels[0]['count'] == 0:
        pixels[0]['prob'] = 0.
        pixels[0]['area'] = 0.

    stats = {
        "prob. coverage (%)": float(f"{pixels[0]['prob']:0.2f}"),
        "90% area coverage (squ.deg.)": float(f"{pixels[0]['area']:0.2f}"),
    }

    log.debug('completed the ``get_atlas_exposures_covering_map`` function')
    return atlasExps, stats


def get_ps_skycells_covering_map(
        log,
        dbConn,
        mapId,
        pixelArea,
        mjdUpper=700000000):
    """*Get all of the panstarrs skycells covering map*

    **Key Arguments:**

    - `log` -- logger
    - `dbConn` -- mysql database connection
    - `mapId` -- the primaryId of the map in database 
    - `pixelArea` -- healpix nside pixel area in square deg    
    - `mjdUpper` -- return exposures taken before this mjd     
    """
    log.debug('starting the ``get_ps_skycells_covering_map`` function')

    from fundamentals.mysql import readquery
    sqlQuery = f"""
        SELECT distinct
            mjd,
            imageID,
            skycell,
            raDeg,
            decDeg,
            exp_time,
            filter,
            limiting_mag
        FROM
            exp_ps e,
            ps1_skycell_map s,
            alert_pixels_128 p
        WHERE
            s.skycell_id = e.skycell
                AND e.primaryId = p.exp_ps_id
                AND p.mapId = {mapId} and e.mjd < {mjdUpper}
        ORDER BY mjd ASC;
    """

    psExps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )

    sqlQuery = f"""
        select count(*) as count, sum(p.prob)*100 as prob, count(*)*{pixelArea} as area from exp_ps e, ps1_skycell_map s,alert_pixels_128 p where s.skycell_id=e.skycell and e.primaryId = p.exp_ps_id and p.mapId = {mapId} and e.mjd < {mjdUpper};
    """
    pixels = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )

    if pixels[0]['count'] == 0:
        pixels[0]['prob'] = 0.
        pixels[0]['area'] = 0.

    stats = {
        "prob. coverage (%)": float(f"{pixels[0]['prob']:0.2f}"),
        "90% area coverage (squ.deg.)": float(f"{pixels[0]['area']:0.2f}"),
    }

    log.debug('completed the ``get_ps_skycells_covering_map`` function')
    return psExps, stats


def get_patches(
        log,
        exposures,
        pointingSide):
    """*Convert the exposures/skycells to matplotlib patches to be added to the plot*

    **Key Arguments:**

    - `log` -- logger
    - `exposures` -- atlas or panstarrs exposures or skycells
    - `pointingSide` -- the pointing side in degrees        
    """
    log.debug('starting the ``get_patches`` function')

    expPatches = []
    for e in exposures:

        raDeg = e['raDeg']
        raDeg = -raDeg + 180
        if raDeg > 180.:
            raDeg -= 360
        raDeg = -raDeg
        decDeg = e['decDeg']

        deltaDeg = pointingSide / 2
        if decDeg < 0:
            deltaDeg = -deltaDeg

        widthRadTop = np.deg2rad(pointingSide) / np.cos(np.deg2rad(decDeg + deltaDeg))
        widthRadBottom = np.deg2rad(pointingSide) / np.cos(np.deg2rad(decDeg - deltaDeg))
        heightRad = np.deg2rad(pointingSide)
        llx = -(np.deg2rad(raDeg) - widthRadBottom / 2)
        lly = np.deg2rad(decDeg) - (heightRad / 2)
        ulx = -(np.deg2rad(raDeg) - widthRadTop / 2)
        uly = np.deg2rad(decDeg) + (heightRad / 2)
        urx = -(np.deg2rad(raDeg) + widthRadTop / 2)
        ury = uly
        lrx = -(np.deg2rad(raDeg) + widthRadBottom / 2)
        lry = lly
        Path = mpath.Path
        path_data = [
            (Path.MOVETO, [llx, lly]),
            (Path.LINETO, [ulx, uly]),
            (Path.LINETO, [urx, ury]),
            (Path.LINETO, [lrx, lry]),
            (Path.CLOSEPOLY, [llx, lly])
        ]
        codes, verts = zip(*path_data)
        path = mpath.Path(verts, codes)
        expPatches.append(patches.PathPatch(path))

    log.debug('completed the ``get_patches`` function')
    return expPatches


if __name__ == '__main__':
    main()
