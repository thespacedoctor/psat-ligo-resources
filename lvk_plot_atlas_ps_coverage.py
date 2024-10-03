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
    lvk_plot_atlas_ps_coverage -g [<gid>]

Options:
    <daysAgo>             only replot events from last N days
    <gid>                 the gravity event ID

    -h, --help            show this help message
    -v, --version         show version
    -s, --settings        the settings file
    -g, --gevent          the ID of a gravity event
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
    maps = list_maps_to_be_plotted(dbConn=dbConn, log=log, daysAgo=a["daysAgo"], gid=a["gid"])

    print(f"Generating {len(maps)} x 4 plots")
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
        atlasExps, atlasStats = get_atlas_exposures_covering_map(log=log, dbConn=dbConn, mapId=mmap["mapId"], pixelArea=pixelArea, mjdLower=mapMjd, mjdUpper=mapMjd + 14, allSkycells=True)
        psExps, psStats = get_ps_skycells_covering_map(log=log, dbConn=dbConn, mapId=mmap["mapId"], pixelArea=pixelArea, mjdLower=mapMjd, mjdUpper=mapMjd + 14, allSkycells=True)

        outputFolder = os.path.dirname(mmap["map"])
        df = pd.DataFrame(atlasExps)
        df = df.round({'mjd': 6, 'mjd_t0': 6, 'limiting_magnitude': 2, 'raDeg': 6, 'decDeg': 6, 'area_90': 5, 'prob_90': 5, 'distmu_90': 2, 'distsigma_90': 2, 'distnorm_90': 7})
        df.rename(columns={"limiting_magnitude": "mag5sig"}, inplace=True)
        df.to_csv(outputFolder + "/atlas_exposures.csv", index=False)
        df = pd.DataFrame(psExps)
        df = df.round({'mjd': 6, 'mjd_t0': 6, 'limiting_magnitude': 2, 'raDeg': 6, 'decDeg': 6, 'area_90': 5, 'prob_90': 5, 'distmu_90': 2, 'distsigma_90': 2, 'distnorm_90': 7})

        if len(df.index):
            mask = (df["stacked"] == 1)
            this = df.loc[mask].copy()
            this.drop(columns=['stacked'], inplace=True)
            this.to_csv(outputFolder + "/ps_skycells_stacks.csv", index=False)
            this = df.loc[~mask].copy()
            this.drop(columns=['stacked'], inplace=True)
            this.to_csv(outputFolder + "/ps_skycells_warps.csv", index=False)
        else:
            df.to_csv(outputFolder + "/ps_skycells_stacks.csv", index=False)
            df.to_csv(outputFolder + "/ps_skycells_warps.csv", index=False)

        coverageStats = []
        for rangeDays in [1, 3, 7, 14]:
            atlasExps, atlasStats = get_atlas_exposures_covering_map(log=log, dbConn=dbConn, mapId=mmap["mapId"], mjdLower=mapMjd, mjdUpper=mapMjd + rangeDays, pixelArea=pixelArea)
            psExps, psStats = get_ps_skycells_covering_map(log=log, dbConn=dbConn, mapId=mmap["mapId"], mjdLower=mapMjd, mjdUpper=mapMjd + rangeDays, pixelArea=pixelArea)

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

            if rangeDays == 14:
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

        from tabulate import tabulate
        print(tabulate(coverageStats, headers='keys', tablefmt='psql'))

        if "ALERT" in meta:
            header = f"# {meta['ALERT']['superevent_id']}, {meta['ALERT']['alert_type']} Alert (issued {meta['ALERT']['time_created'].replace('Z','')} UTC)\n"
        else:
            header = f"# archive map\n"
        coverageStats = tabulate(coverageStats, headers='keys', tablefmt='psql', showindex=False)
        with open(outputFolder + "/map_coverage.txt", "w") as myFile:
            myFile.write(header)
            myFile.write(coverageStats)

        with open(outputFolder + "/README.txt", "w") as myFile:
            content = readme_content()
            myFile.write(content)

    return


def list_maps_to_be_plotted(
        dbConn,
        log,
        daysAgo=False,
        gid=False):
    """*Generate a list of maps needing maps to be plotted*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger
    - `daysAgo` -- plot all events within the last 'daysAgo' days
    - ``gid`` -- single out a gravity event to plot
    """
    log.debug('starting the ``list_maps_to_be_plotted`` function')

    extra = ""
    if daysAgo and not gid:
        utcnow = datetime.utcnow()
        mjdnow = Time([utcnow], scale='utc').mjd[0]
        extra = f"and mjd_obs > {mjdnow} - {daysAgo}"

    if gid:
        extra = f"and superevent_id = '{gid}'"

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
        mjdLower,
        mjdUpper=700000000,
        allSkycells=False):
    """*Get all of the atlas exposures covering map*

    **Key Arguments:**

    - `log` -- logger
    - `dbConn` -- mysql database connection
    - `mapId` -- the primaryId of the map in database
    - `pixelArea` -- healpix nside pixel area in square deg
    - `mjdLower` -- the mjd of the event
    - `mjdUpper` -- return exposures taken before this mjd
    """
    log.debug('starting the ``get_atlas_exposures_covering_map`` function')

    from fundamentals.mysql import readquery

    if allSkycells:
        sqlQuery = f"""
            SELECT expname, m.mjd, m.mjd_t0, filter, exp_time, limiting_magnitude, raDeg, decDeg, area_90, prob_90, distmu_90, distsigma_90, distnorm_90 FROM lvk.exp_atlas_alert_map_matches m, lvk.exp_atlas e  where m.mapId = {mapId} and m.expId=e.primaryId order by m.mjd;
        """
    else:
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
                    and e.mjd < {mjdUpper} and e.mjd > {mjdLower}
            ORDER BY mjd;
        """
    print(sqlQuery)
    atlasExps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )

    sqlQuery = f"""
        select count(*) as count, sum(p.prob)*100 as prob, count(*)*{pixelArea} as area from exp_atlas e, alert_pixels_128 p where p.mapId = {mapId} and e.primaryId = p.exp_atlas_id and e.mjd < {mjdUpper} and e.mjd > {mjdLower};
    """
    pixels = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )

    print(sqlQuery)
    print("\n\n")
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
        mjdLower,
        mjdUpper=700000000,
        allSkycells=False):
    """*Get all of the panstarrs skycells covering map*

    **Key Arguments:**

    - `log` -- logger
    - `dbConn` -- mysql database connection
    - `mapId` -- the primaryId of the map in database
    - `pixelArea` -- healpix nside pixel area in square deg
    - `mjdLower` -- the mjd of the event
    - `mjdUpper` -- return exposures taken before this mjd
    - `allSkycells` -- return all skycells
    """
    log.debug('starting the ``get_ps_skycells_covering_map`` function')

    from fundamentals.mysql import readquery

    if allSkycells:
        sqlQuery = f"""
            SELECT
                imageid AS 'expname',
                skycell,
                m.mjd,
                m.mjd_t0,
                filter,
                exp_time,
                limiting_mag AS 'limiting_magnitude',
                raDeg,
                decDeg,
                area_90,
                prob_90,
                distmu_90,
                distsigma_90,
                distnorm_90,
                m.stacked
            FROM
                lvk.exp_ps_alert_map_matches m,
                lvk.exp_ps e,
                lvk.ps1_skycell_map s
            WHERE
                m.mapId = {mapId} AND m.expId = e.primaryId
                    AND s.skycell_id = skycell
                    and e.mjd < {mjdUpper} and e.mjd > {mjdLower}
            ORDER BY m.mjd ASC;
        """
    else:
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
                    AND p.mapId = {mapId} and e.mjd < {mjdUpper} and e.mjd > {mjdLower}
            ORDER BY mjd ASC;
        """

    psExps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )

    sqlQuery = f"""
        select count(*) as count, sum(p.prob)*100 as prob, count(*)*{pixelArea} as area from exp_ps e, ps1_skycell_map s,alert_pixels_128 p where s.skycell_id=e.skycell and e.primaryId = p.exp_ps_id and p.mapId = {mapId} and e.mjd < {mjdUpper} and e.mjd > {mjdLower};
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


def readme_content():
    """*get the content for the readme file*
    ```           
    """

    content = """
# Exposure Exports

MJD-ordered exports of the ATLAS exposures and PS stack-warp and stack-stack skycells overlapping the 90% map regions are written to the files `atlas_exposures.csv`, `ps_skycells_warps.csv`and `ps_skycells_stacks.csv`.

The ATLAS exposure limiting magnitudes are reported at 5𝜎, whereas the PS warps are 3.5𝜎. Until May 2024, the PS stack limiting mags were reported at 5𝜎, but from May 2024 onward, they are reported at 3.5𝜎, in line with the warps.

Each high-significance event alert map is converted to a level 7 (nside = 128) HEALPix map. The equal-area map pixels are then sorted from the highest-to-lowest likelihood of containing the actual location of the gravity event. A cumulative probability is assigned to each pixel, starting with the highest probability pixel and summing up to 100% once the lowest probability pixel is reached. All pixels with a cumulative probability of ≤90% are cached in a database. Note that for well-constrained event maps, the fraction of pixels within the 90% contour region is tiny compared to the whole sky.

ATLAS exposures and PS skycells taken within 0-14 days of the gravity event are projected onto the level-7 HEALPix sky, and their overlapping pixels are calculated. This allows the matching of the exposures and skycells HEALpix pixels against those within the 90% contours of the event maps. Finally, this allows for the following values to be calculated for each exposure/skycell:

`area_90`: sky area with the 90% contour of the event map covered by the exposure/skycell.

`prob_90`: total probability covered within `area_90` 

`distmu_90`: the mean distance covered within `area_90` (Mpc). For many exposures/skycells, this will be blank as the distance recorded in the maps is infinite at these sky locations.

`distsigma_90`: the mean distance scale parameter covered within `area_90` (Mpc). For many exposures/skycells, this will be blank as the distance recorded in the maps is infinite at these sky locations.

`distnorm_90`: the mean distance normalisation coefficient (Mpc^-2).

The distance posterior (probability per distance interval) at a given location is: 
dp(r)/dr = distnorm_90 * Gaussian(distmu_90, distsigma_90)

See here for more info on distance parameters: https://arxiv.org/pdf/1605.04242
"""

    return content

# use the tab-trigger below for new function
# xt-def-function


if __name__ == '__main__':
    main()
