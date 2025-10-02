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

from astrocalc.times import now
from fundamentals.logs import emptyLogger
log = emptyLogger()
mjd = now(
    log=log
).get_mjd()
mjdLimit = mjd - 15


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

        atExps, atTDOExps, psExps = get_exposures_in_maps_temporal_window(
            log=log, dbConn=dbConn, mmap=mmap, windowDays=14)

        if len(atExps.index) or len(psExps.index) or len(atTDOExps.index):

            print(len(atExps.index), len(atTDOExps.index), len(psExps.index))

            mapDF = get_the_map_as_healpix_dataframe(
                log=log, dbConn=dbConn, mapId=mmap["mapId"])
            match_exp_to_map_pixels(log=log, dbConn=dbConn, exps=atExps,
                                    mapId=mmap["mapId"], survey="atlas", nside=nside, pointingSideRA=5.46, pointingSideDec=5.46, mapDF=mapDF, settings=settings)
            match_exp_to_map_pixels(log=log, dbConn=dbConn, exps=atTDOExps,
                                    mapId=mmap["mapId"], survey="atlas", nside=nside, pointingSideRA=3.34096, pointingSideDec=2.22451556, mapDF=mapDF, settings=settings)
            match_exp_to_map_pixels(log=log, dbConn=dbConn, exps=psExps,
                                    mapId=mmap["mapId"], survey="ps", nside=nside, pointingSideRA=0.4, pointingSideDec=0.4, mapDF=mapDF, settings=settings)

        if index > 0:
            # Cursor up one line and clear line
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")

        percent = (float(index + 1) / float(len(maps))) * 100.
        print(f'{index+1}/{len(maps)} ({percent:1.1f} done)')

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
        select primaryId as mapId, map, mjd_obs from alerts where map is not null and significant = 1 and primaryId = 12045;
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

    - `atExps` -- atlas exposures as pandas dataframe (not tdo)
    - `atTDOExps` -- atlas tdo exposures as pandas dataframe
    - `psExps` -- panstarrs exposures as pandas dataframe
    """
    log.debug('starting the ``get_exposures_in_maps_temporal_window`` function')

    from fundamentals.mysql import readquery
    import pandas as pd
    start = mmap["mjd_obs"]

    sqlQuery = f"""
        SELECT primaryId as expname, raDeg, decDeg, mjd, mjd-{start} as 'mjd_t0' FROM lvk.exp_atlas where mjd > {start} and mjd < {start}+{windowDays} and (processed = 0 or mjd > {mjdLimit}) and expname not like "05%" order by mjd asc;
    """

    atExps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )
    atExps = pd.DataFrame(atExps)

    sqlQuery = f"""
        SELECT primaryId as expname, raDeg, decDeg, mjd, mjd-{start} as 'mjd_t0' FROM lvk.exp_atlas where mjd > {start} and mjd < {start}+{windowDays} and (processed = 0 or mjd > {mjdLimit}) and expname  like "05%" order by mjd asc;
    """

    atTDOExps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )
    atTDOExps = pd.DataFrame(atTDOExps)

    sqlQuery = f"""
        SELECT e.primaryId as expname, raDeg, decDeg, stacked, mjd, mjd-{start} as 'mjd_t0' FROM lvk.exp_ps e, lvk.ps1_skycell_map m where e.skycell=m.skycell_id and mjd > {start} and mjd < {start}+{windowDays} and (processed = 0 or mjd > {mjdLimit}) order by mjd asc;
    """
    psExps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )
    psExps = pd.DataFrame(psExps)

    log.debug('completed the ``get_exposures_in_maps_temporal_window`` function')
    return atExps, atTDOExps, psExps


def match_exp_to_map_pixels(
        dbConn,
        log,
        exps,
        mapId,
        survey,
        nside,
        pointingSideRA,
        pointingSideDec,
        mapDF,
        settings):
    """*match the exposures to the event maps pixels and record in lvk database*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger
    - `exps` -- exposures dataframe
    - `mapId` -- primaryId in the alerts table
    - `survey` -- atlas or ps
    - `nside` -- size of healpix pixel being used
    - ``pointingSideRA`` -- the width of the side of the square exposure, in degrees
    - ``pointingSideDec`` -- the height of the side of the square exposure, in degrees
    - ``mapDF`` -- the map as a dataframe ... one row per ipix
    """
    log.debug('starting the ``match_exp_to_map_pixels`` function')

    import pandas as pd
    from fundamentals.mysql import insert_list_of_dictionaries_into_database_tables
    import healpy as hp
    import numpy as np

    pixelArea = float(hp.nside2pixarea(nside, degrees=True))

    if not len(exps.index):
        return

    tmpDf = exps.copy()

    tmpDf["decCorner1"] = tmpDf["decDeg"] - pointingSideDec / 2
    tmpDf["decCorner2"] = tmpDf["decDeg"] + pointingSideDec / 2

    tmpDf.loc[(tmpDf['decCorner1'] > 90.), 'decCorner1'] = 180. - \
        tmpDf.loc[(tmpDf['decCorner1'] > 90.)]
    tmpDf.loc[(tmpDf['decCorner1'] < -90.), 'decCorner1'] = - \
        180. - tmpDf.loc[(tmpDf['decCorner1'] < -90.)]
    tmpDf.loc[(tmpDf['decCorner2'] > 90.), 'decCorner2'] = 180. - \
        tmpDf.loc[(tmpDf['decCorner2'] > 90.)]
    tmpDf.loc[(tmpDf['decCorner2'] < -90.), 'decCorner2'] = - \
        180. - tmpDf.loc[(tmpDf['decCorner2'] < -90.)]

    tmpDf['decMean'] = tmpDf[['decCorner1', 'decCorner2']].mean(axis=1)

    tmpDf['raCorner1'] = tmpDf["raDeg"] - \
        (pointingSideRA / 2) / np.cos(np.deg2rad(tmpDf['decMean']))
    tmpDf['raCorner2'] = tmpDf["raDeg"] + \
        (pointingSideRA / 2) / np.cos(np.deg2rad(tmpDf['decMean']))

    tmpDf.loc[(tmpDf['raCorner1'] > 360.), 'raCorner1'] = 720. - \
        tmpDf.loc[(tmpDf['raCorner1'] > 360.)]
    tmpDf.loc[(tmpDf['raCorner1'] < 0.), 'raCorner1'] = 360. + \
        tmpDf.loc[(tmpDf['raCorner1'] < 0.)]
    tmpDf.loc[(tmpDf['raCorner2'] > 360.), 'raCorner2'] = 720. - \
        tmpDf.loc[(tmpDf['raCorner2'] > 360.)]
    tmpDf.loc[(tmpDf['raCorner2'] < 0.), 'raCorner2'] = 360. + \
        tmpDf.loc[(tmpDf['raCorner2'] < 0.)]

    one = hp.ang2vec(tmpDf['raCorner1'].values,
                     tmpDf['decCorner1'].values, lonlat=True)
    two = hp.ang2vec(tmpDf['raCorner2'].values,
                     tmpDf['decCorner1'].values, lonlat=True)
    three = hp.ang2vec(tmpDf['raCorner1'].values,
                       tmpDf['decCorner2'].values, lonlat=True)
    four = hp.ang2vec(tmpDf['raCorner2'].values,
                      tmpDf['decCorner2'].values, lonlat=True)

    bigList = []
    # 1,2,4,3 IS NOT A BUG ... HEALPY NEEDS THIS ORDER
    bigList[:] = [[o, t, f, th] for o, t, th, f in zip(one, two, three, four)]
    tmpDf['corners'] = bigList

    ipix = []
    ipix[:] = [hp.query_polygon(nside, np.array(c), nest=True)
               for c in bigList]

    if pointingSideDec < 1:
        exps["ipix"] = ipix
        print(ipix)

    exps.dropna(axis='index', how='any', subset=['ipix'], inplace=True)
    # EXPLODE THE DF TO ONE ROW PER IPIX
    exps = exps.explode('ipix')

    from tabulate import tabulate
    print(tabulate(exps.head(1000), headers='keys', tablefmt='psql'))

    expMapDf = pd.merge(exps, mapDF, how='inner', on=['ipix'])
    expMapDf['area'] = pixelArea

    print(tabulate(expMapDf.head(1000), headers='keys', tablefmt='psql'))
    print(pointingSideDec)

    # SORT BY COLUMN NAME
    expMapDf.sort_values(['mjd'], inplace=True)

    firstIpixCoverage = expMapDf.drop_duplicates(subset=['ipix']).copy()

    # RENAME SOME INDIVIDUALLY
    firstIpixCoverage[f"exp_{survey}_id"] = firstIpixCoverage["expname"]
    firstIpixCoverage = firstIpixCoverage[[
        f'exp_{survey}_id', 'mapId', 'ipix']]

    print(len(firstIpixCoverage.index))
    print(survey)
    print(survey)
    firstIpixCoverage = firstIpixCoverage.to_dict('records')

    # USE dbSettings TO ACTIVATE MULTIPROCESSING - INSERT LIST OF DICTIONARIES INTO DATABASE
    insert_list_of_dictionaries_into_database_tables(
        dbConn=dbConn,
        log=log,
        dictList=firstIpixCoverage,
        dbTableName="alert_pixels_128",
        dateModified=False,
        dateCreated=False,
        batchSize=200000,
        replace=False,
        dbSettings=settings["database settings"]
    )

    # GENERATE EXPOSURE STATS
    if survey == "ps":
        expStats = expMapDf.groupby(f"expname").agg({'prob': 'sum', 'distmu': 'mean', 'distsigma': 'mean', 'distnorm': 'mean',
                                                     'mjd': 'first', 'mjd_t0': 'first', 'area': 'sum', 'mapId': 'first', 'stacked': 'first'}).reset_index()
    else:
        expStats = expMapDf.groupby(f"expname").agg({'prob': 'sum', 'distmu': 'mean', 'distsigma': 'mean',
                                                     'distnorm': 'mean', 'mjd': 'first', 'mjd_t0': 'first', 'area': 'sum', 'mapId': 'first'}).reset_index()

    # RENAME COLUMNS
    renames = {
        "prob": "prob_90",
        "distmu": "distmu_90",
        "distsigma": "distsigma_90",
        "distnorm": "distnorm_90",
        "area": "area_90",
        "expname": "expId"
    }
    expStats.rename(columns=renames, inplace=True)

    # REPLACE NANS
    mask = (expStats['distsigma_90'].isnull())
    expStats.loc[mask, "distmu_90"] = None
    expStats.loc[mask, "distnorm_90"] = None
    expStats.loc[mask, "distsigma_90"] = None
    expStats = expStats.replace({np.nan: None})

    expStats = expStats.to_dict('records')

    # USE dbSettings TO ACTIVATE MULTIPROCESSING - INSERT LIST OF DICTIONARIES INTO DATABASE
    insert_list_of_dictionaries_into_database_tables(
        dbConn=dbConn,
        log=log,
        dictList=expStats,
        dbTableName=f"exp_{survey}_alert_map_matches",
        dateModified=False,
        dateCreated=False,
        batchSize=200000,
        replace=False,
        dbSettings=settings["database settings"]
    )

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
      `area_90` double DEFAULT NULL,
      `mjd` double DEFAULT NULL,
      `mjd_t0` double DEFAULT NULL,
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
      `area_90` double DEFAULT NULL,
      `mjd` double DEFAULT NULL,
      `mjd_t0` double DEFAULT NULL,
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


def get_the_map_as_healpix_dataframe(
        dbConn,
        log,
        mapId):
    """*summary of function*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger,
        mapId

    **Usage:**

    ```eval_rst
    .. todo::

            add usage info
            create a sublime snippet for usage
    ```

    ```python
    usage code 
    ```           
    """
    log.debug('starting the ``get_the_map_as_healpix_dataframe`` function')

    from fundamentals.mysql import readquery
    import pandas as pd

    sqlQuery = f"""
        select * from alert_pixels_128 where mapId = {mapId};
    """
    rows = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )
    mapDF = pd.DataFrame(rows)

    log.debug('completed the ``get_the_map_as_healpix_dataframe`` function')
    return mapDF

# use the tab-trigger below for new function
# xt-def-function


if __name__ == '__main__':
    main()
