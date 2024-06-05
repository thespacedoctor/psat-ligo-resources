#!/usr/bin/env python
# encoding: utf-8
"""
*Convert lvk skymaps to single resolution nside=128 healpix grids and import them into the lvk database*

:Author:
    David Young

:Date Created:
    May 31, 2023


Before using this script in a gocart conda environment, you will need to:

```bash
conda install pymysql -c conda-forge
```

You will also need to add database credientials to the gocart.yaml file.

Usage:
    lvk_map_pixels_to_db 

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

    maps = list_maps_to_be_imported(dbConn=dbConn, log=log)
    import_maps_to_db(log=log, dbConn=dbConn, maps=maps, settings=settings)

    return


def list_maps_to_be_imported(
        dbConn,
        log):
    """*Generate a list of maps still to be imported into the database*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger  
    """
    log.debug('starting the ``list_maps_to_be_imported`` function')

    from fundamentals.mysql import writequery
    from fundamentals.mysql import readquery

    sqlQuery = f"""CREATE TABLE IF NOT EXISTS `alert_pixels_128` (
      `mapId` int(11) DEFAULT NULL,
      `ipix` int(11) DEFAULT NULL,
      `cumprob` double DEFAULT NULL,
      `prob` double DEFAULT NULL,
      `exp_ps_id` bigint(11) DEFAULT NULL,
      `exp_atlas_id` bigint(11) DEFAULT NULL,
      `distmu` double DEFAULT NULL,
      `distsigma` double DEFAULT NULL,
      `distnorm` double DEFAULT NULL,
      UNIQUE KEY `mapid_ipix` (`mapId`,`ipix`)
    ) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_swedish_ci;
    """
    writequery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn
    )

    sqlQuery = f"""
        select primaryId, map from alerts where significant = 1 and primaryId not in (select distinct mapId from alert_pixels_128 where mapId is not null);
    """
    maps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )

    log.debug('completed the ``list_maps_to_be_imported`` function')
    return maps


def convert_map_to_list_of_dicts(
        log,
        mapId,
        mapPath,
        nside=128):
    """*convert an lvk map to a list of dictionaries, each row containing a nside=128 pixel*

    **Key Arguments:**

    - `log` -- logger
    - `mapId` -- primaryId from the alert table
    - `mapPath` -- path to the map
    - `nside` -- resolution of the map pixels to return

    """
    log.debug('starting the ``convert_map_to_list_of_dicts`` function')

    from gocart.commonutils import flatten_healpix_map
    import numpy as np

    if not os.path.exists(mapPath):
        log.error(f'map does not exist at path : `{mapPath}`' % locals())
        return None

    skymap = flatten_healpix_map(
        log=log,
        mapPath=mapPath,
        nside=nside
    )

    print(skymap.columns())

    # SORT BY PROB, CALCULATE CUMULATIVE PROB AND RESORT BY INDEX
    skymap.sort_values(["PROB"],
                       ascending=[False], inplace=True)
    skymap["CUMPROB"] = np.cumsum(skymap['PROB'])

    # DROP PIXELS OUTSIDE 90 CONTOUR
    mask = (skymap["CUMPROB"] <= 0.9)
    skymap = skymap.loc[mask]

    # SET 0 DIST TO NULL
    mask = (skymap["DISTSIGMA"].isnull())
    skymap.loc[mask, "DISTMU"] = None

    # REMOVE COLUMN FROM DATA FRAME
    # try:
    #     skymap.drop(columns=['DISTMU', 'DISTSIGMA', 'DISTNORM'], inplace=True)
    # except:
    #     pass

    skymap["mapId"] = mapId
    skymap["ipix"] = skymap.index

    listOfDicts = skymap.to_dict('records')

    log.debug('completed the ``convert_map_to_list_of_dicts`` function')
    return listOfDicts


def import_maps_to_db(
        dbConn,
        log,
        maps,
        settings):
    """*import map pixels within the 90% contours to the lvk database*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger
    - `maps` -- list of dictionaries (containing `mapId` and `map` path)       
    - `settings` -- dict of gocart settings 
    """
    log.debug('starting the ``import_maps_to_db`` function')

    nside = 128
    from fundamentals.mysql import insert_list_of_dictionaries_into_database_tables

    mapPixelTableName = f"alert_pixels_{nside}"

    for m in maps:
        listOfDicts = convert_map_to_list_of_dicts(log=log, mapId=m["primaryId"], mapPath=m["map"], nside=nside)

        # USE dbSettings TO ACTIVATE MULTIPROCESSING - INSERT LIST OF DICTIONARIES INTO DATABASE
        insert_list_of_dictionaries_into_database_tables(
            dbConn=dbConn,
            log=log,
            dictList=listOfDicts,
            dbTableName=mapPixelTableName,
            uniqueKeyList=["mapId", "ipix"],
            dateModified=False,
            dateCreated=False,
            batchSize=2500,
            replace=True,
            dbSettings=settings["database settings"]
        )

    log.debug('completed the ``import_maps_to_db`` function')
    return None


if __name__ == '__main__':
    main()
