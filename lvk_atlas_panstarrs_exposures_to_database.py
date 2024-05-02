#!/usr/bin/env python
# encoding: utf-8
"""
*Import CSV files of ATLAS exposures and PanSTARRS skycells lvk database*

:Author:
    David Young

:Date Created:
    May 31, 2023


Before using this script in a gocart conda environment, you will need to:

```bash
conda install multiprocess
conda install pymysql -c conda-forge
```

You will also need to add database credientials to the gocart.yaml file.

Usage:
    lvk_atlas_panstarrs_exposures_to_database <pathToExports>

Options:
    pathToExports         path to the directory of CSV exports

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
    *The main function used when ``lvk_atlas_panstarrs_exposures_to_database.py`` is run as a single script from the cl*
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

    pathToExports = a["pathToExports"]
    home = expanduser("~")
    pathToExports = pathToExports.replace("~", home)

    create_tables_if_not_exist(log=log, dbConn=dbConn)
    tableNames, csvContents = read_csv_files(log=log, pathToExports=pathToExports)
    import_csv_content_to_database(log=log, dbConn=dbConn, tableNames=tableNames, csvContents=csvContents, settings=settings, pastDays=5)

    return


def read_csv_files(
        log,
        pathToExports):
    """*read in all exposure CSV files*

    **Key Arguments:**

    - ``log`` -- logger
    - ``pathToExports`` -- path to the CSV export folder

    **Return:**

    - ``tableNames`` -- a list of database tablenames
    - ``csvContents`` -- a list of the CSV connects (list of dictionaries). List equal in length to `tableNames`       
    """
    log.debug('starting the ``read_csv_files`` function')

    import csv

    # LIST CSV FILES
    allowedCsvs = [
        "ps13piWSExps.csv",
        "pso4WSExps.csv",
        "pso4SSExps.csv",
        "ps13piSSExps.csv",
        "atlas4Exps.csv"
    ]

    tableNames = []
    csvContents = []
    for d in os.listdir(pathToExports):
        filepath = os.path.join(pathToExports, d)
        if os.path.isfile(filepath) and os.path.splitext(filepath)[1] == ".csv":
            if d in allowedCsvs:
                if "SS" in d:
                    stacked = {"stacked": 1}
                else:
                    stacked = {}

                if "ps1" in d or "pso" in d or "ps2" in d:
                    tableNames.append("exp_ps")
                if "atlas" in d:
                    tableNames.append("exp_atlas")
                with open(filepath, 'r') as csvFile:
                    csvReader = csv.DictReader(csvFile, dialect='excel', delimiter=',', quotechar='"')
                    dictList = []
                    dictList[:] = [{**d, **stacked} for d in csvReader]
                    csvContents.append(dictList)
                csvFile.close()

    log.debug('completed the ``read_csv_files`` function')
    return tableNames, csvContents


def import_csv_content_to_database(
        log,
        dbConn,
        tableNames,
        csvContents,
        settings,
        pastDays):
    """*import csv content intp exp_ps and exp_atlas database talbes*

    **Key Arguments:**

    - `log` -- logger
    - `dbConn` -- mysql database connection
    - `tableNames` -- a list of database tablenames
    - `csvContents` -- a list of the CSV connects (list of dictionaries). List equal in length to `tableNames`
    - `settings` -- the settings dict    
    - `pastDays` -- import only the last N days of exposures  
    """
    log.debug('starting the ``import_csv_content_to_database`` function')

    from fundamentals.mysql import insert_list_of_dictionaries_into_database_tables
    from astropy.time import Time
    from datetime import datetime
    utcnow = datetime.utcnow()
    mjdnow = Time([utcnow], scale='utc').mjd[0]
    # USE dbSettings TO ACTIVATE MULTIPROCESSING - INSERT LIST OF DICTIONARIES INTO DATABASE

    if pastDays:
        print(f"Ingesting the last {pastDays} days worth of exposures")

    assert len(tableNames) == len(csvContents)

    for t, c in zip(tableNames, csvContents):

        if pastDays:
            c = [i for i in c if float(i["mjd"]) > mjdnow - pastDays]

        if not len(c):
            continue

        if "imageid" in c[0]:
            uniqueKeyList = ["imageid", "skycell"]
        elif "expname" in c[0]:
            uniqueKeyList = ["expname"]

        insert_list_of_dictionaries_into_database_tables(
            dbConn=dbConn,
            log=log,
            dictList=c,
            dbTableName=t,
            uniqueKeyList=uniqueKeyList,
            batchSize=200000,
            dateModified=False,
            dateCreated=False,
            replace=False,
            dbSettings=settings["database settings"]

        )

    log.debug('completed the ``import_csv_content_to_database`` function')
    return None


def create_tables_if_not_exist(
        dbConn,
        log):
    """*create the exp_ps and exp_atlas tables if they do not exist*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger     
    """
    log.debug('starting the ``create_tables_if_not_exist`` function')

    from fundamentals.mysql import writequery

    sqlQuery = f"""CREATE TABLE IF NOT EXISTS `exp_atlas` (
      `primaryId` bigint(20) NOT NULL AUTO_INCREMENT COMMENT 'An internal counter',
      `expname` varchar(50) DEFAULT NULL,
      `obj` varchar(20) DEFAULT NULL,
      `mjd` double DEFAULT NULL,
      `raDeg` double DEFAULT NULL,
      `decDeg` double DEFAULT NULL,
      `filter` varchar(5) DEFAULT NULL,
      `exp_time` float DEFAULT NULL,
      `limiting_magnitude` double(6,3) DEFAULT NULL,
      `dateCreated` datetime DEFAULT current_timestamp(),
      `dateLastModified` datetime DEFAULT current_timestamp() ON UPDATE current_timestamp(),
      `processed` TINYINT(1) NULL DEFAULT 0,
      PRIMARY KEY (`primaryId`),
      UNIQUE KEY `expname` (`expname`)
    ) ENGINE=MyISAM AUTO_INCREMENT=0 DEFAULT CHARSET=latin1 COLLATE=latin1_swedish_ci;
    """

    writequery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
    )

    sqlQuery = f"""CREATE TABLE `exp_ps` (
      `primaryId` bigint(20) NOT NULL AUTO_INCREMENT COMMENT 'An internal counter',
      `exp_time` float DEFAULT NULL,
      `filter` varchar(5) DEFAULT NULL,
      `imageid` bigint(11) DEFAULT NULL,
      `limiting_mag` double(6,3) DEFAULT NULL,
      `mjd` double DEFAULT NULL,
      `skycell` varchar(20) DEFAULT NULL,
      `dateCreated` datetime DEFAULT current_timestamp(),
      `dateLastModified` datetime DEFAULT current_timestamp() ON UPDATE current_timestamp(),
      `processed` tinyint(1) DEFAULT 0,
      `stacked` tinyint(4) DEFAULT 0,
      PRIMARY KEY (`primaryId`),
      UNIQUE KEY `imageid_skycell` (`imageid`,`skycell`),
      KEY `processed` (`processed`),
      KEY `skycell` (`skycell`),
      KEY `mjd` (`mjd`)
    ) ENGINE=InnoDB AUTO_INCREMENT=0 DEFAULT CHARSET=latin1 COLLATE=latin1_swedish_ci;
    """

    writequery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
    )

    log.debug('completed the ``create_tables_if_not_exist`` function')
    return None


if __name__ == '__main__':
    main()
