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
import yaml


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
    maps = list_maps_to_be_plotted(dbConn=dbConn, log=log)

    for mmap in maps:
        atlasExps = get_atlas_exposures_covering_map(log=log, dbConn=dbConn, mapId=mmap["mapId"])

        outputFolder = os.path.dirname(mmap["map"])

        # GRAB META
        try:
            yamlFilePath = outputFolder + "/meta.yaml"
            with open(yamlFilePath, 'r') as stream:
                meta = yaml.safe_load(stream)
        except:
            meta = {}

        pointingSide = 0.4
        pointingSide = 5.46

        import matplotlib.path as mpath
        import matplotlib.patches as patches

        atlasPatches = []
        for e in atlasExps:

            raDeg = e['raDeg']
            raDeg = -raDeg + 180
            if raDeg > 180.:
                raDeg -= 360
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
            atlasPatches.append(patches.PathPatch(path))

        from gocart.convert import aitoff
        converter = aitoff(
            log=log,
            mapPath=mmap["map"],
            outputFolder=outputFolder,
            settings=settings,
            plotName="atlas_coverage.png",
            meta=meta,
            patches=atlasPatches
        )
        converter.convert()

    return


def list_maps_to_be_plotted(
        dbConn,
        log):
    """*Generate a list of maps needing maps to be plotted*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger  
    """
    log.debug('starting the ``list_maps_to_be_plotted`` function')

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

    log.debug('completed the ``list_maps_to_be_plotted`` function')
    return maps


def get_atlas_exposures_covering_map(
        log,
        dbConn,
        mapId):
    """*get all of the atlas exposures covering map*

    **Key Arguments:**

    - `dbConn` -- mysql database connection
    - `log` -- logger

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
    log.debug('starting the ``get_atlas_exposures_covering_map`` function')

    from fundamentals.mysql import readquery
    sqlQuery = f"""
        select e.raDeg, e.decDeg from exp_atlas e, alert_pixels_128 p where p.mapId = {mapId} and e.primaryId = p.exp_atlas_id
    """
    atlasExps = readquery(
        log=log,
        sqlQuery=sqlQuery,
        dbConn=dbConn,
        quiet=False
    )

    print(sqlQuery)

    log.debug('completed the ``get_atlas_exposures_covering_map`` function')
    return atlasExps

# use the tab-trigger below for new function
# xt-def-function


if __name__ == '__main__':
    main()
