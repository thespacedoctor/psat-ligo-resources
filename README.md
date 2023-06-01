# PanSTARRS and ATLAS (PSAT) LVK Script and Resources

**WARNING**. This repo is public so do not add sensitive credentials.

## .monitrc

The `.monitrc` is used by monit on ligo@db3 to periodically run commands and monitor the status of the gocart daemon. It is symlinked to `/home/ligo/.monitrc`

## gocart plugins

All scripts beginning with `gp_` are gocart plugins and symlinked to the `/home/ligo/.config/gocart/plugins/` directory.

* `gp_alerts_to_db.py`: parses each alert's metadata and adds to lvk database

## Scripts

`lvk_atlas_panstarrs_exposures_to_database.py`: imports the CSV files of ATLAS and PanSTARRS exposures found in `/home/ligo/coverage` to the lvk database.
