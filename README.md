# README

This is a BAT fit of the Alpha Background.

### Usage

`./runGPIIBackgroundAlpha -Z masterconf.json`

### MasterConf

* `verbosity` : 1, # 0: normal 1: debugging
* `reprocess-data` : { `force` : false, `filename` : `./data/data.root` },
* `output-directory` : `./results`,
* `fitoverflow` : false,
* `skipblinded` : { `do` : false, `from` : 2000.0, `to` : 2150.0 },
* `histo` : { `min` : 4850.0, `max` : 5250.0, `binning` : 10.0 },
* `MCMC-fill-histograms` : true,
* `MCMC-write-chain` : true,
* `ROI` : { `from` : 2014.0, `to` : 2064.0 },
* `precision` : `kLow`,
* `update-parameters` : false,
* `parconf` : `./config/all-parconf.json`,
* `runconf` : `./config/all-runconf.json`,
* `detconf` : `./config/all-detlist.json`,
* `envconf` : `./config/mpik-env.json`

See `config/default-masterconf.json` as an example.

### Dependencies

* jsoncpp : https://github.com/open-source-parsers/jsoncpp
* progressbar : https://github.com/gipert/progressbar


### Create shared library for ProgressBar

```
cd $PROGRESS_BAR_DIR
gcc -c -Wall -Werror -fpic ProgressBar.cc
gcc -shared -o libProgressBar.so ProgressBar.cc
```

### Set ENVIRONMENT variables

* `JSONCPP_BASE_DIR`
* `PROGRESS_BAR_DIR`

### Missing Features

* Skip blinding window in fit
* Draw Poisson plots automatically when the fit is finished
* Get posterior distributions for counts in ROI
* Read MC for single detectors and build the MC pdf by scaling with the detector live time

### Contact

vonsturm@pd.infn.it
