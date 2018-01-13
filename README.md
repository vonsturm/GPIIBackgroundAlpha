# README

This is a BAT fit of the Alpha Background.

### Dependencies

* jsoncpp : https://github.com/open-source-parsers/jsoncpp
* progressbar : https://github.com/gipert/progressbar

Create shared library for ProgressBar

```
cd $PROGRESS_BAR_DIR
gcc -c -Wall -Werror -fpic ProgressBar.cc
gcc -shared -o libProgressBar.so ProgressBar.cc
```

Set ENVIRONMENT variables

* `JSONCPP_BASE_DIR`
* `PROGRESS_BAR_DIR`


### Contact

vonsturm@pd.infn.it
