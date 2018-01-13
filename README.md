# README #

This is a BAT fit of the Alpha Background.

### Dependencies ###

* jsoncpp : https://github.com/open-source-parsers/jsoncpp
* progressbar : https://github.com/gipert/progressbar

Create shared library for ProgressBar
gcc -c -Wall -Werror -fpic ProgressBar.cc
gcc -shared -o libProgressBar.so ProgressBar.cc

### ENVIRONMENT ###

set JSONCPP_BASE_DIR
set PROGRESS_BAR_DIR

### Contact ###

vonsturm@pd.infn.it
