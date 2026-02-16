####
########
################
################################
################################################################
################################################################################################################################
# Package Management
#
# not active part of the package, just a file to manage running / updating the package including package testing, environment
# management, package function tests, etc...

################################################################
# Renv Management
# details from:
# https://rstudio.github.io/renv/articles/packages.html / https://cran.r-project.org/web/packages/renv/vignettes/packages.html
#
# SETUP & initialize renv for this package:
#
# to setup renv in this project - it found the DESCRIPTION file and installed those packages! make sure to run from package root...
renv::init()
#
# not going to be found in the DESCRIPTION as a requirement BUT still needed...
renv::install('devtools')
#
# see where the renv installs are going - since this is a github connected package/project none of these should be in this directory!
.libPaths()
# there are still some files that I manually added to the .gitignore file... see obsidian for more detailed notes...
#
# have this file ALSO be ignored by the R build / development process (not running this will throw a warning when running check())
usethis::use_build_ignore("PACKAGE_MANAGEMENT.R")

################################
# maintenance of renv within this package / project...
renv::install() # install a package...
renv::update() # check for updates to installed packages
renv::status() # see if the lockfile needs to be updated
renv::snapshot() # update the lockfile


################################################################
# package development tools...
library(devtools)

# add a new package to this package...
usethis::use_package("ggplot2") # will add ggplot2 to the Imports field of the DESCRIPTION file...

################################
# automated checks:
devtools::load_all() # load all the packages within this project / package...
devtools::check() # run through a package check ALWAYS DO THIS BEFORE COMITTING CHANGES!
devtools::install() # install THIS package after testing everything...

################################
# function specific tests:
usethis::use_testthat() # setup a tests folder to hold the individual function tests... this enables running test()
usethis::use_test("CTtoRE") # this will open (or generate if none exists) a test file for the CTtoRE function...
devtools::test() # to run through ALL of the created tests...

################################
# readme & help docs
devtools::build_readme() # generate the md file from the Rmd file...
devtools::document() # update the help pages from the roxygen2 headers...
