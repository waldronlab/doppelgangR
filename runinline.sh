#!/bin/sh

cd ..
R --vanilla <<RSCRIPT
library(inlinedocs);
package.skeleton.dx("doppelgangR", excludePattern="AllClasses.R")
RSCRIPT

R CMD build doppelgangR

