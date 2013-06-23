#!/bin/sh

cd ..
R --vanilla <<RSCRIPT
library(inlinedocs);
package.skeleton.dx("doppelgangR")
RSCRIPT

R CMD build doppelgangR

