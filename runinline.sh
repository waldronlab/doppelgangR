#!/bin/sh

cd ..
R --vanilla <<RSCRIPT
library(inlinedocs);
package.skeleton.dx("doppelgangR", excludePattern="AllClasses|nonexports|sn")
RSCRIPT

R CMD build doppelgangR

## Note: After doing ./runinline.sh, please remove the line \alias{doppelgangR} from man/doppelgangR-package.Rd.
