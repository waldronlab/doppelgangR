library(knitr)
opts_chunk$set(cache = TRUE)
knit("doc/doppelgangR.Rnw")
system("pdflatex doppelgangR.tex")
