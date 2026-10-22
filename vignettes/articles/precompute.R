# conjoint.Rmd is precomputed because irw_fetch() needs a Redivis login, which the
# pkgdown build does not have. To rebuild it, run from the package root:
#   Rscript vignettes/articles/precompute.R
old <- setwd("vignettes/articles")
knitr::knit("conjoint.Rmd.orig", "conjoint.Rmd")
# drop knitr's default "plot of chunk ..." captions
out <- readLines("conjoint.Rmd")
writeLines(out[!grepl('^<p class="caption">plot of chunk', out)], "conjoint.Rmd")
setwd(old)
