# Deploy the EpiTracer web app to shinyapps.io.
#
# One-time setup
# --------------
# install.packages(c("rsconnect", "remotes"))
# # so the server can install EpiTracer from GitHub (it is not on CRAN):
# remotes::install_github("alhafidzhamdan/EpiTracer")
# # connect your shinyapps.io account (Account -> Tokens on shinyapps.io):
# rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
#
# Then:  source("inst/shiny/deploy.R")
#
# rsconnect inspects the library() calls in app.R / helpers.R and reinstalls those
# packages on the server. Because EpiTracer was installed from GitHub above,
# rsconnect records that source and rebuilds it remotely (Bioconductor deps and
# all). The first deploy takes a while.

library(rsconnect)

rsconnect::deployApp(
  appDir  = "inst/shiny",
  appName = "epitracer",
  appTitle = "EpiTracer — plot & call amplicon mechanisms",
  appFiles = c("app.R", "helpers.R",
               "example/DEMO1_CN_segments.tsv", "example/DEMO1_SV_bedpe.tsv"),
  forceUpdate = TRUE,
  launch.browser = TRUE
)
