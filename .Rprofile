# renv 0.16.0 + pak 0.9.0 can fail with
# "Cannot parse package: ." when installing packages.
# Use base renv installs for stability.
options(renv.config.pak.enabled = FALSE)

source("renv/activate.R")
# Also source user-level R profile
if (file.exists("~/.Rprofile")) {
  source("~/.Rprofile")
}
