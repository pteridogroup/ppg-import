source("R/packages.R")
source("R/functions.R")

Sys.setenv(TAR_PROJECT = "deploy")

tar_plan(
  tar_target(
    app_file,
    "app.R",
    format = "file"
  ),
  tar_target(
    wf_ppg_genus_plus_csv,
    "_targets/user/results/wf_ppg_genus_plus.csv",
    format = "file"
  ),
  tar_target(
    wf_ppg_data_versions_csv,
    "_targets/user/results/wf_ppg_data_versions.csv",
    format = "file"
  ),
  deploy_files = get_shiny_deploy_files(
    app_file = app_file,
    genus_csv = wf_ppg_genus_plus_csv,
    versions_csv = wf_ppg_data_versions_csv
  ),
  deploy_shiny_app = deploy_shiny_app_shinyapps(
    account = Sys.getenv("SHINYAPPS_ACCOUNT", unset = "pteridogroup"),
    app_name = Sys.getenv("SHINYAPPS_APP_NAME", unset = "ppg-wf-explorer"),
    app_files = deploy_files
  )
)
