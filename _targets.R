# Setup ----
source("R/packages.R")
source("R/functions.R")

# Note that dwctatxon options are set in .Rprofile

tar_plan(
  # Process World Ferns data into Darwin Core (DWC) format ----
  # - World Ferns data including synonyms
  tar_file_read(
    wf_with_syn,
    "data_raw/001.csv",
    load_raw_wf(path = !!.x)
  ),
  # - Split out only the synonyms
  wf_syns = split_out_syns(wf_with_syn),
  # - Split out only data with no parentage info
  wf_dwc_no_parentage = remove_parentage(wf_with_syn, wf_syns),
  # - Split out numbered ranks
  rank_by_num = get_rank_by_num(wf_with_syn, wf_syns, wf_dwc_no_parentage),
  # - Final clean data frame
  wf_dwc = clean_wf(wf_with_syn, wf_syns, rank_by_num, wf_dwc_no_parentage),
  # Produce report ----
  tar_quarto(
    wf_report,
    "reports/dwc_to_md.Qmd",
    quiet = FALSE,
    packages = c("gluedown", "glue", "tidyverse", "assertr")
  )
)
