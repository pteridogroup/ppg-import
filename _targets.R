# Setup ----
source("R/packages.R")
source("R/functions.R")

Sys.setenv(TAR_PROJECT = "main")

# Set number of workers to run in parallel differently for local vs GH actions
targets::tar_option_set(
  controller = crew_controller_local(
    workers = if (nzchar(Sys.getenv("GITHUB_ACTIONS"))) 2 else 10
  )
)

# Note that dwctatxon options are set in .Rprofile

tar_plan(
  # Process World Ferns data into Darwin Core (DWC) format ----
  # - World Ferns data including synonyms
  # - Input path resolved from WORLD_FERNS_INPUT env var; falls back to
  #   newest WorldFerns_ver_*.csv in _targets/user/data_raw/
  tar_file_read(
    wf_with_syn,
    resolve_wf_input(),
    load_raw_wf(path = !!.x)
  ),
  # - Split out only the synonyms
  wf_syns = split_out_syns(wf_with_syn),
  # - Split out only data with no parentage info
  wf_dwc_no_parentage = remove_parentage(wf_with_syn, wf_syns),
  # - Split out numbered ranks
  rank_by_num = get_rank_by_num(wf_with_syn, wf_syns, wf_dwc_no_parentage),
  # - Cleanup into dwctaxon format. Still with original author names from WF.
  #   Data frame at species level, sorted by scientificName
  wf_dwc_auth_orig = clean_wf(
    wf_with_syn,
    wf_syns,
    rank_by_num,
    wf_dwc_no_parentage
  ),

  # Keep original WF authors in main pipeline.
  # IPNI enrichment has moved to _targets_ipni.R.
  wf_dwc = wf_dwc_auth_orig,

  # Count number of taxa at various ranks
  wf_taxa_count = count_taxa_in_wf(wf_with_syn),

  # WF vs PPG comparison ----
  # - Pin PPG data version for reuse in outputs
  ppg_version = "0.0.0.9008",
  # - Load PPG data
  ppg_full = load_ppg(ver = ppg_version),
  # - Build source version metadata for app display
  wf_ppg_data_versions = build_data_versions(ppg_version),
  tar_file(
    wf_ppg_data_versions_csv,
    write_csv_tar(
      wf_ppg_data_versions,
      "_targets/user/results/wf_ppg_data_versions.csv",
      na = ""
    )
  ),
  # - Comparison at genus level (doesn't use author names, so no need for
  # standardization to IPNI)
  wf_ppg_genus_plus = compare_wf_ppg_genus_plus(wf_dwc, ppg_full),
  # - Write genus-level comparison CSV for deployment
  tar_file(
    wf_ppg_genus_plus_csv,
    write_csv_tar(
      wf_ppg_genus_plus,
      "_targets/user/results/wf_ppg_genus_plus.csv",
      na = ""
    )
  )
) |>
  tar_hook_before(
    hook = conflicted::conflict_prefer("filter", "dplyr"),
    names = everything()
  ) |>
  tar_hook_before(
    hook = dwctaxon::dct_options(
      # - won't error on duplicated sci names
      check_sci_name = FALSE,
      valid_tax_status = "accepted, synonym, ambiguous synonym, variant",
      skip_missing_cols = TRUE,
      extra_cols = c(
        "ipniURL",
        "tribe",
        "modified",
        "modifiedBy",
        "modifiedByID"
      )
    ),
    names = everything()
  )
