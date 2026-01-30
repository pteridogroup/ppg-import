# Setup ----
source("R/packages.R")
source("R/functions.R")

targets::tar_option_set(
  controller = crew_controller_local(workers = 10)
)

# Note that dwctatxon options are set in .Rprofile

tar_plan(
  # Process World Ferns data into Darwin Core (DWC) format ----
  # - World Ferns data including synonyms
  tar_file_read(
    wf_with_syn,
    "_targets/user/data_raw/WorldFerns_ver_25-12.csv",
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
  wf_dwc_auth_orig_with_dups = clean_wf(
    wf_with_syn,
    wf_syns,
    rank_by_num,
    wf_dwc_no_parentage
  ),
  # - Load manually curated list of duplicate names
  tar_file_read(
    dups_exclude_raw,
    "_targets/user/data_raw/Duplicates_fromPPG_rev_MHa.xlsx",
    readxl::read_xlsx(!!.x)
  ),
  # Remove duplicates
  wf_dwc_auth_orig = remove_dups(
    wf_dwc_auth_orig_with_dups,
    dups_exclude_raw
  ),

  # Lookup author names in IPNI ----
  tar_group_size(
    ipni_query,
    prep_ipni_query(wf_dwc_auth_orig),
    size = 1000
  ),
  tar_target(
    ipni_results,
    search_ipni(ipni_query),
    pattern = map(ipni_query)
  ),
  ipni_results_summary = summarize_ipni_results(
    wf_dwc_auth_orig,
    ipni_query,
    ipni_results
  ),

  # - Replace World Ferns author names and publications with IPNI data
  # when available
  wf_dwc = convert_to_ipni_names(wf_dwc_auth_orig, ipni_results_summary),

  # Count number of taxa at various ranks
  wf_taxa_count = count_taxa_in_wf(wf_with_syn),

  # Produce CSV file ---
  tar_file(
    wf_dwc_csv,
    write_csv_tar(
      wf_dwc,
      "_targets/user/results/wf_dwc.csv",
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
