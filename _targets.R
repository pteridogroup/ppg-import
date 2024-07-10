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
    "_targets/user/data_raw/ferns_v19-4.csv",
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
    wf_with_syn, wf_syns, rank_by_num, wf_dwc_no_parentage),

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
    wf_dwc_auth_orig, ipni_query, ipni_results),

  # - Replace World Ferns author names and publications with IPNI data
  # when available (final PPG dataframe)
  ppg = make_ppg(wf_dwc_auth_orig, ipni_results_summary),

  # - Data frame at genus and higher
  ppg_gen = filter_to_genus(ppg),

  # Produce report ----

  tar_quarto_rep(
    wf_report,
    "reports/ppg.Qmd",
    execute_params = tibble(
      tax_level = c("species", "genus"),
      output_file = c(
        "_targets/user/results/ppg.md", "_targets/user/results/ppg_gen.md")
    ),
    quiet = FALSE,
    packages = c("gluedown", "glue", "tidyverse", "assertr")
  ),
  # Check the status of newly approved taxa in World Ferns ---
  #
  # The newly approved taxa are appended to the beginning of the data
  # as `new_taxon` and `new_rank`. If the rest of the columns are `NA`,
  # it indicates that these taxa may not be in the WF data.
  name_check = check_new_taxa(ppg_gen),

  # Produce CSV file ---
  # - Data frame at genus and higher, sorted by rank
  tar_file(
    ppg_gen_csv,
    write_csv_tar(
      ppg_gen,
      "_targets/user/results/ppg_gen.csv"
    )
  ),
  tar_file(
    ppg_csv,
    write_csv_tar(
      ppg,
      "_targets/user/results/ppg.csv"
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
        "ipniURL", "tribe", "modified", "modifiedBy", "modifiedByID")
    ),
    names = everything()
  )
