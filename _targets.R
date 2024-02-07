# Setup ----
source("R/packages.R")
source("R/functions.R")

tar_option_set(
  packages = c("dwctaxon"),
  imports = c("dwctaxon")
)

# Note that dwctatxon options are set in .Rprofile

tar_plan(
  # Process World Ferns data into Darwin Core (DWC) format ----
  # - World Ferns data including synonyms
  tar_file_read(
    wf_with_syn,
    "_targets/user/data_raw/ferns_v18-1.csv",
    load_raw_wf(path = !!.x)
  ),
  # - Split out only the synonyms
  wf_syns = split_out_syns(wf_with_syn),
  # - Split out only data with no parentage info
  wf_dwc_no_parentage = remove_parentage(wf_with_syn, wf_syns),
  # - Split out numbered ranks
  rank_by_num = get_rank_by_num(wf_with_syn, wf_syns, wf_dwc_no_parentage),
  # - Data frame at species level
  wf_dwc_unsorted = clean_wf(
    wf_with_syn, wf_syns, rank_by_num, wf_dwc_no_parentage),
  # - Data frame at species level, sort by rank
  wf_dwc = sort_wf_dwc(wf_dwc_unsorted, rank_by_num),
  # - Data frame at genus and higher, sorted by rank
  wf_dwc_gen = filter_to_genus(wf_dwc_unsorted) %>%
    sort_wf_dwc(rank_by_num),
  # Produce report ----
  tar_quarto_rep(
    wf_report,
    "reports/ppg.Qmd",
    execute_params = tibble(
      tax_level = c("species", "genus"),
      output_file = c(
        "_targets/user/results/ppg_sp.md", "_targets/user/results/ppg.md")
    ),
    quiet = FALSE,
    packages = c("gluedown", "glue", "tidyverse", "assertr")
  ),
  # Check the status of newly approved taxa in World Ferns ---
  #
  # The newly approved taxa are appended to the beginning of the data
  # as `new_taxon` and `new_rank`. If the rest of the columns are `NA`,
  # it indicates that these taxa may not be in the WF data.
  name_check = check_new_taxa(wf_dwc_gen),
  # Lookup author names in IPNI ----
  ipni_query = prep_ipni_query(wf_dwc_gen),
  tar_target(
    ipni_results,
    search_ipni(ipni_query),
    pattern = map(ipni_query)
  ),
  # Produce CSV file ---
  # - v2 of CSV file
  wf_dwc_gen_v2 = make_ppg_v2(wf_dwc_gen, ipni_results),
  # - Data frame at genus and higher, sorted by rank
  tar_file(
    wf_dwc_gen_csv,
    write_csv_tar(
      # Parse HTML expressions back into plain text first
      unescape_html_df(wf_dwc_gen) %>%
        select(-order_rank, -higher_tax, -sort_order),
      "_targets/user/results/ppg.csv"
    )
  ),
  tar_file(
    wf_dwc_csv,
    write_csv_tar(
      # Parse HTML expressions back into plain text first
      unescape_html_df(wf_dwc) %>%
        select(-order_rank, -higher_tax, -sort_order),
      "_targets/user/results/ppg_sp.csv"
    )
  )
)
