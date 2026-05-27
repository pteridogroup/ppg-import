source("R/packages.R")
source("R/functions.R")

Sys.setenv(TAR_PROJECT = "ipni")

# Keep IPNI workers modest to reduce transient API failures.
targets::tar_option_set(
  controller = crew_controller_local(
    workers = if (nzchar(Sys.getenv("GITHUB_ACTIONS"))) 2 else 10
  )
)

tar_plan(
  # Process World Ferns data into Darwin Core format with WF authors ----
  tar_file_read(
    wf_with_syn,
    resolve_wf_input(),
    load_raw_wf(path = !!.x)
  ),
  wf_syns = split_out_syns(wf_with_syn),
  wf_dwc_no_parentage = remove_parentage(wf_with_syn, wf_syns),
  rank_by_num = get_rank_by_num(wf_with_syn, wf_syns, wf_dwc_no_parentage),
  wf_dwc_auth_orig = clean_wf(
    wf_with_syn,
    wf_syns,
    rank_by_num,
    wf_dwc_no_parentage
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
  wf_dwc = convert_to_ipni_names(wf_dwc_auth_orig, ipni_results_summary),

  tar_file(
    wf_dwc_ipni_csv,
    write_csv_tar(
      wf_dwc,
      "_targets/user/results/wf_dwc_ipni.csv",
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
