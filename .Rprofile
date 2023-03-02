source("renv/activate.R")
# Also source user-level R profile
if (file.exists("~/.Rprofile")) {
  source("~/.Rprofile")
}

# Set dwctaxon options
dwctaxon::dct_options(
  # - won't error on duplicated sci names
  check_sci_name = FALSE,
  valid_tax_status = "accepted, synonym, ambiguous synonym, variant",
  skip_missing_cols = TRUE
  )

# Use conflicted
conflicted::conflict_prefer("filter", "dplyr")
