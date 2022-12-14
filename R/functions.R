digest_any <- function(...) {
  digest::digest(c(...))
}

make_taxon_id <- function(df, ...) {
  df %>%
    dplyr::select(...) %>%
    purrr::pmap_chr(digest_any) %>%
    substr(1, 12)
}

end_of_next_val <- function(x) {
  c(x[c(2:length(x))] - 1, Inf)
}

genus_from_sp <- function(x, sep = " ") {
  stringr::str_split(x, sep) %>%
    purrr::map_chr(1)
}

epithet_from_sp_single <- function(x, sep = " ") {
  res <- stringr::str_split(x, sep) %>%
    purrr::map_chr(2)
  if (nchar(res) == 1) {
    res <- stringr::str_split(x, sep) %>%
    purrr::map_chr(3)
  }
  res
}

epithet_from_sp <- function(x, sep = " ") {
  purrr::map2_chr(x, sep, ~epithet_from_sp_single(x = .x, sep = .y))
}

is_uniq_2 <- function(...) {
  is_uniq(..., allow.na = TRUE)
}

not_missing <- function(x) {
  if (length(x) == 0 | length(x) == 1) {
    if (is.null(x)) {
      return(FALSE)
    }
    if (is.na(x)) {
      return(FALSE)
    }
  }
  return(TRUE)
}

# Function to arrange df so accepted names are first, then synonyms
arrange_acc_syn <- function(tax_dat) {
  syns <- tax_dat %>%
  filter(!is.na(acceptedNameUsage)) %>%
  mutate(name_temp = acceptedNameUsage) %>%
  group_by(name_temp) %>%
  arrange(scientificName) %>%
  nest() %>%
  ungroup()

  acc <- tax_dat %>%
    filter(is.na(acceptedNameUsage)) %>%
    mutate(name_temp = scientificName) %>%
    group_by(name_temp) %>%
    nest() %>%
    ungroup()
  
  bind_rows(acc, syns) %>%
    arrange(name_temp) %>%
    unnest(cols = data) %>%
    select(-name_temp)
}
