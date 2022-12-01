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
