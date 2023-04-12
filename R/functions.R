#' Load raw data from World Ferns
#'
#' @param path Path to delimited data from World Ferns.
#'
#' @return Dataframe
load_raw_wf <- function(path) {
  readr::read_delim(
    path,
    delim = "|",
    col_types = cols(.default = col_character())
  ) %>%
  janitor::clean_names() %>%
  janitor::remove_empty("cols") %>%
  # Remove duplicated names
  filter(
    !name %in% c(
      "Cyathea lewisii (C. V. Morton &amp; Proctor) Proctor",
      "Anemia candidoi Brade"
    )
  ) %>%
  # Fix some mis-classifications
  mutate(
    number = case_when(
      # Change Pleurosoriopsis from subfam Microsoroideae to 
      # Polypodioideae (039.500x)
      name == "Pleurosoriopsis Fomin" ~ "039.5006",
      TRUE ~ number
    )
  ) %>%
  # Delete non-DWC columns
  select(-photo, -orientation, -author) %>%
  # convert taxonRank to full name
  mutate(taxon = case_when(
    taxon == "O" ~ "order",
    taxon == "F" ~ "family",
    taxon == "SF" ~ "subfamily",
    taxon == "T" ~ "tribe",
    taxon == "G" ~ "genus",
    taxon == "S" ~ "species",
    taxon == "FM" ~ "form",
    taxon == "SS" ~ "subspecies",
    taxon == "V" ~ "variety",
    TRUE ~ NA_character_
  )) %>%
  assert(not_na, taxon) %>%
  # fill in "number" for species and below, to enable mapping to genus
  # ! requires same row order as in raw file
  fill(number, .direction = "down") %>%
  rename(
    scientificName = name,
    vernacularName = trivial_name,
    taxonRank = taxon,
    namePublishedIn = literature,
    taxonRemarks = remarks
    ) %>%
  mutate(taxonID = make_taxon_id(., scientificName, namePublishedIn)) %>%
    # all sci names are unique
  assert(is_uniq, scientificName, taxonID)
}

split_out_syns <- function(wf_with_syn) {
  wf_with_syn %>%
  filter(!is.na(synonyms)) %>%
  select(taxonRank, number, acceptedNameUsageID = taxonID, synonyms) %>%
  separate_rows(synonyms, sep = "=") %>%
  filter(!synonyms == "") %>%
  mutate(synonyms = str_squish(synonyms)) %>%
  separate(
    synonyms,
    c("scientificName", "namePublishedIn"),
    sep = "\\[",
    remove = TRUE,
    extra = "merge",
    fill = "right") %>%
  mutate(
    scientificName = str_squish(scientificName),
    namePublishedIn = str_remove_all(namePublishedIn, "\\]$"),
  ) %>%
  unique() %>%
  mutate(
    taxonID = make_taxon_id(
      ., acceptedNameUsageID, scientificName, namePublishedIn)
  ) %>%
  add_count(scientificName) %>%
  mutate(
    taxonomicStatus = case_when(
      n == 1 ~ "synonym",
      n > 1 ~ "ambiguous synonym"
    )
  ) %>%
  select(-n) %>%
  assert(is_uniq, taxonID)
}

remove_parentage <- function(wf_with_syn, wf_syns) {
  # Combine accepted and synonyms
  # - doesn't incldue parentage mapping yet
  wf_with_syn %>%
      select(-synonyms) %>%
      select(any_of(dct_terms$term)) %>%
      mutate(taxonomicStatus = "accepted") %>%
      bind_rows(wf_syns) %>%
      select(any_of(dct_terms$term)) %>%
      dct_validate()
}
get_rank_by_num <- function(wf_with_syn, wf_syns, wf_dwc_no_parentage) {
  wf_with_syn %>%
  select(number, taxonID) %>%
  left_join(
    select(wf_dwc_no_parentage, taxonID, scientificName, taxonRank),
    by = "taxonID"
  ) %>%
  tidyr::extract(number, "major_num", "(\\d{3})\\.", remove = FALSE) %>%
  mutate(major_num = parse_number(major_num)) %>%
  tidyr::extract(number, "minor_num", "\\.(\\d+)", remove = FALSE) %>%
  mutate(
    minor_num = parse_number(minor_num),
    taxonomicStatus = "accepted"
  )

}

clean_wf <- function(wf_with_syn, wf_syns, rank_by_num, wf_dwc_no_parentage) {
    
  # Split out each high-level taxon and parse the number system
  order_by_num <- filter(rank_by_num, taxonRank == "order") %>%
    arrange(number) %>%
    mutate(
      start = major_num,
      end = end_of_next_val(major_num)
    )
  
  family_by_num <- filter(rank_by_num, taxonRank == "family") %>%
    arrange(number) %>%
    mutate(
      start = minor_num,
      end = minor_num + 999
    )
  
  subfamily_by_num <- filter(rank_by_num, taxonRank == "subfamily") %>%
    arrange(number) %>%
    mutate(
      start = minor_num,
      end = minor_num + 999
    )
  
  tribe_by_num <- filter(rank_by_num, taxonRank == "tribe") %>%
    arrange(number) %>%
    mutate(
      start = minor_num,
      end = minor_num + 99
    )
  
  genus_by_num <- filter(rank_by_num, taxonRank == "genus") %>%
    arrange(number) %>%
    mutate(
      start = minor_num,
      end = minor_num + 999
    )
  
  # Map lower to higher taxa
  # (family to order, subfamily to family, genus to subfamily or family)
  family_mapped <-
  family_by_num %>%
    left_join(
      select(
        order_by_num,
        parentNameUsageID = taxonID,
        parentNameUsage = scientificName,
        start,
        end
      ),
      by = join_by(
        major_num >= start, major_num <= end
      )
    ) %>%
    select(any_of(dct_terms$term)) %>%
    assert(not_na, parentNameUsageID) %>%
    dct_check_taxon_id()
  
  subfamily_mapped <-
  subfamily_by_num %>%
    left_join(
      select(
        family_by_num,
        parentNameUsageID = taxonID,
        parentNameUsage = scientificName,
        major_num
      ),
      by = "major_num"
    ) %>%
    select(any_of(dct_terms$term)) %>%
    assert(not_na, parentNameUsageID) %>%
    dct_check_taxon_id()
  
  # tribe maps to subfamily
  # (only have 5 tribes, all in Microsoroideae)
  tribe_mapped <-
    tribe_by_num %>%
    left_join(
      select(
        subfamily_by_num,
        parentNameUsageID = taxonID,
        parentNameUsage = scientificName,
        major_num,
        start,
        end
      ),
      by = join_by(
        major_num, minor_num >= start, minor_num < end
      )
    ) %>%
    select(any_of(dct_terms$term)) %>%
    assert(not_na, parentNameUsageID) %>%
    dct_check_taxon_id()
  
  genus_mapped_to_tribe <-
  genus_by_num %>%
    left_join(
      select(
        tribe_by_num,
        parentNameUsageID = taxonID,
        parentNameUsage = scientificName,
        major_num,
        start,
        end
      ),
      by = join_by(
        major_num, minor_num >= start, minor_num < end
      )
    ) %>%
    select(any_of(dct_terms$term)) %>%
    filter(!is.na(parentNameUsageID)) %>%
    dct_check_taxon_id()
  
  genus_mapped_to_subfamily <-
    genus_by_num %>%
    # tribe is immediately above genus, so exclude those with tribe
    anti_join(genus_mapped_to_tribe, by = "taxonID") %>%
    left_join(
      select(
        subfamily_by_num,
        parentNameUsageID = taxonID,
        parentNameUsage = scientificName,
        major_num,
        start,
        end
      ),
      by = join_by(
        major_num, minor_num >= start, minor_num < end
      )
    ) %>%
    select(any_of(dct_terms$term)) %>%
    filter(!is.na(parentNameUsageID)) %>%
    dct_check_taxon_id()
  
  genus_mapped_to_family <-
  genus_by_num %>%
    anti_join(genus_mapped_to_tribe, by = "taxonID") %>%
    anti_join(genus_mapped_to_subfamily, by = "taxonID") %>%
    left_join(
      select(
        family_by_num,
        parentNameUsageID = taxonID,
        parentNameUsage = scientificName,
        major_num,
        start,
        end
      ),
      by = join_by(
        major_num, minor_num >= start, minor_num < end
      )
    ) %>%
    select(any_of(dct_terms$term)) %>%
    assert(not_na, parentNameUsageID) %>%
    dct_check_taxon_id()
  
  # Make tibble of accepted species for looking up parent taxon of
  # infrasp taxa
  acc_sp <-
    rank_by_num %>%
    filter(taxonRank == "species") %>%
    mutate(
      genericName = genus_from_sp(scientificName),
      specificEpithet = epithet_from_sp(scientificName)
      ) %>%
    # need to match on genus + epithet, so make sure that combination is unique
    assert_rows(col_concat, is_uniq, genericName, specificEpithet)
  
  # Make autonyms (accepted names only) for infraspecific taxa that lack species
  autonyms_acc <-
  wf_dwc_no_parentage %>%
    filter(taxonomicStatus == "accepted") %>%
    filter(taxonRank %in% c("subspecies", "form", "variety")) %>%
    mutate(
      taxon = gn_parse_tidy(scientificName)$canonicalsimple,
      author = gn_parse_tidy(scientificName)$authorship
      ) %>%
    mutate(
      genericName = genus_from_sp(scientificName),
      specificEpithet = epithet_from_sp(scientificName)
      ) %>%
    anti_join(acc_sp, by = c("genericName", "specificEpithet")) %>%
    select(-genericName, -specificEpithet) %>%
    separate(
      taxon,
      c("genericName", "specificEpithet", "infraspecificEpithet"), sep = " ") %>%
    assert(not_na, genericName, specificEpithet, infraspecificEpithet) %>%
    filter(specificEpithet == infraspecificEpithet) %>%
    left_join(
      select(rank_by_num, number, minor_num, major_num, scientificName),
      by = "scientificName"
    ) %>%
    transmute(
      number,
      minor_num,
      major_num,
      scientificName = paste(genericName, specificEpithet, author),
      namePublishedIn = namePublishedIn,
      taxonRank = "species",
      taxonomicStatus = "accepted",
      taxonRemarks = "autonym"
    ) %>%
    mutate(taxonID = make_taxon_id(., scientificName, namePublishedIn))
  
  # Map species to genus
  # - those that are not autonyms first
  species_mapped_no_autonyms <-
    rank_by_num %>%
    filter(taxonRank == "species") %>%
    left_join(
      select(
        genus_by_num,
        parentNameUsageID = taxonID,
        parentNameUsage = scientificName,
        number,
      ),
      by = join_by(number)
    ) %>%
    select(any_of(dct_terms$term)) %>%
    assert(not_na, parentNameUsageID) %>%
    dct_check_taxon_id()
  
  # - then autonyms
  autonyms_mapped <-
    autonyms_acc %>%
    left_join(
      select(
        genus_by_num,
        parentNameUsageID = taxonID,
        parentNameUsage = scientificName,
        number,
      ),
      by = join_by(number)
    ) %>%
    select(any_of(dct_terms$term)) %>%
    assert(not_na, parentNameUsageID) %>%
    dct_check_taxon_id()
  
  species_mapped <- bind_rows(
    species_mapped_no_autonyms,
    autonyms_mapped
    )%>%
    dct_check_taxon_id()
  
  # Map infrasp to species
  infrasp_mapped <-
  rank_by_num %>%
    filter(taxonRank %in% c("subspecies", "form", "variety")) %>%
    mutate(
      genericName = genus_from_sp(scientificName),
      specificEpithet = epithet_from_sp(scientificName)
      ) %>%
    left_join(
      transmute(
        species_mapped,
        parentNameUsageID = taxonID,
        genericName = genus_from_sp(scientificName),
        specificEpithet = epithet_from_sp(scientificName)
      ),
      by = c("genericName", "specificEpithet")
    ) %>%
    select(any_of(dct_terms$term)) %>%
    assert(not_na, parentNameUsageID) %>%
    dct_check_taxon_id()
  
  # Combine the mapped higher taxa
  parent_mapping <-
    order_by_num %>%
    select(any_of(dct_terms$term)) %>%
    bind_rows(family_mapped) %>%
    bind_rows(subfamily_mapped) %>%
    bind_rows(tribe_mapped) %>%
    bind_rows(genus_mapped_to_subfamily) %>%
    bind_rows(genus_mapped_to_family) %>%
    bind_rows(genus_mapped_to_tribe) %>%
    bind_rows(species_mapped) %>%
    bind_rows(infrasp_mapped) %>%
    dct_validate() %>%
    select(taxonID, parentNameUsageID) %>%
    filter(!is.na(parentNameUsageID))
  
  wf_dwc_no_higher_tax <-
    wf_dwc_no_parentage %>%
    # add autonyms
    bind_rows(select(autonyms_acc, any_of(dct_terms$term))) %>%
    # add parent usage for accepted names
    left_join(
      rename(parent_mapping, parentNameUsageID_1 = parentNameUsageID),
      by = "taxonID") %>%
    # add parent usage for synonyms (where parent is that of the accepted name)
    left_join(
      rename(parent_mapping, parentNameUsageID_2 = parentNameUsageID),
      by = c(acceptedNameUsageID = "taxonID")
    ) %>%
    mutate(
      parentNameUsageID = coalesce(parentNameUsageID_1, parentNameUsageID_2)) %>%
    select(-parentNameUsageID_1, -parentNameUsageID_2) %>%
    dct_validate()
  
  wf_dwc_acc_only_no_higher_tax <-
    wf_dwc_no_higher_tax %>%
    filter(taxonomicStatus == "accepted")
  
  # Add higher level taxa columns
  # things are simple above genus: can join order -> family -> subfamily
  higher_taxa_above_genus <-
  wf_dwc_acc_only_no_higher_tax %>%
    filter(taxonRank == "order") %>%
    select(order = scientificName, order_id = taxonID) %>%
    left_join(
      filter(wf_dwc_acc_only_no_higher_tax, taxonRank == "family") %>%
      select(
        family_id = taxonID,
        order_id = parentNameUsageID,
        family = scientificName),
      by = "order_id",
      multiple = "all") %>%
    assert(is_uniq_2, family) %>%
    left_join(
      filter(wf_dwc_acc_only_no_higher_tax, taxonRank == "subfamily") %>%
      select(
        subfamily_id = taxonID,
        family_id = parentNameUsageID,
        subfamily = scientificName),
      by = "family_id",
      multiple = "all") %>%
    assert(is_uniq_2, subfamily) %>%
    left_join(
      filter(wf_dwc_acc_only_no_higher_tax, taxonRank == "tribe") %>%
      select(
        tribe_id = taxonID,
        subfamily_id = parentNameUsageID,
        tribe = scientificName),
      by = "subfamily_id",
      multiple = "all") %>%
    assert(is_uniq_2, tribe)
  
  # at genus level things get tricky
  # three scenarios
  # family -> genus (no subfamily)
  # family -> subfamily -> genus (no tribe)
  # family -> subfamiy -> tribe -> genus (most nested)
  
  higher_taxa_genus_in_tribe <-
    higher_taxa_above_genus %>%
    left_join(
      filter(wf_dwc_acc_only_no_higher_tax, taxonRank == "genus") %>%
      select(
        genus_id = taxonID,
        tribe_id = parentNameUsageID,
        genus = scientificName),
      by = "tribe_id",
      multiple = "all") %>%
    assert(is_uniq_2, genus) %>%
    filter(!is.na(genus))
  
  higher_taxa_genus_in_subfamily <-
    higher_taxa_above_genus %>%
    filter(is.na(tribe)) %>%
    left_join(
      filter(wf_dwc_acc_only_no_higher_tax, taxonRank == "genus") %>%
      select(
        genus_id = taxonID,
        subfamily_id = parentNameUsageID,
        genus = scientificName),
      by = "subfamily_id",
      multiple = "all") %>%
    assert(is_uniq_2, genus) %>%
    filter(!is.na(genus))
  
  higher_taxa_genus_in_family <-
    higher_taxa_above_genus %>%
    left_join(
      filter(wf_dwc_acc_only_no_higher_tax, taxonRank == "genus") %>%
      anti_join(higher_taxa_genus_in_tribe, by = c(taxonID = "genus_id")) %>%
      anti_join(higher_taxa_genus_in_subfamily, by = c(taxonID = "genus_id")) %>%
      select(
        genus_id = taxonID,
        family_id = parentNameUsageID,
        genus = scientificName),
      by = "family_id",
      multiple = "all") %>%
    assert(is_uniq_2, genus) %>%
    filter(!is.na(genus))
  
  # Combine mapped higher taxa
  higher_taxa <-
    bind_rows(
      higher_taxa_genus_in_tribe,
      higher_taxa_genus_in_subfamily,
      higher_taxa_genus_in_family
    ) %>%
    assert(is_uniq_2, genus) %>%
    select(genus, tribe, subfamily, family, order) %>%
    mutate(across(everything(), ~str_split(.x, " ") %>% purrr::map_chr(1))) %>%
    unique() %>%
    assert(not_na, genus) %>%
    assert(is_uniq, genus)
  
  # Map higher taxa to taxonID of accepted names
  higher_taxa_taxon_id_map_acc <- bind_rows(
    filter(wf_dwc_acc_only_no_higher_tax, taxonRank == "order") %>%
      mutate(order = genus_from_sp(scientificName)),
    filter(wf_dwc_acc_only_no_higher_tax, taxonRank == "family") %>%
      mutate(family = genus_from_sp(scientificName)) %>%
      left_join(unique(select(higher_taxa, family, order)), by = "family") %>%
      assert(not_na, family, order),
    filter(wf_dwc_acc_only_no_higher_tax, taxonRank == "subfamily") %>%
      mutate(subfamily = genus_from_sp(scientificName)) %>%
      left_join(
        unique(
          select(higher_taxa, subfamily, family, order)), by = "subfamily"
        ) %>%
      assert(not_na, subfamily, family, order),
    filter(wf_dwc_acc_only_no_higher_tax, taxonRank == "tribe") %>%
      mutate(tribe = genus_from_sp(scientificName)) %>%
      left_join(
        unique(
          select(higher_taxa, tribe, subfamily, family, order)), by = "tribe") %>%
      assert(not_na, tribe, subfamily, family, order),
    filter(
      wf_dwc_acc_only_no_higher_tax,
      taxonRank %in% c("genus", "species", "subspecies", "form", "variety")) %>%
      mutate(genus = genus_from_sp(scientificName)) %>%
      left_join(higher_taxa, by = "genus") %>%
      assert(not_na, genus, family, order)
  ) %>%
    verify(all(taxonID %in% wf_dwc_acc_only_no_higher_tax$taxonID)) %>%
    left_join(
      select(wf_dwc_acc_only_no_higher_tax, taxonID),
      .,
      by = "taxonID"
    ) %>%
    dct_validate(check_col_names = FALSE) %>%
    select(taxonID, genus, tribe, subfamily, family, order) %>%
    assert(not_na, taxonID) %>%
    assert(is_uniq, taxonID)
  
  # Map higher taxa to taxonID of synonyms
  higher_taxa_taxon_id_map_syns <- wf_dwc_no_parentage %>%
    filter(str_detect(taxonomicStatus, "synonym")) %>%
    select(taxonID, acceptedNameUsageID) %>%
    left_join(
      higher_taxa_taxon_id_map_acc, by = c(acceptedNameUsageID = "taxonID")) %>%
    select(-acceptedNameUsageID) %>%
    assert(not_na, taxonID) %>%
    assert(is_uniq, taxonID)
  
  higher_taxa_taxon_id_map <- bind_rows(
    higher_taxa_taxon_id_map_acc,
    higher_taxa_taxon_id_map_syns
  ) %>%
    assert(not_na, taxonID) %>%
    assert(is_uniq, taxonID)
  
  wf_dwc <-
  wf_dwc_no_higher_tax %>%
    left_join(higher_taxa_taxon_id_map, by = "taxonID") %>%
    assert(not_na, order) %>%
    dct_fill_col(
      fill_to = "acceptedNameUsage",
      fill_from = "scientificName",
      match_to = "taxonID",
      match_from = "acceptedNameUsageID"
    ) %>%
    dct_fill_col(
      fill_to = "parentNameUsage",
      fill_from = "scientificName",
      match_to = "taxonID",
      match_from = "parentNameUsageID"
    ) %>%
    mutate(modified = NA) %>%
    dct_validate(check_col_names = FALSE)
  
  # Verify that all names have parentNameUsageID except for 'order' which is the
  # highest rank
  wf_dwc %>%
    filter(is.na(parentNameUsageID)) %>%
    count(taxonRank) %>%
    verify(.$taxonRank == "order", success_fun = success_logical)
  
  wf_dwc
}


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

#' Parse HTML expressions back into plain text
#'
#' For example, the ampersand symbol (&) is indicated in HTML like "&amp;"
#' @param str A character vector
#' @return A character vector
#' @example
#' unescape_html("You &amp; I")
unescape_html <- function(str) {
  unescape_html_single <- function(str){
    xml2::xml_text(xml2::read_html(paste0("<x>", str, "</x>")))
  }
  purrr::map_chr(str, unescape_html_single)
}

#' Filter taxonomic data to only levels at genus or higher and convert HTML
#' expressions to plain text
#'
#' @param wf_dwc Dataframe; World Ferns taxonomic data in Darwin Core format
filter_to_genus <- function(wf_dwc) {
  wf_dwc %>%
    filter(!taxonRank %in% c("form", "species", "subspecies", "variety")) %>%
    # allow one non-valid column name through: tribe
    dct_validate(check_col_names = FALSE) %>%
    filter(taxonomicStatus == "accepted") %>%
    mutate(across(everything(), unescape_html))
}
