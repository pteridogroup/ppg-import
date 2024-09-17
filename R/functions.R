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
  # Delete non-DWC columns, drop vernacular name
  select(
    -photo, -orientation, -author, -trivial_name,
    -conservation_status, -distribution) %>%
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
    taxonRank = taxon,
    namePublishedIn = literature,
    taxonRemarks = remarks
    ) %>%
  # Unescape HTML
  mutate(
    across(
      c(scientificName, namePublishedIn,
        synonyms, taxonRemarks),
      ~unescape_html(.x)
    )) %>%
  # Create unique taxonID
  mutate(taxonID = make_taxon_id(., scientificName, namePublishedIn)) %>%
  # Check that all sci names are unique
  assert(is_uniq, scientificName, taxonID)
}

get_last_infrasp_marker <- function(names) {
  str_match_all(names, " var\\. | f\\. | ssp\\.| subvar\\. | monstr\\. ") %>%
    map_chr(last) %>%
    str_squish()
}

split_out_syns <- function(wf_with_syn) {

  # Split up synonyms from multiple per row to one per row
  wf_syn_split <- wf_with_syn %>%
    filter(!is.na(synonyms)) %>%
    select(number, acceptedNameUsageID = taxonID, synonyms) %>%
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

  # Parse names to help identify rank
  names_parsed <- 
    wf_syn_split %>%
    parse_pterido_names() %>%
    select(scientificName = verbatim, taxon) %>%
    unique() %>%
    assert(is_uniq, scientificName) %>%
    assert(not_na, scientificName)

  # Taxon rank is unknown for synonyms in World Ferns data.
  # Try to ascertain here by number of spaces in taxon name (no author)
  # and presence of 'var.', 'f.', and 'ssp.'
  wf_syn_split %>%
    left_join(
      names_parsed,
      by = join_by(scientificName),
      relationship = "many-to-one") %>%
    assert(not_na, taxon) %>%
    assert(is_uniq, taxonID) %>%
    assert(not_na, taxonID) %>%
    mutate(n_spaces = str_count(taxon, " ")) %>%
    mutate(last_infra = get_last_infrasp_marker(scientificName)) %>%
    mutate(taxonRank = case_when(
      n_spaces == 0 & str_detect(taxon, "phyta$") ~ "division",
      n_spaces == 0 & str_detect(taxon, "phytina$") ~ "subdivision",
      n_spaces == 0 & str_detect(taxon, "opsida$") ~ "class",
      n_spaces == 0 & str_detect(taxon, "idae$") ~ "subclass",
      n_spaces == 0 & str_detect(taxon, "ales$") ~ "order",
      n_spaces == 0 & str_detect(taxon, "ineae$") ~ "suborder",
      n_spaces == 0 & str_detect(taxon, "aceae$") ~ "family",
      n_spaces == 0 & str_detect(taxon, "oideae$") ~ "subfamily",
      n_spaces == 0 & str_detect(taxon, "eae$") ~ "tribe",
      n_spaces == 0 & str_detect(taxon, "inae$") ~ "subtribe",
      last_infra == "var." ~ "variety",
      last_infra == "f." ~ "form",
      last_infra == "ssp." ~ "subspecies",
      last_infra == "subvar." ~ "subvariety",
      last_infra == "monstr." ~ "monstrosity",
      n_spaces == 1 ~ "species",
      # avoid detecting xAsplenosorus x boydstoniae K. S. Walter as genus
      #   parses to Asplenosorus
      n_spaces == 0 & str_detect(scientificName, "^[A-Z][a-z]") ~ "genus",
      .default = NA_character_
    )) %>%
    select(
      taxonRank, number, acceptedNameUsageID, scientificName,
      namePublishedIn, taxonID, taxonomicStatus
    )
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
  
  # Verify that all accepted names have parentNameUsageID except for 'order',
  # which is the highest rank
  wf_dwc %>%
    filter(is.na(parentNameUsageID), taxonomicStatus == "accepted") %>%
    count(taxonRank) %>%
    verify(.$taxonRank == "order", success_fun = success_logical)

  # Same as above, but all synonyms wihout a parentNameUsageID should
  # be synonyms of orders
  wf_dwc %>%
    filter(is.na(parentNameUsageID), taxonomicStatus == "synonym") %>%
    select(scientificName, acceptedNameUsageID) %>%
    left_join(
      select(wf_dwc, acceptedNameUsageID = taxonID, taxonRank)
    ) %>%
    count(taxonRank) %>%
    verify(.$taxonRank == "order", success_fun = success_logical)
  
  # Drop some bad taxa
  # names that are listed as both 'accepted' and 'synonym'
  # MH has verified that synonyms are OK to drop
  wf_dwc %>%
   filter(!(scientificName == "Hypodematiaceae Ching" & taxonomicStatus == "synonym")) %>%
   filter(!(scientificName == "Acrostichum L." & taxonomicStatus == "synonym")) %>%
   filter(!(scientificName == "Lomariopsis Fée" & taxonomicStatus == "synonym")) %>%
   filter(!(scientificName == "Hymenophylloideae Burnett" & taxonomicStatus == "synonym")) %>%
    # Arrange by sciname
    arrange(scientificName)
}

digest_any <- function(...) {
  digest::digest(c(...))
}

make_taxon_id <- function(df, ...) {
  df %>%
    dplyr::select(...) %>%
    purrr::pmap_chr(digest_any)
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
  purrr::map_chr(str, unescape_html_single) %>%
    dplyr::na_if("NA")
}

#' Parse HTML expressions back into plain text for an entire dataframe
#'
#' For example, the ampersand symbol (&) is indicated in HTML like "&amp;"
#' @param df A dataframe
#' @return A dataframe
unescape_html_df <- function (df) {
  mutate(df, across(everything(), unescape_html))
}

#' Filter taxonomic data to only levels at genus or higher
#'
#' @param wf_dwc Dataframe; World Ferns taxonomic data in Darwin Core format
filter_to_genus <- function(wf_dwc) {
  # Remove all names at species and below
  wf_dwc_no_subsp <- wf_dwc %>%
    filter(!taxonRank %in% c(
      "species", "form", "subvariety", "subspecies", "variety"))

  # There are still some rows with NA for taxonRank
  #  - remove any that don't have an acceptedNameUsageID
  wf_dwc_no_subsp_missing_rank <-
    wf_dwc_no_subsp %>%
    filter(!is.na(acceptedNameUsageID)) %>%
    anti_join(wf_dwc_no_subsp, join_by(acceptedNameUsageID == taxonID))

  wf_dwc_no_subsp %>%
    anti_join(wf_dwc_no_subsp_missing_rank, join_by(taxonID)) %>%
    # filter out a few more synonyms with unknown taxonrank
    arrange(scientificName) %>%
    dct_validate()
}

sort_wf_dwc <- function(wf_dwc_unsorted, rank_by_num) {
  # Load ranking of orders for sorting
order_rank <-
  rank_by_num %>%
  filter(taxonRank == "order") %>%
  transmute(
    order = genus_from_sp(scientificName),
    order_rank = str_pad(major_num, 2, "left", "0"))

# Sort data
wf_dwc_unsorted %>%
  arrange_acc_syn() %>%
  left_join(order_rank, by = "order") %>%
  assert(not_na, order_rank) %>%
  mutate(is_order = if_else(taxonRank == "order", 0, 1)) %>%
  mutate(is_family = if_else(taxonRank == "family", 0, 1)) %>%
  mutate(is_subfamily = if_else(taxonRank == "subfamily", 0, 1)) %>%
  mutate(is_tribe = if_else(taxonRank == "tribe", 0, 1)) %>%
  mutate(is_genus = if_else(taxonRank == "genus", 0, 1)) %>%
  mutate(
    higher_tax = paste(
      order_rank, is_order, family, is_family, subfamily, is_subfamily, 
      tribe, is_tribe, genus, is_genus, sep = " | ")) %>%
  arrange(higher_tax) %>%
  mutate(sort_order = 1:nrow(.)) %>%
  select(-contains("is_"))
}

# Write a csv and return the output path.
# for targets pipelines.
write_csv_tar <- function(x, file, ...) {
  write_csv(x = x, file = file, ...)
  file
}


# Extract useful information to dataframe
# TODO: this is used in multiple repos (ppg-import, ppg-voting), so should
# be put into a package
fetch_issues <- function(repo) {

  issues_json <-
    glue::glue("https://api.github.com/repos/{repo}/issues?state=all") |>
    jsonlite::fromJSON()

  # Create initial tibble of issues (may include PRs)
  issues_df <- tibble::tibble(
    number = issues_json$number,
    title = issues_json$title,
    url = issues_json$url,
    created_at = issues_json$created_at,
    user = issues_json$user$login,
    state = issues_json$state,
    body = issues_json$body
  )

  # If any PRs exist, remove them
  if (!is.null(issues_json$draft)) {
    issues_df <-
      issues_df |>
      dplyr::mutate(draft = issues_json$draft) |>
      dplyr::filter(is.na(draft)) |>
      dplyr::select(-draft)
  }

  # Format final data frame
  issues_df |>
    dplyr::mutate(
    url = stringr::str_replace_all(
      url, "https://api.github.com/repos/", "https://github.com/"),
    name = stringr::str_match(body, "Name of taxon[\r|\n]*(.*)[\r|\n]*") |>
             magrittr::extract(, 2),
    rank = stringr::str_match(body, "Rank of taxon[\r|\n]*(\\w+)[\r|\n]*") |>
             magrittr::extract(, 2),
    no_species = stringr::str_match(
      body, "number of species affected[\r|\n]*(.*)") |>
       magrittr::extract(, 2),
    description = stringr::str_match(
      body, "Description of change[\r|\n]*(.*)") |>
        magrittr::extract(, 2)
  ) |>
    dplyr::select(-body)
}

# Check the status of newly approved taxa in World Ferns
check_new_taxa <- function(wf_dwc_gen) {

# Filter github issues to only names that have passed voting
  new_taxa <-
    fetch_issues("pteridogroup/ppg") %>%
    filter(str_detect(title, "\\[PASSED\\]")) %>%
    select(name, rank) %>%
    mutate(
      name = str_replace_all(name, "and", ",")
    ) %>%
    separate_longer_delim(name, delim = ",") %>%
    mutate(name = str_squish(name)) %>%
    rename(taxon = name)

# Left join to World Ferns data to see status of names that passed voting
  wf_dwc_gen_taxa <-
    wf_dwc_gen %>%
    mutate(taxon = gn_parse_tidy(scientificName)$canonicalsimple)
  
  left_join(new_taxa, wf_dwc_gen_taxa) %>%
    rename(new_taxon = taxon, new_rank = rank)
}

parse_pterido_names <- function(wf_dwc) {
  wf_dwc %>%
    pull(scientificName) %>%
    rgnparser::gn_parse_tidy() %>%
    select(
      verbatim,
      taxon = canonicalsimple,
      scientificNameAuthorship = authorship) %>%
    unique() %>%
    # need to manually fix Ericetorum (Jermy) Li Bing Zhang & X. M. Zhou
    mutate(
      taxon = case_when(
         verbatim == "Ericetorum (Jermy) Li Bing Zhang & X. M. Zhou" ~
           "Ericetorum",
         .default = taxon
      )
    ) %>%
    mutate(
      scientificNameAuthorship = case_when(
         verbatim == "Ericetorum (Jermy) Li Bing Zhang & X. M. Zhou" ~
           "(Jermy) Li Bing Zhang & X. M. Zhou",
         .default = scientificNameAuthorship
      )
    )
}

prep_ipni_query <- function(wf_dwc) {

  names_parsed <- parse_pterido_names(wf_dwc)

  names_parsed %>%
    left_join(
      select(wf_dwc, verbatim = "scientificName", taxonRank, taxonID),
      by = "verbatim"
    ) %>%
    filter(
      taxonRank %in% c(
        "family", "subfamily", "genus", "species",
        "subspecies", "variety", "form") | is.na(taxonRank)) %>%
    mutate(ipni_filter = case_when(
      taxonRank == "family" ~ "families",
      taxonRank == "subfamily" ~ "infrafamilies",
      taxonRank == "genus" ~ "genera",
      taxonRank == "species" ~ "species",
      taxonRank %in% c("subspecies", "variety", "form") ~ "infraspecies",
      .default = NA_character_
    )) %>%
    select(taxonID, taxon, ipni_filter) %>%
    unique() %>%
    assertr::assert(assertr::not_na, c(taxonID, taxon))
}

pluck_dat <- function(ipni_res, field) {
    purrr::pluck(ipni_res, "results", 1, field, .default = NA_character_)
  }

search_ipni_single <- function(taxon, ipni_filter = NULL, taxonID) {
  if(is.na(ipni_filter)) {
    ipni_filter <- NULL
  }
  suppressMessages(
    ipni_res <- kewr::search_ipni(
      query = taxon, filters = ipni_filter)
  )
  tibble(
    taxonID = taxonID,
    name = pluck_dat(ipni_res, "name"),
    authors = pluck_dat(ipni_res, "authors"),
    publishingAuthor = pluck_dat(ipni_res, "publishingAuthor"),
    rank = pluck_dat(ipni_res, "rank"),
    publication = pluck_dat(ipni_res, "publication"),
    publicationId = pluck_dat(ipni_res, "publicationId"),
    reference = pluck_dat(ipni_res, "reference"),
    referenceRemarks = pluck_dat(ipni_res, "referenceRemarks"),
    typeName = pluck_dat(ipni_res, "typeName"),
    basionymStr = pluck_dat(ipni_res, "basionymStr"),
    basionymAuthorStr = pluck_dat(ipni_res, "basionymAuthorStr"),
    basionymId = pluck_dat(ipni_res, "basionymId"),
    url = pluck_dat(ipni_res, "url"),
    ipni_id = pluck_dat(ipni_res, "id"),
    fq_id = pluck_dat(ipni_res, "fqId"),
    wfo_id = pluck_dat(ipni_res, "wfoId")
  )
}

search_ipni <- function(ipni_query) {
  purrr::pmap_df(
    dplyr::select(
      ipni_query,
      taxon,
      ipni_filter,
      taxonID),
    search_ipni_single,
    .progress = TRUE)
}

summarize_ipni_results <- function(wf_dwc_auth_orig, ipni_query, ipni_results) {
  wf_dwc_auth_orig |>
    select(taxonID, taxonomicStatus, wf_sci_name = scientificName) |>
    inner_join(ipni_results, join_by(taxonID)) |>
    add_count(name, name = "n_ipni_hits") |>
    mutate(ipni_sci_name = paste(name, authors)) |>
    mutate(
      keep = case_when(
        is.na(name) ~ FALSE,
        wf_sci_name == ipni_sci_name ~ TRUE,
        n_ipni_hits > 1 ~ FALSE,
        n_ipni_hits == 1 ~ TRUE,
        .default = NA
      )
    )
}

# Convert the author names in World Ferns data to those matched from IPNI
make_ppg <- function(wf_dwc_auth_orig, ipni_results_summary) {

  wf_dwc_auth_orig <- wf_dwc_auth_orig %>%
    select(
      -genus, -tribe, -subfamily, -family, -order)

  ipni_results_complete <- ipni_results_summary %>%
    filter(keep) %>%
    transmute(
      taxonID = taxonID,
      scientificName = paste(name, authors),
      namePublishedIn = reference,
      ipniURL = paste0("https://www.ipni.org", url)
    )

  wf_dwc_names_in_ipni <-
    wf_dwc_auth_orig %>%
    select(-scientificName, -namePublishedIn) %>%
    inner_join(ipni_results_complete, by = "taxonID")

  wf_dwc_names_missing_from_ipni <-
    wf_dwc_auth_orig %>%
    anti_join(wf_dwc_names_in_ipni, by = "taxonID")

  wf_dwc_names_in_ipni %>%
    bind_rows(wf_dwc_names_missing_from_ipni) %>%
    # get rid of "subfam." before names in WF data
    mutate(scientificName = str_remove_all(scientificName, "^subfam\\. ")) %>%
    # Update NameUsage cols to match new names from IPNI
    select(-c(acceptedNameUsage, parentNameUsage)) %>%
    dwctaxon::dct_fill_col(
      fill_to = "acceptedNameUsage",
      fill_from = "scientificName",
      match_to = "taxonID",
      match_from = "acceptedNameUsageID"
    ) %>%
    dwctaxon::dct_fill_col(
      fill_to = "parentNameUsage",
      fill_from = "scientificName",
      match_to = "taxonID",
      match_from = "parentNameUsageID"
    ) %>%
    # WFO does not want parent taxa for synonyms, so remove these
    # matches both normal synonyms and ambiguous synonyms
    mutate(
      parentNameUsage = case_when(
        stringr::str_detect(taxonomicStatus, "synonym") ~ NA_character_,
        .default = parentNameUsage
      ),
      parentNameUsageID = case_when(
        stringr::str_detect(taxonomicStatus, "synonym") ~ NA_character_,
        .default = parentNameUsageID
      )
    ) %>%
    select(
      taxonID,
      scientificName,
      taxonRank,
      taxonomicStatus,
      acceptedNameUsage,
      acceptedNameUsageID,
      parentNameUsage,
      parentNameUsageID,
      namePublishedIn,
      taxonRemarks,
      ipniURL
    ) %>%
    # Set final order by sci name, break ties on taxonID
    arrange(scientificName, taxonID) %>%
    dwctaxon::dct_validate(
      extra_cols = "ipniURL"
    )
}

# Helper function for join_higher_taxa()
add_parent_info <- function(df, parent_df, number) {
  parent_name_col <- sym(glue("parent_{number}_name"))
  parent_rank_col <- sym(glue("parent_{number}_rank"))
  highest_parent_name_col <- sym(glue("parent_{number - 1}_name"))

  df %>%
    left_join(
      select(
        parent_df,
        scientificName,
        !!parent_name_col := parentNameUsage,
        !!parent_rank_col := parentNameRank
      ),
      by = setNames("scientificName", as.character(highest_parent_name_col))
    )
}

# Helper function for join_higher_taxa()
pivot_higher_tax <- function(df, names_from, values_from) {
  ranks <- c("order", "family", "subfamily", "tribe", "genus", "species")
  df %>%
    pivot_wider(
      names_from = {{ names_from }},
      values_from = {{ values_from }},
      names_prefix = "."
    ) %>%
    select(-`.NA`)
}

# Helper function for join_higher_taxa()
combine_columns <- function(df, col_name) {
  col_sym <- sym(col_name)
  temp_col_sym <- sym(paste0(".", col_name))

  df %>%
    mutate(!!col_sym := coalesce(!!col_sym, !!temp_col_sym)) %>%
    select(-!!temp_col_sym)
}

# Helper function for join_higher_taxa()
fill_higher_taxon <- function(df, col_name) {
  col_sym <- sym(col_name)

  df %>%
    mutate(
      !!col_sym := case_when(
        taxonRank == col_name ~ scientificName,
        .default = !!col_sym
      )
    )
}

# Join higher taxa (genus, tribe, subfamily, family, order) to PPG dataframe
join_higher_taxa <- function(wf_dwc) {

  # Prepare initial dataframe for joining parent taxa
  # Adds non-DWC 'parentNameRank' column
  wf_dwc_p <-
    wf_dwc %>%
    filter(taxonomicStatus == "accepted") %>%
    select(
      taxonID, scientificName, taxonRank,
      parentNameUsage
    )

  wf_dwc_p <-
    wf_dwc_p %>%
    left_join(
      select(wf_dwc_p, scientificName, parentNameRank = taxonRank),
      join_by(parentNameUsage == scientificName)
    )

  # Make dataframe of taxonID for accepted taxa mapped to their higher taxa  
  accepted_with_higher_taxa <- wf_dwc_p %>%
    # progressively map on higher taxa
    rename(parent_1_name = parentNameUsage, parent_1_rank = parentNameRank) %>%
    add_parent_info(wf_dwc_p, 2) %>%
    add_parent_info(wf_dwc_p, 3) %>%
    add_parent_info(wf_dwc_p, 4) %>%
    add_parent_info(wf_dwc_p, 5) %>%
    add_parent_info(wf_dwc_p, 6) %>%
    add_parent_info(wf_dwc_p, 7) %>%  
    # done when there are no more names to add
    # for species, should be done after adding 6 levels.
    # confirm that 7th doesn't add new information.
    verify(!all(is.na(parent_6_name))) %>%
    verify(all(is.na(parent_7_name))) %>%
    select(-contains("_7_")) %>%
    # start pivoting wider
    pivot_higher_tax(parent_6_rank, parent_6_name) %>%
    rename(order = `.order`) %>%
    pivot_higher_tax(parent_5_rank, parent_5_name) %>%
    rename(family = `.family`) %>%
  # continue to pivot wider and combine duplicated col names
    combine_columns("order") %>%
    pivot_higher_tax(parent_4_rank, parent_4_name) %>%
    rename(subfamily = `.subfamily`) %>%
    combine_columns("order") %>%
    combine_columns("family") %>%
    pivot_higher_tax(parent_3_rank, parent_3_name) %>%
    rename(tribe = `.tribe`) %>%
    combine_columns("order") %>%
    combine_columns("family") %>%
    combine_columns("subfamily") %>%
    pivot_higher_tax(parent_2_rank, parent_2_name) %>%
    rename(genus = `.genus`) %>%
    combine_columns("order") %>%
    combine_columns("family") %>%
    combine_columns("subfamily") %>%
    combine_columns("tribe") %>%
    pivot_higher_tax(parent_1_rank, parent_1_name) %>%
    rename(species = `.species`) %>%
    combine_columns("order") %>%
    combine_columns("family") %>%
    combine_columns("subfamily") %>%
    combine_columns("tribe") %>%
    combine_columns("genus") %>%
    select(-species) %>%
    # Fill in own columns (genus for genus, etc)
    fill_higher_taxon("order") %>%
    fill_higher_taxon("family") %>%
    fill_higher_taxon("subfamily") %>%
    fill_higher_taxon("tribe") %>%
    fill_higher_taxon("genus") %>%
    select(taxonID, genus, tribe, subfamily, family, order)

  syns_with_higher_taxa <-
    wf_dwc %>%
    filter(taxonomicStatus != "accepted") %>%
    left_join(
      accepted_with_higher_taxa,
      join_by(acceptedNameUsageID == taxonID),
      relationship = "many-to-one"
    ) %>%
    select(taxonID, genus, tribe, subfamily, family, order)

  higher_taxa <-
    bind_rows(
      accepted_with_higher_taxa,
      syns_with_higher_taxa
    ) %>%
    select(taxonID, genus, tribe, subfamily, family, order) %>%
    # drop author
    pivot_longer(names_to = "rank", values_to = "taxon", -taxonID) %>%
    mutate(taxon = genus_from_sp(taxon)) %>%
    pivot_wider(names_from = "rank", values_from = "taxon")
    
  wf_dwc %>%
    left_join(higher_taxa, join_by(taxonID), relationship = "one-to-one")
}
