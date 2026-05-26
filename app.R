library(shiny)
library(jsonlite)
library(DT)

genus_csv_path <- "_targets/user/results/wf_ppg_genus_plus.csv"
versions_csv_path <- "_targets/user/results/wf_ppg_data_versions.csv"

if (!file.exists(genus_csv_path)) {
  stop(
    "Missing CSV: ",
    genus_csv_path,
    "\nRun targets::tar_make(names = wf_ppg_genus_plus_csv) first."
  )
}

dat <- readr::read_csv(genus_csv_path, show_col_types = FALSE)

versions <- if (file.exists(versions_csv_path)) {
  readr::read_csv(versions_csv_path, show_col_types = FALSE)
} else {
  tibble::tibble(source = c("World Ferns", "PPG"), version = NA_character_)
}

version_lookup <- versions |>
  dplyr::mutate(source_key = tolower(trimws(source)))

wf_version <- version_lookup |>
  dplyr::filter(source_key == "world ferns") |>
  dplyr::pull(version)
wf_version <- if (length(wf_version)) wf_version[[1]] else NA_character_

ppg_version <- version_lookup |>
  dplyr::filter(source_key == "ppg") |>
  dplyr::pull(version)
ppg_version <- if (length(ppg_version)) ppg_version[[1]] else NA_character_

version_text <- paste0(
  "World Ferns version: ",
  ifelse(is.na(wf_version) || !nzchar(wf_version), "unknown", wf_version),
  " | PPG version: ",
  ifelse(is.na(ppg_version) || !nzchar(ppg_version), "unknown", ppg_version)
)

hide_redundant_rank_values <- function(df) {
  if (!"rank" %in% names(df)) {
    return(df)
  }

  rank_norm <- tolower(trimws(as.character(df$rank)))
  higher_ranks <- c(
    "tribe",
    "subfamily",
    "family",
    "suborder",
    "order",
    "subclass",
    "class"
  )

  for (rank_name in higher_ranks) {
    rank_cols <- grep(
      paste0("^", rank_name, "_(wf|ppg)$"),
      names(df),
      value = TRUE
    )

    if (length(rank_cols) == 0) {
      next
    }

    same_rank_row <- !is.na(rank_norm) & rank_norm == rank_name
    if (!any(same_rank_row)) {
      next
    }

    for (col_name in rank_cols) {
      df[[col_name]][same_rank_row] <- NA_character_
    }
  }

  df
}

dat <- hide_redundant_rank_values(dat)

col_labels <- names(dat) |>
  stringr::str_replace("_wf$", " (WF)") |>
  stringr::str_replace("_ppg$", " (PPG)") |>
  stringr::str_replace_all("_", " ") |>
  tools::toTitleCase()

same_treatment_to_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  if (is.numeric(x)) {
    return(x != 0)
  }

  x_chr <- tolower(trimws(as.character(x)))
  ifelse(
    x_chr %in% c("true", "t", "1", "yes"),
    TRUE,
    ifelse(x_chr %in% c("false", "f", "0", "no"), FALSE, NA)
  )
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML(
      ".top-controls { margin-bottom: 14px; }\n",
      ".version-note { margin: 2px 0 14px 0; color: #4b5563; }\n",
      ".dt-buttons .dt-button { margin-right: 8px !important; }\n",
      ".dt-buttons { margin-bottom: 10px; }"
    ))
  ),
  titlePanel("WF vs PPG: Genus-level and higher comparison"),
  div(class = "version-note", strong("Data versions: "), version_text),
  fluidRow(
    column(
      width = 12,
      actionButton(
        "toggle_diff",
        "Show All Taxa",
        class = "top-controls btn-primary"
      )
    )
  ),
  fluidRow(
    column(
      width = 12,
      DTOutput("table")
    )
  )
)

server <- function(input, output, session) {
  filter_diff <- reactiveVal(TRUE)

  observeEvent(input$toggle_diff, {
    show_diff <- !filter_diff()
    filter_diff(show_diff)

    updateActionButton(
      session,
      "toggle_diff",
      label = if (show_diff) {
        "Show All Taxa"
      } else {
        "Show Different Treatment Only"
      }
    )
  })

  filtered_dat <- reactive({
    if (!filter_diff()) {
      return(dat)
    }

    if (!"same_treatment" %in% names(dat)) {
      return(dat[0, , drop = FALSE])
    }

    same_treatment <- same_treatment_to_logical(dat$same_treatment)
    dat[!is.na(same_treatment) & !same_treatment, , drop = FALSE]
  })

  output$table <- renderDT(
    filtered_dat(),
    filter = "top",
    rownames = FALSE,
    extensions = "Buttons",
    colnames = col_labels,
    options = list(
      pageLength = 25,
      scrollX = TRUE,
      autoWidth = FALSE,
      dom = "Blfrtip",
      buttons = list(
        list(extend = "csv", text = "Download CSV"),
        list(extend = "excel", text = "Download Excel")
      )
    ),
    class = "stripe hover compact"
  )
}

shinyApp(ui, server)
