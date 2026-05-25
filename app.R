library(shiny)
library(DT)

genus_csv_path <- "_targets/user/results/wf_ppg_genus_plus.csv"

if (!file.exists(genus_csv_path)) {
  stop(
    "Missing CSV: ",
    genus_csv_path,
    "\nRun targets::tar_make(names = wf_ppg_genus_plus_csv) first."
  )
}

dat <- readr::read_csv(genus_csv_path, show_col_types = FALSE)

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
      ".dt-buttons .dt-button { margin-right: 8px !important; }\n",
      ".dt-buttons { margin-bottom: 10px; }"
    ))
  ),
  titlePanel("WF vs PPG: Genus-level Comparison"),
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
