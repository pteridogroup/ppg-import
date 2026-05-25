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

ui <- fluidPage(
  titlePanel("WF vs PPG: Genus-level Comparison"),
  fluidRow(
    column(
      width = 12,
      DTOutput("table")
    )
  )
)

server <- function(input, output, session) {
  output$table <- renderDT(
    dat,
    filter = "top",
    rownames = FALSE,
    extensions = "Buttons",
    colnames = col_labels,
    options = list(
      pageLength = 25,
      scrollX = TRUE,
      autoWidth = FALSE,
      dom = "Blfrtip",
      buttons = c("csv", "excel")
    ),
    class = "stripe hover compact"
  )
}

shinyApp(ui, server)
