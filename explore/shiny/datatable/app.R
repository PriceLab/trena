library(shiny)
library(DT)
tbl <- get(load("~/github/trena/inst/unitTests/tbl.RData"))
app <- shinyApp(
    ui = fluidPage(
        fluidRow(
            column(12,
                   DTOutput('table')
            )
        )
    ),
    server = function(input, output) {
        output$table <- renderDT(tbl,
                                 filter = "top",
                                 options = list(
                                     pageLength = 5
                                 )
        )
    }
)

runApp(app)
