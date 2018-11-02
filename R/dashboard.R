library(shinydashboard)
library(shiny)
library(DBI)
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
tables <- dbGetQuery(con, "SELECT table_schema, table_name FROM information_schema.tables WHERE table_schema != 'pg_catalog' AND table_schema !='information_schema'")
ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(
    sidebarMenu(
      lapply(unique(tables$table_schema), function(schema) menuItem(schema, 
                                                                    lapply(tables$table_name[tables$table_schema==schema],
                                                                           function(table) menuSubItem(table, tabName = table, icon = icon("th"))), 
                                                                    tabName = schema, icon = icon("dashboard")))
    )
  ),
  dashboardBody(
    do.call(tabItems,
      lapply(tables$table_name, function(table) {
        tabItem(tabName = table,
                h2(table), DT::dataTableOutput(outputId = table))
      })
    )
  )
)
server <- function(input, output) {
  lapply(tables$table_name, function(table) {
    output[[table]] <- DT::renderDataTable(dbReadTable(con, Id(schema=tables$table_schema[tables$table_name==table], table=table)))
  })
}
runApp(shinyApp(ui, server), host = "10.7.0.151", port = 2018)