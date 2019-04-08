library(shinydashboard)
library(shiny)
library(DBI)
library(tidyverse)
con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")
tables <- dbGetQuery(con, "SELECT table_schema, table_name, table_type FROM information_schema.tables WHERE table_schema != 'pg_catalog' AND table_schema !='information_schema'")
mviews <- dbGetQuery(con, "SELECT schemaname AS table_schema, matviewname AS table_name FROM pg_matviews")
if(nrow(mviews) > 0)
  tables <- rbind(tables, cbind(mviews, table_type = "Materialized View"))
tables <- tables %>% mutate(table_type = fct_recode(table_type, "View" = "VIEW", "Table" = "BASE TABLE")) %>% arrange(table_name)

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(
    sidebarMenu(
      #lapply(unique(tables$table_schema), function(schema) menuItem(schema, 
      #                                                              lapply(tables$table_name[tables$table_schema==schema],
      #                                                                     function(table) menuSubItem(table, tabName = table, icon = icon("th"))), 
      #                                                              tabName = schema, icon = icon("dashboard")))
      lapply(unique(tables$table_schema), function(schema) menuItem(schema,
        lapply(unique(tables$table_type[tables$table_schema==schema]), function(type) menuItem(type,
          lapply(tables$table_name[tables$table_schema==schema & tables$table_type==type], function(table) menuSubItem(table, 
            tabName = table, icon = icon("th"))), tabName = schema, icon = icon("dashboard")))))
    )
  ),
  dashboardBody(
    do.call(tabItems, 
      lapply(tables$table_name, function(table) {
        tabItem(tabName = table, 
                h2(sprintf("%s > %s", tables$table_schema[tables$table_name==table], table)), DT::dataTableOutput(outputId = table),
                downloadButton(sprintf("d%s", table), "Download"))
      })
    )
  )
)
server <- function(input, output) {
  lapply(tables$table_name, function(table) {
    selectedData <- reactive(dbGetQuery(con, sprintf("SELECT * FROM %s.%s LIMIT %s", tables$table_schema[tables$table_name==table], table, 50000)))
    output[[table]] <- DT::renderDataTable(selectedData())
    output[[sprintf("d%s", table)]] <- downloadHandler(
        filename = function() {
          paste(tables$table_schema[tables$table_name==table],".",table,"_",strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y-%m-%dT%H:%M:%S%z"),".csv",sep="")
        },
        content = function(file) {
          write.csv(selectedData(), file, row.names = FALSE)
        }
      )
  })
}
runApp(shinyApp(ui, server), host = "10.7.0.151", port = 2018)
