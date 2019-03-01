library(tidyverse)
library(DBI)
library(ggthemes)
library(ggplot2)
library(shiny)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q <- "SELECT gene_symbol FROM ref.driver_genes"
genes <- c(dbGetQuery(con, q)$gene_symbol,"TERT")

q <- "
    WITH
/*
Define a list of 'selected tumor pairs': select a single primary/recurrent pair per subject/case after removing any aliquots from the blocklist and sorting samples by surgical interval,
so that pairs with the highest surgical interval are prioritized, and pairs with WGS are prioritized over WXS when both are available
*/
selected_tumor_pairs AS
(
  SELECT
  tumor_pair_barcode,
  row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
  FROM analysis.tumor_pairs ps
  LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
  LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
  LEFT JOIN analysis.mutation_freq mf1 ON mf1.aliquot_barcode = ps.tumor_barcode_a  -- join with mutation freq to remove hypermutators
  LEFT JOIN analysis.mutation_freq mf2 ON mf2.aliquot_barcode = ps.tumor_barcode_b
  WHERE
  comparison_type = 'longitudinal' AND
  sample_type_b <> 'M1' AND 													-- exclude metastatic samples here because this is outside the scope of our study
  --b1.cnv_exclusion = 'allow' AND b2.cnv_exclusion = 'allow' AND
  b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' AND
  b1.clinical_exclusion <> 'block' AND b2.clinical_exclusion <> 'block' AND
  mf1.coverage_adj_mut_freq < 10 AND mf2.coverage_adj_mut_freq < 10			-- filter hypermutators
),
  /*
  Aggregate counts over tumor pairs and genes
  Performs an INNER JOIN with selected_tumor_pairs from above so that only those pairs are included in the analysis.
  Restrict to events with coverage >= 15 in both A and B
  Variants for each tumor pair/gene combination are ordered according to variant_classification_priority (see new table analysis.variant_classifications) and whether or not the mutation was called in a/b and finally based on read_depth
  Note that I am adding `mutect2_call_a` to `mutect2_call_b` (true = 1, false = 0) as to avoid prioritizing mutations in either A or B over the other
  The row_number() function asigns a row number to each row within each group of gene_symbol and case_barcode, after ordering by the given parameters
  */
  variants_by_case_and_gene AS
  (
  SELECT
  gtc.gene_symbol,
  gtc.case_barcode,
  gtc.variant_classification,
  alt_count_a::decimal / (alt_count_a + ref_count_a) AS vaf_a,
  alt_count_b::decimal / (alt_count_b + ref_count_b) AS vaf_b,
  row_number() OVER (PARTITION BY gtc.gene_symbol, gtc.case_barcode ORDER BY vc.variant_classification_priority, mutect2_call_a::integer + mutect2_call_b::integer DESC, read_depth_a + read_depth_b DESC) AS priority
  FROM analysis.master_genotype_comparison gtc
  INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode AND stp.priority = 1
  LEFT JOIN analysis.variant_classifications vc ON gtc.variant_classification = vc.variant_classification
  WHERE
  gene_symbol = ? AND
  (mutect2_call_a OR mutect2_call_b) AND
  (alt_count_a + ref_count_a) >= 5 AND (alt_count_b + ref_count_b) >= 5
  AND (variant_classification_priority IS NOT NULL OR -- this removes any variant types we don't care about, eg. Silent and Intronic mutations, see the analysis.variant_classifications table for more details
  (gene_symbol = 'TERT' AND gtc.variant_classification = '5''Flank'))
  )
  SELECT *
  FROM variants_by_case_and_gene vg
  WHERE priority = 1"

vaf_res <- dbSendQuery(con, q)
#vaf_bind <- dbBind(vaf_res, list("EGFR"))
#qres <- dbFetch(vaf_bind)

ui <- fluidPage(
  titlePanel("GLASS VAF comparison"),
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      selectizeInput('genes', label = 'Genes', choices = genes, selected = "TP53")
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      plotOutput(outputId = "distPlot")
    )
  )
)

server <- function(input, output, session) {
  
  output$distPlot <- renderPlot({
    vaf_bind <- dbBind(vaf_res, list(input$genes))
    qres <- dbFetch(vaf_bind)
    
    ggplot(qres, aes(vaf_a, vaf_b)) +
      geom_point(aes(color=variant_classification)) + 
      geom_abline(slope=1, alpha=0.2, linetype=2) +
      labs(x="Primary", y="Recurrence", color = "Variant Classification") +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      theme_bw(base_size = 18) +
      theme(axis.text=element_text(size=10))
    
  })
  
  #updateSelectizeInput(session, 'genes', choices = genes, server = TRUE)
  
}

runApp(shinyApp(ui, server), host = "10.7.0.151", port = 2019)