library(tidyverse)
library(odbc)
library(DBI)
library(ggthemes)
library(ggplot2)
library(shiny)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q <- "SELECT gene_symbol,COUNT(DISTINCT aliquot_barcode) AS ct
FROM analysis.called_genotypes
WHERE gene_symbol <> 'Unknown' AND variant_classification IN ('In_Frame_Ins','In_Frame_Del','Frame_Shift_Del','Frame_Shift_Ins','Missense_Mutation','Nonstop_Mutation','Nonsense_Mutation')
GROUP BY gene_symbol
HAVING COUNT(DISTINCT aliquot_barcode) > 9
ORDER BY ct DESC"
genes <- dbGetQuery(con, q)$gene_symbol

q <- "SELECT ts.case_barcode, case_sex, s.idh_status, s.codel_status, v.gene_symbol, ts.sample_type, ts.sample_barcode, 
v.chrom, v.start, v.end, v.alt, v.variant_classification, v.variant_type, 
gt.aliquot_barcode, v.hgvs_p, gt.ref_count, gt.alt_count, gt.read_depth, sift, polyphen, gt.called, mf.coverage_adj_mut_freq
FROM analysis.snvs v
FULL JOIN analysis.snv_genotypes gt ON v.chrom = gt.chrom AND v.start = gt.start AND v.end = gt.end AND v.alt = gt.alt
INNER JOIN biospecimen.aliquots ta ON ta.aliquot_barcode = gt.aliquot_barcode
INNER JOIN biospecimen.samples ts ON ts.sample_barcode = ta.sample_barcode
LEFT JOIN analysis.mutation_freq mf ON mf.aliquot_barcode = gt.aliquot_barcode
LEFT JOIN clinical.surgeries s ON s.sample_barcode = ts.sample_barcode
INNER JOIN clinical.cases ca ON ts.case_barcode = ca.case_barcode
WHERE ts.sample_type IN ('TP','R1') 
AND v.variant_classification IN ('In_Frame_Ins','In_Frame_Del','Frame_Shift_Del','Frame_Shift_Ins','Missense_Mutation','Nonstop_Mutation','Nonsense_Mutation')
AND v.gene_symbol = ?"

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
    
    df = qres %>% 
      filter(coverage_adj_mut_freq < 8, complete.cases(idh_status, codel_status)) %>% ## (1) Filter out hypermutator samples
      group_by(case_barcode) %>%
      mutate(idh_codel_grp = ifelse(any(idh_status == "IDHmut") && any(codel_status == "codel"), "IDHmut-codel", 
                                    ifelse(any(idh_status == "IDHmut"), "IDHmut",
                                           ifelse(any(idh_status == "IDHwt"), "IDHwt", NA)))) %>%
      ungroup() %>%
      mutate(var = sprintf("%s:%s-%s_%s", chrom, start, end, alt),
             called = called == "1",
             severity_score = case_when(variant_classification == "Nonsense" ~ 0,
                                        variant_classification %in% c("Frame_Shift_Del","Frame_Shift_Del") ~ 1,
                                        variant_classification %in% c("In_Frame_Del","In_Frame_Ins") ~ 2,
                                        variant_classification == "Missense_Mutation" ~ 3,
                                        variant_classification == "5'Flank" ~ 4)) %>%
      group_by(sample_barcode, var) %>%
      mutate(optimal_variant = order(read_depth, decreasing = T)) %>%
      ungroup() %>%
      filter(optimal_variant == 1) %>% ## (2) For each sample/variant combination, select the variant with the highest read depth, eg. when a variant was profiled across multple sectors or with both WGS and WES
      group_by(case_barcode, var) %>%
      mutate(sufficient_dp = all(ref_count + alt_count > 14)) %>%
      ungroup() %>%
      filter(sufficient_dp) %>% ## (3) For each patient/variant combination filter out variants that do not have >14x coverage across all subsamples
      group_by(case_barcode,gene_symbol) %>%
      mutate(any_called = any(called, na.rm=T),
             num_samples = n_distinct(sample_barcode)) %>%
      ungroup() %>%
      filter(num_samples > 1) %>% ## (4) filter out singletons (unpaired samples)
      group_by(gene_symbol) %>%
      mutate(num_patient = n_distinct(case_barcode),
             gene_symbol_label = sprintf("%s (n=%s)", gene_symbol, num_patient)) %>%
      ungroup() %>%
      filter(any_called) %>% ## (4) a. filter out gene/patient combinations where that gene was not called as mutant in any subsample and b. filter out singletons (unpaired samples)
      arrange(severity_score) %>%
      group_by(case_barcode, gene_symbol) %>%
      mutate(keep = chrom == chrom[which(called)[1]] & 
               start == start[which(called)[1]] & 
               end == end[which(called)[1]] & 
               alt == alt[which(called)[1]]) %>%
      ungroup() %>%
      filter(keep) ## (5) For each patient/gene combination, keep only those specific variants that were called in at least one subsample
    
    ggdat = df %>%
      mutate(vaf = alt_count/(alt_count+ref_count)) %>%
      select(case_barcode, idh_codel_grp, case_sex, gene_symbol_label, chrom, start, end, hgvs_p, variant_type, variant_classification, var, sample_type, vaf, dp = read_depth) %>%
      gather(variable, value, vaf, dp) %>%
      unite(temp, sample_type, variable) %>%
      spread(temp, value) 
    
    ggplot(ggdat, aes(TP_vaf, R1_vaf)) +
      geom_point(aes(color=variant_classification)) + 
      geom_abline(slope=1, alpha=0.2, linetype=2) +
      facet_wrap(~gene_symbol_label) +
      labs(x="Primary", y="First Recurrence", color = "Variant Classification", shape = "Glioma subtype") +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      theme_bw(base_size = 18) +
      theme(axis.text=element_text(size=10))
    
  })
  
  #updateSelectizeInput(session, 'genes', choices = genes, server = TRUE)
  
}

runApp(shinyApp(ui, server), host = "10.7.0.151", port = 2019)