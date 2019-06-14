#Code to make several figures that compare the changes in ccf of shared neoantigens and shared variants when going from initial to recurren tumors
#Adapted from one of Floris's previous exploratory figures
#Last part of this script makes Extended Data Figure 12A and is labelled clearly
#-----------------------------------------------------

library(DBI)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")  
res <- dbGetQuery(con, read_file("/projects/varnf/GLASS/GLASS/sql/neoag_ccf_shared.sql"))

############################################################################################################################################
## Plot ladder plot comparing initial and recurrence CCF
## For all genes w/ colored driver genes
## Faceted by driver gene
############################################################################################################################################

## Munge data
tmp <- res %>% filter(rank==1) %>% 
  group_by(case_barcode)  %>%
  mutate(w = n()) %>%
  ungroup() %>%
  select(case_barcode, w, idh_codel_subtype, gene_symbol, variant_classification_vep, cellular_prevalence_a, cellular_prevalence_b,is_neoag) %>%
  gather(c(cellular_prevalence_a, cellular_prevalence_b), key = "sample_type", value = "ccf") %>%
  mutate(sample_type = factor(sample_type, levels = c("cellular_prevalence_a", "cellular_prevalence_b"), labels = c("P","R")),
         plot_id = paste(case_barcode,gene_symbol),
         gene_label = factor(case_when(gene_symbol == "TP53" ~ "TP53",
                                gene_symbol == "IDH1" ~ "IDH1",
                                gene_symbol == "ATRX" ~ "ATRX",
                                gene_symbol == "PTEN" ~ "PTEN",
                                gene_symbol == "PIK3CA" ~ "PIK3CA",
                                gene_symbol == "PIK3R1" ~ "PIK3R1",
                                gene_symbol == "NF1" ~ "NF1",
                                gene_symbol == "EGFR" ~ "EGFR",
                                gene_symbol == "CIC" ~ "CIC",
                                gene_symbol == "FUBP1" ~ "FUBP1",
                                gene_symbol == "TERT" ~ "TERT",
                                gene_symbol == "RB1" ~ "RB1",
                                TRUE ~ NA_character_)))

## Perform statistical testing
testWilcoxGroup <- function(df) {
  wtest = wilcox.test(df$ccf ~ df$sample_type, paired = TRUE, conf.int = TRUE)
  data.frame(n = nrow(df)/2,
             median_a = median(df$ccf[df$sample_type=="P"]),
             median_b = median(df$ccf[df$sample_type=="R"]),
             statistic = wtest$statistic,
             estimate = wtest$estimate,
             lcl = wtest$conf.int[1],
             ucl = wtest$conf.int[2],
             wilcox_p = wtest$p.value,
             test_str = sprintf("n=%s\n%s",
                                nrow(df),
                                case_when(wtest$p.value < 0.0001 ~ "P<0.0001",
                                          wtest$p.value > 0.05 ~ sprintf("P=%s", format(round(wtest$p.value, 2),scientific=F)),
                                          TRUE ~ sprintf("P=%s", format(round(wtest$p.value, 4),scientific=F)))),
             stringsAsFactors = FALSE)
}

test_case <- tmp %>% 
  group_by(case_barcode,idh_codel_subtype) %>%
  do(testWilcoxGroup(.)) %>%
  ungroup()

test_neoag_case <- tmp %>% 
  group_by(case_barcode,idh_codel_subtype,is_neoag) %>%
  do(testWilcoxGroup(.)) %>%
  ungroup()

test_subtype <- tmp %>% 
  group_by(idh_codel_subtype) %>%
  do(testWilcoxGroup(.)) %>%
  ungroup()

test_neoag_subtype <- tmp %>% 
  group_by(is_neoag,idh_codel_subtype) %>%
  do(testWilcoxGroup(.)) %>%
  ungroup()

test_gene <- tmp %>% 
  filter(gene_symbol %in% tmp$gene_label) %>%
  group_by(gene_symbol) %>% 
  do(testWilcoxGroup(.)) %>%
  ungroup() %>%
  mutate(gene_label = gene_symbol)

test_gene_subtype <- tmp %>% 
  filter(gene_symbol %in% c("TP53","IDH1")) %>%
  group_by(gene_symbol,idh_codel_subtype) %>% 
  do(testWilcoxGroup(.)) %>%
  ungroup()

## Plot for all genes
pdf("/projects/varnf/GLASS/Figures/resubmission/CCF_neoag_ladder.pdf")
g1 <- ggplot(tmp %>% filter(complete.cases(gene_label)), aes(x=sample_type, y=ccf, group=plot_id, color = gene_label)) + 
  geom_point(data = tmp %>% filter(!complete.cases(gene_label))) + 
  geom_line(data = tmp %>% filter(!complete.cases(gene_label))) + 
  geom_point() + 
  geom_line(na.rm = TRUE) + 
  geom_text(data = test_neoag_subtype, aes(x=1.5, y=0.05, label = test_str), group = NA, color = "black") +
  facet_wrap(is_neoag~idh_codel_subtype) +
  scale_color_manual(values = c("#B4464B", "#B47846", "#B4AF46", "#82B446", "#4BB446", "#46B478", "#46B4AF", "#4682B4", "#464BB4", "#7846B4", "#AF46B4", "#B44682"), na.value = "gray75") +
  theme_bw(base_size = 12) +
  labs(x = "Sample Type", y = "Cancer Cell Fraction", color = "Gene Symbol")

g1
dev.off()


## Plot CCF wilcoxon results for each patient
test_neoag_case<- test_neoag_case %>%
  mutate(direction = case_when(wilcox_p < 0.05 & median_b > median_a ~ "CCF increase",
                               wilcox_p < 0.05 & median_b < median_a ~ "CCF decrease",
                               wilcox_p > 0.05 ~ "CCF stable",
                               TRUE ~ "???"))

test_neoag_case_counts <- test_neoag_case %>% count(idh_codel_subtype,is_neoag,direction) %>%
  group_by(idh_codel_subtype,is_neoag) %>% 
  mutate(n_total=sum(nn),
         prop = nn/sum(nn),
         prop_txt = sprintf("%s%%", round(100*nn/sum(nn),1))) %>% 
  ungroup()

pdf("/projects/varnf/GLASS/Figures/resubmission/CCF_neoag_barplot.pdf",width=6,height=6)
g2d <- ggplot(test_neoag_case, aes(x = idh_codel_subtype, fill = direction)) +
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(data = test_neoag_case_counts, aes(label = prop_txt, y = nn + 2), position = position_dodge(width = 1)) +
  facet_grid(is_neoag~.) +
  theme_bw(base_size = 10) +
  labs(x = "Subtype", y = "Number of Patients", fill = "CCF Change") +
  scale_fill_manual(values = c( "#00A3FF", "#5C00FF", "#FFDC00")) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  coord_cartesian(ylim=c(0,65)) +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_blank()) #+
  #facet_wrap(~idh_codel_subtype)

g2d
dev.off()

## Plot bar plots counting clonality

## Munge data
tmp2 <- res %>% filter(rank==1) %>% 
  select(case_barcode, idh_codel_subtype, variant_classification_vep, gene_symbol, cellular_prevalence_a, cellular_prevalence_b,is_neoag) %>%
  mutate(plot_id = paste(case_barcode,gene_symbol),
         gene_label = factor(case_when(gene_symbol == "TP53" ~ "TP53",
                                       gene_symbol == "IDH1" ~ "IDH1",
                                       gene_symbol == "ATRX" ~ "ATRX",
                                       gene_symbol == "PTEN" ~ "PTEN",
                                       gene_symbol == "PIK3CA" ~ "PIK3CA",
                                       gene_symbol == "PIK3R1" ~ "PIK3R1",
                                       gene_symbol == "NF1" ~ "NF1",
                                       gene_symbol == "EGFR" ~ "EGFR",
                                       gene_symbol == "CIC" ~ "CIC",
                                       gene_symbol == "FUBP1" ~ "FUBP1",
                                       gene_symbol == "TERT" ~ "TERT",
                                       gene_symbol == "RB1" ~ "RB1",
                                       TRUE ~ NA_character_)))

#Check to see if frameshift indels are enriched in CCF changes
wilcox_tab <- tmp2[which(tmp2[,"is_neoag"]==1),]
g1 <- nrow(wilcox_tab[which(wilcox_tab[,"cellular_prevalence_b"]<wilcox_tab[,"cellular_prevalence_a"] & 
		   wilcox_tab[,"variant_classification_vep"]=="Frame_Shift_Del" | wilcox_tab[,"variant_classification_vep"]=="Frame_Shift_Ins"),])
g2 <- nrow(wilcox_tab[which(wilcox_tab[,"cellular_prevalence_b"]<wilcox_tab[,"cellular_prevalence_a"] & 
		   wilcox_tab[,"variant_classification_vep"]!="Frame_Shift_Del" & wilcox_tab[,"variant_classification_vep"]!="Frame_Shift_Ins"),])
g3 <- nrow(wilcox_tab[which(wilcox_tab[,"cellular_prevalence_b"]>wilcox_tab[,"cellular_prevalence_a"] & 
		   wilcox_tab[,"variant_classification_vep"]=="Frame_Shift_Del" | wilcox_tab[,"variant_classification_vep"]=="Frame_Shift_Ins"),])
g4 <- nrow(wilcox_tab[which(wilcox_tab[,"cellular_prevalence_b"]>wilcox_tab[,"cellular_prevalence_a"] & 
		   wilcox_tab[,"variant_classification_vep"]!="Frame_Shift_Del" & wilcox_tab[,"variant_classification_vep"]!="Frame_Shift_Ins"),])
fish <- matrix(c(g1,g2,g3,g4),nrow=2,ncol=2)
fisher.test(fish)

## Plot scatter plot facetted per gene
pdf("/projects/varnf/GLASS/Figures/resubmission/CCF_gene_neoag_scatter.pdf",width=6,height=6)
g3 <- tmp2 %>% filter(complete.cases(gene_label)) %>%
  ggplot(aes(x=cellular_prevalence_a, y=cellular_prevalence_b, color = variant_classification_vep)) + 
  geom_point() + 
  facet_wrap(is_neoag~gene_label) +
  theme_bw(base_size = 12) +
  geom_abline(slope=1, alpha=0.2, linetype=2) +
  labs(x="Initial CCF", y="Recurrence CCF", color = "Variant Classification") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
  scale_color_manual(values=c("5'Flank" = brewer.pal(9, "Paired")[9],
                              "Frame_Shift_Del" = brewer.pal(7, "Paired")[1],
                              "Frame_Shift_Ins" = brewer.pal(7, "Paired")[2],
                              "In_Frame_Del" = brewer.pal(7, "Paired")[3],
                              "In_Frame_Ins" = brewer.pal(7, "Paired")[4],
                              "Missense_Mutation" = brewer.pal(7, "Paired")[5],
                              "Nonsense_Mutation" = brewer.pal(7, "Paired")[6],
                              "Splice_Site" = brewer.pal(9, "Paired")[7],
                              "Translation_Start_Site" = brewer.pal(9, "Paired")[8]))

g3
dev.off()

############################################################################################################################################
## Plot bar plots counting clonality
## For all genes w/ colored driver genes
## Faceted by driver gene
############################################################################################################################################

tmp3 <- res %>% filter(rank==1) %>% 
  select(case_barcode, idh_codel_subtype, gene_symbol, clonality_a, clonality_b, is_neoag) %>%
  gather(c(clonality_a, clonality_b), key = "sample_type", value = "clonality") %>%
  mutate(sample_type = factor(sample_type, levels = c("clonality_a", "clonality_b"), labels = c("P","R")),
         plot_id = paste(case_barcode,gene_symbol),
         gene_label = factor(case_when(gene_symbol == "TP53" ~ "TP53",
                                       gene_symbol == "IDH1" ~ "IDH1",
                                       gene_symbol == "ATRX" ~ "ATRX",
                                       gene_symbol == "PTEN" ~ "PTEN",
                                       gene_symbol == "PIK3CA" ~ "PIK3CA",
                                       gene_symbol == "NF1" ~ "NF1",
                                       gene_symbol == "EGFR" ~ "EGFR",
                                       TRUE ~ NA_character_))) %>%
  filter(complete.cases(clonality))

pdf("/projects/varnf/GLASS/Figures/resubmission/CCF_clonality_bar_bad.pdf",width=6,height=6)
g5 <- ggplot(tmp3, aes(x=clonality, fill = sample_type)) + 
  geom_bar(position = "dodge") +
  facet_wrap(is_neoag~idh_codel_subtype, scales = "free_y") + 
  labs(x = "Clonality", y = "Number of Mutations", fill = "Sample Type") +
  theme_bw(base_size = 12) + 
  scale_fill_manual(values = c("#B47846", "#4682B4"))

g5
dev.off()

tmp4 <- res %>% 
  filter(rank==1, clonality_a != "ND", clonality_b != "ND") %>% 
  select(case_barcode, idh_codel_subtype, gene_symbol, clonality_a, clonality_b,is_neoag) %>%
  mutate(clonality = sprintf("%s-%s", clonality_a, clonality_b)) %>%
  filter(complete.cases(clonality_a, clonality_b)) 
  
tmp4[,"is_neoag"] <- as.character(tmp4[,"is_neoag"])
tmp4[which(tmp4[,"is_neoag"]=="0"),"is_neoag"] <- "Non-immunogenic"
tmp4[which(tmp4[,"is_neoag"]=="1"),"is_neoag"] <- "Immunogenic"

tmp4_counts <- tmp4 %>% 
  count(idh_codel_subtype, is_neoag, clonality) %>% 
  group_by(is_neoag,idh_codel_subtype) %>% 
  mutate(n_total=sum(n),
         prop = n/sum(n),
         prop_txt = sprintf("%s%%", round(100*n/sum(n),1))) %>% 
  ungroup() %>%
  add_row(idh_codel_subtype="IDHmut-codel",
  		 is_neoag="Immunogenic",
  		 clonality="S-S",
  		 n = 0,
  		 n_total=76,
  		 prop=0,
  		 prop_txt="0%",
  		 .after=7)

subtypes <- (unique(as.data.frame(tmp4_counts[,"idh_codel_subtype"])))[,1]
pval <- rep(0,3)
for(i in 1:length(subtypes))
{
	tmp_test <- as.data.frame(tmp4_counts[which(tmp4_counts[,"idh_codel_subtype"]==subtypes[i]),])
	cc <- unique(tmp_test[,"clonality"])
	imm <- c("Immunogenic","Non-immunogenic")
	test <- matrix(0,nrow=length(cc),ncol=length(imm))
	rownames(test) <- cc
	colnames(test) <- imm
	test[,1] <- tmp_test[which(tmp_test[,"is_neoag"]==colnames(test)[1]),"n"]
	test[,2] <- tmp_test[which(tmp_test[,"is_neoag"]==colnames(test)[2]),"n"]
	
	pval[i] <- chisq.test(test)$p.value
}

#**********************************************************
#Extended Data Figure 12A: CCF clonality bar
#**********************************************************
gtsize = 6/(14/5)
pdf("/projects/varnf/GLASS/Figures/resubmission/final/EDF_CCF_clonality_bar.pdf", width=6, height = 4)
g6 <- ggplot(tmp4, aes(x=1, fill=clonality)) + 
  geom_bar(position = position_dodge(width = 1,preserve="single")) +
  geom_text(data = tmp4_counts, aes(label = prop_txt, y = n + 0.02 * n_total), position = position_dodge(width = 1),size=gtsize) +
  facet_wrap(is_neoag~idh_codel_subtype, scales = "free") + 
  labs(x = "Subtype", y = "Number of Mutations", fill = "Clonality") +
  theme_bw(base_size = 7)  + 
  scale_fill_manual(values = c("#CA6720", "#2ECA20", "#2083CA", "#BC20CA")) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size=7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_blank()) #+
g6
dev.off()