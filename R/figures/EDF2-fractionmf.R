library(DBI)
library(tidyverse)
library(ggbeeswarm)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")


## Grab mutation paired mutation frequencies (using gold set)
## Seperated by private and shared variants
q <- "SELECT * FROM analysis.tumor_mut_comparison_anno"
res <- dbGetQuery(con,q)

## Plot beeswarm plot of fractionated mutation frequencies
tmp <- res %>% select(case_barcode, idh_codel_subtype, hypermutator_status, time_initial, time_recurrence, mf_private_a, mf_private_b, mf_shared) %>%
  gather(mf_private_a, mf_private_b, mf_shared, key = "fraction", value = "mf") %>%
  mutate(fraction = case_when(fraction == "mf_private_a" ~ "P",
                              fraction == "mf_private_b" ~ "R",
                              fraction == "mf_shared" ~ "S",
                              TRUE ~ NA_character_),
         hypermutator_status = case_when(hypermutator_status == "1" ~ "Hypermutator",
                                         hypermutator_status == "0" ~ "Non-hypermutator",
                                         TRUE ~ NA_character_))

testAOV <- function(df) {
  fit <- aov(log10(df$mf) ~ df$idh_codel_subtype)
  p <- summary(fit)[[1]][["Pr(>F)"]][1]
  p <- case_when(p < 0.0001 ~ "P<0.0001",
                 p > 0.05 ~ sprintf("P=%s", format(round(p, 2),scientific=F)),
                 TRUE ~ sprintf("P=%s", format(round(p, 4),scientific=F)))
  return(data.frame(p=p))
}

testres <- tmp %>% group_by(hypermutator_status, fraction) %>% do(testAOV(.))

g <- ggplot(tmp, aes(x=idh_codel_subtype, y = log10(mf), color = fraction)) +
  geom_beeswarm() +
  geom_text(data=testres, aes(x=1, y=-2, label = p)) +
  facet_wrap(hypermutator_status~fraction) +
  scale_color_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  theme_bw(base_size = 12) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed") +
  labs(x = "IDH-codel subtype", y = "log10 coverage adjusted mutation frequency (mut/Mb)", color = "Fraction") +
  guides(color = FALSE)

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/fractionmf.pdf", height = 6, width = 10, useDingbats = FALSE)
plot(g)
dev.off()


## Plot primary vs recurrent MF scatter plot

tmp2 <- res %>% mutate(received_alk = factor(case_when(received_alk == "1" ~ "Yes",
                                                received_alk == "0" ~ "No",
                                                TRUE ~ "Unknown"), levels = c("Yes", "No", "Unknown")))

g2 <- ggplot(tmp2, aes(x = mf_initial, y = mf_recurrence)) + 
  geom_vline(xintercept = 10, linetype=2, alpha=0.6) +
  geom_hline(yintercept = 10, linetype=2, alpha=0.6) +
  geom_abline(slope=1, linetype=2, alpha=0.6) +
  geom_point(aes(color = received_alk)) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = c(0.1,100), ylim = c(0.1,100)) +
  labs(title = expression(paste("Coverage adjusted mutation frequency ", italic("mutations/covered MB"))),
       x="Initial", y="Recurrence", color = "Alylating Agent Treatment") +
  theme_bw(base_size = 12) +
  theme(axis.text=element_text(size=10)) +
  scale_color_manual(values = c("#4BB446", "#AF48B4", "#B47846"))

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/paired_mf_scatter.pdf", height = 6, width = 8, useDingbats = FALSE)
plot(g2)
dev.off()
