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

g <- ggplot(tmp, aes(x=idh_codel_subtype, y = mf, color = fraction)) +
  geom_beeswarm() +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = 10^..y.., ymin = 10^..y..), width = .75, linetype = "dashed") +
  scale_y_log10(breaks = c(0.1,1,10,100)) +
  geom_text(data=testres, aes(x=1, y=10^(-2), label = p)) +
  facet_wrap(hypermutator_status~fraction) +
  scale_color_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  theme_bw(base_size = 12) +
  labs(x = "IDH-codel subtype", y = "Mutations per Megabase", color = "Fraction") +
  guides(color = FALSE)

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/EDF2/c-fractionmf.pdf", height = 6, width = 10, useDingbats = FALSE)
plot(g)
dev.off()


## Private MF vs age

lm_stats = function(df){
  m = lm(log10(mf) ~ idh_codel_subtype + time_initial, data = df)
  n = nrow(df)
  
  eq = sprintf("Slope: %s mut/Mb per year", round(10^(coef(m)['time_initial']),2))
  praw = coef(summary(m))[,4]['time_initial']
  p = ifelse(praw < 0.0001, "P<0.0001", ifelse(praw > 0.05, sprintf("P=%s", round(praw, 2)), sprintf("P=%s", round(praw, 4))))
  
  cortxt = sprintf("n=%s\n%s\n%s", n, p, ifelse(praw < 0.05, eq, ''))
  
  if(is.na(p))
    p = "P = 1.00"
  
  return(data.frame(eq,cortxt,n,p))
}

testres <- tmp %>% group_by(hypermutator_status, fraction) %>% do(lm_stats(.))


g <- ggplot(tmp, aes(x=time_initial, y = mf)) +
  geom_point(aes(color = idh_codel_subtype), alpha = 0.3) +
  geom_smooth(method = "lm", aes(color = fraction)) +
  scale_y_log10(breaks = c(0.1,1,10,100)) +
  geom_text(data=testres, aes(x=50, y=10^(-1.8), label = cortxt, color = fraction)) +
  facet_wrap(hypermutator_status~fraction) +
  scale_color_manual(values=c("IDHmut-codel" = "#F8766D", "IDHmut-noncodel" = "#00BA38", "IDHwt" = "#619CFF", "R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  theme_bw(base_size = 12) +
  labs(x = "Age (years)", y = "Mutations per Megabase", color = "Fraction") +
  guides(color = FALSE)

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/EDF2/mf-age.pdf", height = 7, width = 10, useDingbats = FALSE)
plot(g)
dev.off()

## Private MF vs interval

tmp$interval = (tmp$time_recurrence - tmp$time_initial)*12

lm_stats = function(df){
  m = lm(log10(mf) ~ idh_codel_subtype + log10(interval), data = df)
  n = nrow(df)
  
  eq = sprintf("Slope: %s mut/Mb per log10(month)", round((coef(m)['log10(interval)']),4))
  praw = coef(summary(m))[,4]['log10(interval)']
  p = ifelse(praw < 0.0001, "P<0.0001", ifelse(praw > 0.05, sprintf("P=%s", round(praw, 2)), sprintf("P=%s", round(praw, 4))))
  
  cortxt = sprintf("n=%s\n%s", n, p)
  
  if(is.na(p))
    p = "P = 1.00"
  
  return(data.frame(eq,cortxt,n,p))
}

testres <- tmp %>% group_by(hypermutator_status, fraction) %>% do(lm_stats(.))


g <- ggplot(tmp, aes(x=interval, y = mf)) +
  geom_point(aes(color = idh_codel_subtype), alpha = 0.3) +
  geom_smooth(method = "lm", aes(color = fraction)) +
  scale_y_log10(breaks = c(0.1,1,10,100)) +
  scale_x_log10(breaks = c(1,10,100)) +
  geom_text(data=testres, aes(x=10, y=10^(-1.8), label = cortxt, color = fraction)) +
  facet_wrap(hypermutator_status~fraction) +
  scale_color_manual(values=c("IDHmut-codel" = "#F8766D", "IDHmut-noncodel" = "#00BA38", "IDHwt" = "#619CFF", "R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  theme_bw(base_size = 12) +
  labs(x = "Surgical Interval (months)", y = "Mutations per Megabase", color = "Fraction") +
  guides(color = FALSE)

g

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/EDF2/mf-interval.pdf", height = 7, width = 10, useDingbats = FALSE)
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
