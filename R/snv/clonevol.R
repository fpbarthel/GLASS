library(clonevol)
library(tidyverse)

## - founding cluster has to be cluster 1

loci <- read_tsv("results/pyclone/run/GLSS-SF-0018-WXS/tables/loci.tsv") 
locidf <- loci %>% 
  filter(cluster_id != 0) %>%
  mutate(sample_id = gsub("-", "_", substr(sample_id, 14, 18)),
         cluster = ifelse(cluster_id==2,1,ifelse(cluster_id>2,cluster_id,cluster_id + 1)),
         variant_allele_frequency = variant_allele_frequency*100) %>%
  select(-cluster_id) %>%
  mutate(sample_id = gsub("TP", "S1", sample_id)) %>%
  mutate(sample_id = gsub("R1", "S2", sample_id)) %>%
  gather(variable, value, -mutation_id, -sample_id, -cluster) %>%
  unite(temp, sample_id, variable) %>%
  spread(temp, value) %>%
  arrange(cluster)

# shorten vaf column names as they will be
vaf.col.names <- grep('_variant_allele_frequency', colnames(locidf), value=T)
ccf.col.names <- grep('_cellular_prevalence$', colnames(locidf), value=T)

sample.names <- gsub('_variant_allele_frequency', '', vaf.col.names)
locidf[, sample.names] <- locidf[, ccf.col.names]
locidf[, sprintf("%s_VAF", sample.names)] <- locidf[, vaf.col.names]

ccf.col.names <- sample.names
vaf.col.names <- sprintf("%s_VAF", sample.names)

sample.groups <- sample.names
names(sample.groups) <- sample.names

clone.colors <- NULL

pdf('box.pdf', width = 3, height = 3, useDingbats = FALSE, title='')
p1 <- plot.variant.clusters(locidf,
                            cluster.col.name = 'cluster',
                            show.cluster.size = FALSE,
                            cluster.size.text.color = 'blue',
                            vaf.col.names = vaf.col.names,
                            vaf.limits = 70,
                            sample.title.size = 12,
                            base_size = 10,
                            violin = FALSE,
                            box = FALSE,
                            jitter = TRUE,
                            jitter.shape = 1,
                            jitter.color = clone.colors,
                            jitter.size = 1,
                            jitter.alpha = 1,
                            jitter.center.method = 'median',
                            jitter.center.size = 1,
                            jitter.center.color = 'darkgray',
                            jitter.center.display.value = 'none',
                            highlight = 'is.driver',
                            highlight.shape = 1,
                            highlight.color = 'blue',
                            highlight.fill.color = 'green',
                            highlight.note.col.name = 'gene',
                            highlight.note.size = 1,
                            order.by.total.vaf = FALSE)
dev.off()
  
## NB. results plotted to file
plot.pairwise(locidf, col.names = vaf.col.names,
              out.prefix = 'variants.pairwise.plot',
              colors = clone.colors)

pdf('flow.pdf', width=3, height=3, useDingbats=FALSE, title='')
plot.cluster.flow(locidf, vaf.col.names = vaf.col.names,
                  sample.names = sample.names,
                  colors = clone.colors)
dev.off()

y = infer.clonal.models(variants = locidf,
                        cluster.col.name = 'cluster',
                        vaf.col.names = ccf.col.names,
                        ccf.col.names = ccf.col.names,
                        sample.groups = sample.groups,
                        cancer.initiation.model='monoclonal',
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 1000,
                        founding.cluster = 1,
                        cluster.center = 'median',
                        ignore.clusters = c(3),
                        clone.colors = clone.colors,
                        min.cluster.vaf = 0.001,
                        # min probability that CCF(clone) is non-negative
                        sum.p = 0.05,
                        # alpha level in confidence interval estimate for CCF(clone)
                        alpha = 0.05)

#locidf[, sample.names] <- locidf[, sample.names]/100

#y <- transfer.events.to.consensus.trees(y,
#                                        locidf[locidf$cluster==3,],
#                                        cluster.col.name = 'cluster',
#                                        event.col.name = 'mutation_id')


y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')

## plot trees

plot.clonal.models(y,
                   # box plot parameters
                   box.plot = TRUE,
                   # fancy.boxplot = TRUE,
                   # fancy.variant.boxplot.highlight = 'is.driver',
                   # fancy.variant.boxplot.highlight.shape = 21,
                   # fancy.variant.boxplot.highlight.fill.color = 'red',
                   # fancy.variant.boxplot.highlight.color = 'black',
                   # fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   # fancy.variant.boxplot.highlight.note.color = 'blue',
                   # fancy.variant.boxplot.highlight.note.size = 2,
                   # fancy.variant.boxplot.jitter.alpha = 1,
                   # fancy.variant.boxplot.jitter.center.color = 'grey50',
                   # fancy.variant.boxplot.base_size = 12,
                   # fancy.variant.boxplot.plot.margin = 1,
                   # fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1.5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = 'output',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 8,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,2))


pdf('trees.pdf', width = 3, height = 5, useDingbats = FALSE)
plot.all.trees.clone.as.branch(y, branch.width = 0.5,
                               node.size = 1, node.label.size = 0.5)
dev.off()
