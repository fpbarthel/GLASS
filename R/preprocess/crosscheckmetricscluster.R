tmp = read.delim("/fastscratch/verhaak-lab/GLASS-WG/results/fingerprinting/GLASS-WG.crosscheck_metrics", skip = 6)

d = dist(tmp[,c(3,4,7)])
fit = hclust(tmp[,c(3,4,7)])
plot(as.dendrogram(fit), horiz=T)

x=tmp[,c(1,2,5)] %>% spread(RIGHT_GROUP_VALUE, LOD_SCORE)
rownames(x) = x$LEFT_GROUP_VALUE
x$LEFT_GROUP_VALUE = NULL
x = as.matrix(x)

table(is.na(x))

fit = hclust(dist(x))

fit = hclust(tmp[,c(1,2,5)])

plot(fit)
