library(maftools)
library(tidyverse)
library(iterators)
library(ComplexHeatmap)
counter = icount()
setwd('/fastscratch/verhaak-lab/GLASS-WG/')

genedf = openxlsx::read.xlsx("data/ref/glioma_driver_genes.xlsx")
anno = openxlsx::read.xlsx("data/ref/glass_wg_subtypes.xlsx")

codel = read.delim("data/ref/glass_1p19q_codeletion_status.txt", as.is=T) %>% select(aliquot_id=sample_id,codel_status)
anno = anno %>%
  mutate(sample_type = factor(substr(anno$aliquot_id,14,15), levels = c("TP", sprintf("R%s",1:4), "M1"))) %>%
  filter(sample_type != "NB") %>%
  left_join(codel) %>%
  mutate(aliquot_id = gsub("-", ".", aliquot_id),
         IDHCodel = factor(ifelse(codel_status=='codel', 'IDHmut-codel', ifelse(Sequencing.IDH=='Mutant', 'IDHmut', ifelse(Sequencing.IDH=='WT', 'IDHwt', NA))))) %>%
  select(Tumor_Sample_Barcode=aliquot_id, sample_type, IDHCodel)

anno = anno %>% 
  mutate(Patient = substr(Tumor_Sample_Barcode,0,12)) %>%
  arrange(IDHCodel, Patient, sample_type) %>%
  group_by(Patient) %>%
  mutate(group = factor(nextElem(counter) %% 2, levels=c(0,1), labels=c("0","1"))) %>%
  ungroup() %>%
  select(-Patient)

mafdf = read_tsv('results/mutect2/vcf2maf/GLASS.maf',skip=1)
mafdf2 = mafdf %>%
  mutate(Tumor_Sample_Barcode = gsub("-", ".", Tumor_Sample_Barcode)) %>%
  filter(Tumor_Sample_Barcode %in% anno$Tumor_Sample_Barcode,
         Hugo_Symbol %in% genedf$gene) %>%
  filter((Hugo_Symbol == "TERT" & Variant_Classification == '5\'Flank') | !(Variant_Classification %in% c('3\'Flank','3\'UTR','5\'Flank','5\'UTR','Silent','Splice_Region')))

anno = anno %>% filter(Tumor_Sample_Barcode %in% mafdf2$Tumor_Sample_Barcode)

#cndf = read_tsv('data/ref/maf_tools_input_wgs.txt')
#cndf = cndf %>%
#  rename(Tumor_Sample_Barcode = sample_id, Hugo_Symbol = gene) %>%
#  mutate(Tumor_Sample_Barcode = gsub("-", ".", Tumor_Sample_Barcode)) %>%
#  filter(Tumor_Sample_Barcode %in% anno$Tumor_Sample_Barcode,
#         Hugo_Symbol %in% genedf$gene)

maf = read.maf(mafdf2, removeSilent = FALSE, isTCGA = FALSE)
#ssmaf = subsetMaf(maf, genes = genedf$gene, mafObj = TRUE)

fabcolors = list("sample_type" = RColorBrewer::brewer.pal(n = 6,name = 'Accent'), "IDHCodel" = RColorBrewer::brewer.pal(n = 3,name = 'Set2'), "group" = c("black","white"))
names(fabcolors$sample_type) = levels(anno$sample_type)
names(fabcolors$IDHCodel) = levels(anno$IDHCodel)
names(fabcolors$group) = c("0","1")

varcolors = RColorBrewer::brewer.pal(n = 9, name = 'Set1')
names(varcolors) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del', '5\'Flank')



oncoPrint(mat, column_order = anno$Tumor_Sample_Barcode)

pdf('test3.pdf', width = 12, height = 12)
mat = oncoplot(maf, top=20, annotation = as.data.frame(anno), writeMatrix = TRUE, annotationColor = fabcolors, colors = varcolors, sortByAnnotation = TRUE) 
dev.off()

## changing the sort order ##

function (maf, writeMatrix = FALSE, top = 20, genes = NULL, drawRowBar = TRUE, 
          drawColBar = TRUE, showTumorSampleBarcodes = FALSE, annotation = NULL, 
          annotationColor = NULL, genesToIgnore = NULL, removeNonMutated = TRUE, 
          colors = NULL, fontSize = 10, sortByMutation = FALSE, sortByAnnotation = FALSE) 
{
  set.seed(seed = 1024)
  numMat = maf@numericMatrix
  mat_origin = maf@oncoMatrix
  if (ncol(numMat) < 2) {
    stop("Cannot create oncoplot for single sample. Minimum two sample required ! ")
  }
  if (!is.null(genesToIgnore)) {
    numMat = numMat[!rownames(numMat) %in% genesToIgnore, 
                    ]
    mat_origin = mat_origin[!rownames(mat_origin) %in% genesToIgnore, 
                            ]
  }
  if (!is.null(genes)) {
    if (length(genes[!genes %in% rownames(numMat)]) > 0) {
      message("Following genes are not available in MAF:")
      print(genes[!genes %in% rownames(numMat)])
      message("Ignoring them.")
      genes = genes[genes %in% rownames(numMat)]
    }
    if (length(genes) < 2) {
      stop("Provide at least 2 genes.")
    }
    numMat = numMat[genes, , drop = FALSE]
    numMat = sortByMutation(numMat = numMat, maf = maf)
    mat_origin = mat_origin[rownames(numMat), , drop = FALSE]
    mat_origin = mat_origin[, colnames(numMat), drop = FALSE]
    mat = mat_origin
  }
  else {
    if (sortByAnnotation) {
      if (is.null(annotation)) {
        stop("Missing annotation data. Use argument `annotation` to provide annotations.")
      }
      numMat = sortByAnnotation(numMat, maf, annotation)
      mat_origin = mat_origin[rownames(numMat), , drop = FALSE]
      mat_origin = mat_origin[, colnames(numMat), drop = FALSE]
    }
    else {
      if (sortByMutation) {
        numMat = sortByMutation(numMat = numMat, maf = maf)
        mat_origin = mat_origin[rownames(numMat), , drop = FALSE]
        mat_origin = mat_origin[, colnames(numMat), drop = FALSE]
      }
    }
    if (nrow(mat_origin) < top) {
      mat = mat_origin
    }
    else {
      mat = mat_origin[1:top, , drop = FALSE]
    }
  }
  if (removeNonMutated) {
    numMat = numMat[rownames(mat), , drop = FALSE]
    numMat = numMat[, colnames(mat), drop = FALSE]
    tsb = colnames(numMat)
    tsb.exclude = colnames(numMat[, colSums(numMat) == 0, 
                                  drop = FALSE])
    tsb.include = tsb[!tsb %in% tsb.exclude]
    mat = mat[, tsb.include, drop = FALSE]
  }
  if (writeMatrix) {
    write.table(mat, "onco_matrix.txt", sep = "\t", quote = FALSE)
  }
  mat[mat == ""] = "xxx"
  if (is.null(colors)) {
    col = c(RColorBrewer::brewer.pal(12, name = "Paired"), 
            RColorBrewer::brewer.pal(11, name = "Spectral")[1:3], 
            "black", "violet", "royalblue")
    names(col) = names = c("Nonstop_Mutation", "Frame_Shift_Del", 
                           "IGR", "Missense_Mutation", "Silent", "Nonsense_Mutation", 
                           "RNA", "Splice_Site", "Intron", "Frame_Shift_Ins", 
                           "Nonstop_Mutation", "In_Frame_Del", "ITD", "In_Frame_Ins", 
                           "Translation_Start_Site", "Multi_Hit", "Amp", "Del")
  }
  else {
    col = colors
  }
  bg = "#CCCCCC"
  col = c(col, xxx = bg)
  variant.classes = unique(unlist(as.list(apply(mat_origin, 
                                                2, unique))))
  variant.classes = unique(unlist(strsplit(x = variant.classes, 
                                           split = ";", fixed = TRUE)))
  variant.classes = variant.classes[!variant.classes %in% c("xxx")]
  type_col = structure(col[variant.classes], names = names(col[variant.classes]))
  type_col = type_col[!is.na(type_col)]
  type_name = structure(variant.classes, names = variant.classes)
  variant.classes = variant.classes[!variant.classes %in% c("Amp", 
                                                            "Del")]
  vc.mat = unique(unlist(as.list(apply(mat, 2, unique))))
  vc.mat = unique(unlist(strsplit(x = vc.mat, split = ";", 
                                  fixed = TRUE)))
  vc.mat = vc.mat[!vc.mat %in% c("xxx")]
  vc.type_name = structure(vc.mat, names = vc.mat)
  vc.type_col = structure(col[vc.mat], names = names(col[vc.mat]))
  vc.type_col = vc.type_col[!is.na(vc.type_col)]
  if (!is.null(annotation)) {
    annotation[, 1] = gsub(pattern = "-", replacement = ".", 
                           x = annotation[, 1])
    rownames(annotation) = annotation[, 1]
    annotation = annotation[colnames(mat_origin), ]
    annotation = annotation[complete.cases(annotation), ]
    anno.df = data.frame(row.names = annotation[, 1])
    anno.df = cbind(anno.df, annotation[, 2:ncol(annotation)])
    colnames(anno.df) = colnames(annotation)[2:ncol(annotation)]
    if (sortByMutation || sortByAnnotation) {
      sorted.order = colnames(mat)
      anno.df.sorted = as.data.frame(anno.df[sorted.order, 
                                             ])
      rownames(anno.df.sorted) = sorted.order
      colnames(anno.df.sorted) = colnames(anno.df)
      anno.df = anno.df.sorted
    }
    if (!is.null(annotationColor)) {
      bot.anno = ComplexHeatmap::HeatmapAnnotation(df = anno.df, 
                                                   col = annotationColor)
    }
    else {
      bot.anno = ComplexHeatmap::HeatmapAnnotation(anno.df)
    }
  }
  anno_pct = function(index) {
    n = length(index)
    pct = apply(mat_origin[rev(index), ], 1, function(x) sum(!grepl("^\\s*$", 
                                                                    x))/length(x)) * 100
    pct = paste0(round(pct), "%")
    grid::pushViewport(viewport(xscale = c(0, 1), yscale = c(0.5, 
                                                             n + 0.5)))
    grid::grid.text(pct, x = 1, y = seq_along(index), default.units = "native", 
                    just = "right", gp = grid::gpar(fontsize = fontSize))
    grid::upViewport()
  }
  ha_pct = ComplexHeatmap::HeatmapAnnotation(pct = anno_pct, 
                                             width = grid::grobWidth(grid::textGrob("100%", gp = grid::gpar(fontsize = 10))), 
                                             which = "row")
  anno_row_bar = function(index) {
    n = length(index)
    tb = list()
    for (i in nrow(mat):1) {
      x = mat[i, ]
      x = x[x != ""]
      x = x[x != "xxx"]
      x = unlist(strsplit(x, ";"))
      x = sort(x)
      tb[[i]] = table(x)
    }
    tb = rev(tb)
    max_count = max(sapply(tb, sum))
    grid::pushViewport(grid::viewport(xscale = c(0, max_count * 
                                                   1.1), yscale = c(0.5, n + 0.5)))
    for (i in seq_along(tb)) {
      if (length(tb[[i]])) {
        x = cumsum(tb[[i]])
        grid::grid.rect(x = x, i, width = tb[[i]], height = 0.8, 
                        default.units = "native", just = "right", gp = grid::gpar(col = NA, 
                                                                                  fill = type_col[names(tb[[i]])]))
      }
    }
    breaks = grid::grid.pretty(c(0, max_count))
    grid::grid.xaxis(at = breaks, label = breaks, main = FALSE, 
                     gp = grid::gpar(fontsize = 10))
    grid::upViewport()
  }
  ha_row_bar = ComplexHeatmap::HeatmapAnnotation(row_bar = anno_row_bar, 
                                                 width = grid::unit(4, "cm"), which = "row")
  anno_column_bar = function(index) {
    n = length(index)
    tb = apply(mat_origin[, index], 2, function(x) {
      x = unlist(strsplit(x, ";"))
      x = x[!grepl("^\\s*$", x)]
      x = sort(x)
      table(x)
    })
    max_count = max(sapply(tb, sum))
    grid::pushViewport(grid::viewport(yscale = c(0, max_count * 
                                                   1.1), xscale = c(0.5, n + 0.5)))
    for (i in seq_along(tb)) {
      if (length(tb[[i]])) {
        y = cumsum(tb[[i]])
        grid::grid.rect(i, y, height = tb[[i]], width = 0.8, 
                        default.units = "native", just = "top", gp = grid::gpar(col = NA, 
                                                                                fill = type_col[names(tb[[i]])]))
      }
    }
    breaks = grid::grid.pretty(c(0, max_count))
    grid::grid.yaxis(at = breaks, label = breaks, gp = grid::gpar(fontsize = 10))
    grid::upViewport()
  }
  ha_column_bar = ComplexHeatmap::HeatmapAnnotation(column_bar = anno_column_bar, 
                                                    which = "column")
  add_oncoprint = function(type, x, y, width, height) {
    grid::grid.rect(x, y, width - unit(0.5, "mm"), height - 
                      grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = bg))
    for (i in 1:length(variant.classes)) {
      if (any(type %in% variant.classes[i])) {
        grid::grid.rect(x, y, width - unit(0.5, "mm"), 
                        height - grid::unit(1, "mm"), gp = grid::gpar(col = NA, 
                                                                      fill = type_col[variant.classes[i]]))
      }
      else if (any(type %in% "Amp")) {
        grid::grid.rect(x, y, width - unit(0.5, "mm"), 
                        height - grid::unit(1, "mm"), gp = grid::gpar(col = NA, 
                                                                      fill = bg))
        grid::grid.rect(x, y, width - unit(0.5, "mm"), 
                        height - unit(15, "mm"), gp = grid::gpar(col = NA, 
                                                                 fill = type_col["Amp"]))
      }
      else if (any(type %in% "Del")) {
        grid::grid.rect(x, y, width - unit(0.5, "mm"), 
                        height - grid::unit(1, "mm"), gp = grid::gpar(col = NA, 
                                                                      fill = bg))
        grid::grid.rect(x, y, width - unit(0.5, "mm"), 
                        height - grid::unit(15, "mm"), gp = grid::gpar(col = NA, 
                                                                       fill = type_col["Del"]))
      }
    }
  }
  add_oncoprint2 = function(type, x, y, width, height) {
    for (i in 1:length(variant.classes)) {
      if (any(type %in% variant.classes[i])) {
        grid::grid.rect(x, y, width - unit(0.5, "mm"), 
                        height - grid::unit(1, "mm"), gp = grid::gpar(col = NA, 
                                                                      fill = type_col[variant.classes[i]]))
      }
      else if (any(type %in% "Amp")) {
        grid::grid.rect(x, y, width - unit(0.5, "mm"), 
                        height - unit(15, "mm"), gp = grid::gpar(col = NA, 
                                                                 fill = type_col["Amp"]))
      }
      else if (any(type %in% "Del")) {
        grid::grid.rect(x, y, width - unit(0.5, "mm"), 
                        height - grid::unit(15, "mm"), gp = grid::gpar(col = NA, 
                                                                       fill = type_col["Del"]))
      }
    }
  }
  celFun = function(j, i, x, y, width, height, fill) {
    type = mat[i, j]
    if (type != "xxx") {
      typeList = sort(unlist(strsplit(x = as.character(type), 
                                      split = ";")), decreasing = TRUE)
      if (length(typeList) > 1) {
        for (i in 1:length(typeList)) {
          add_oncoprint2(typeList[i], x, y, width, height)
        }
      }
      else {
        for (i in 1:length(typeList)) {
          add_oncoprint(typeList[i], x, y, width, height)
        }
      }
    }
    else {
      add_oncoprint(type, x, y, width, height)
    }
  }
  if (drawColBar) {
    if (is.null(annotation)) {
      ht = ComplexHeatmap::Heatmap(mat, na_col = bg, rect_gp = grid::gpar(type = "none"), 
                                   cell_fun = celFun, row_names_gp = grid::gpar(fontsize = fontSize), 
                                   show_column_names = showTumorSampleBarcodes, 
                                   show_heatmap_legend = FALSE, top_annotation = ha_column_bar, 
                                   top_annotation_height = grid::unit(2, "cm"))
    }
    else {
      ht = ComplexHeatmap::Heatmap(mat, na_col = bg, rect_gp = grid::gpar(type = "none"), 
                                   cell_fun = celFun, row_names_gp = grid::gpar(fontsize = fontSize), 
                                   show_column_names = showTumorSampleBarcodes, 
                                   show_heatmap_legend = FALSE, top_annotation = ha_column_bar, 
                                   top_annotation_height = grid::unit(2, "cm"), 
                                   bottom_annotation = bot.anno)
    }
  }
  else {
    if (is.null(annotation)) {
      ht = ComplexHeatmap::Heatmap(mat, na_col = bg, rect_gp = grid::gpar(type = "none"), 
                                   cell_fun = celFun, row_names_gp = grid::gpar(fontsize = fontSize), 
                                   show_column_names = showTumorSampleBarcodes, 
                                   show_heatmap_legend = FALSE)
    }
    else {
      ht = ComplexHeatmap::Heatmap(mat, na_col = bg, rect_gp = grid::gpar(type = "none"), 
                                   cell_fun = celFun, row_names_gp = grid::gpar(fontsize = fontSize), 
                                   show_column_names = showTumorSampleBarcodes, 
                                   show_heatmap_legend = FALSE, bottom_annotation = bot.anno)
    }
  }
  ha_pct = ComplexHeatmap::HeatmapAnnotation(pct = anno_pct, 
                                             width = grid::grobWidth(grid::textGrob("100%", gp = grid::gpar(fontsize = 10))), 
                                             which = "row")
  ht_list = ha_pct + ht
  if (drawRowBar) {
    ht_list = ht_list + ha_row_bar
  }
  legend = grid::legendGrob(labels = vc.type_name[names(vc.type_col)], 
                            pch = 15, gp = grid::gpar(col = vc.type_col), nrow = 2)
  suppressWarnings(ComplexHeatmap::draw(ht_list, newpage = FALSE, 
                                        annotation_legend_side = "bottom", annotation_legend_list = list(legend)))
}

## END ##
