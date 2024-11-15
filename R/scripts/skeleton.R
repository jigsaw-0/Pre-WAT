# set working space
setwd("/Users/jigsaw-0/Workspace/KMPC/Pre-WAT/R/")

# package handling
packages <- c("tidyverse", "reshape2", "readxl", "writexl",
              "AnnotationHub", "ensembldb", "biomaRt", "fst", "msigdbr",
              "tximport", "DESeq2", "apeglm", "ashr",
              "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db",
              "GSVA", "decoupleR", "OmnipathR", "clusterProfiler",
              "EnhancedVolcano", "circlize", "ggrepel", "ComplexHeatmap")

# install pacakges
BiocManager::install(packages)

# load packages
invisible(lapply(packages, library, character.only = T))



#@@@
# 0. save biomart query
#@@@
CacheBiomaRt()



#@@@
# 1. create DESeq object using salmon output
#@@@

# set base file path 
salmon_output_root <- "data/salmon_output"

# list up sample names
sample_names.orig <- list.files(salmon_output_root)

# change name order (no need to make it as factor)
sample_names.orig <- c(sample_names.orig[27:28],
                       sample_names.orig[23:24],
                       sample_names.orig[13:22],
                       sample_names.orig[25:26],
                       sample_names.orig[9:12],
                       sample_names.orig[1:8])

salmon_output_path <- file.path(salmon_output_root, sample_names.orig, "quant.sf")
names(salmon_output_path) <- sample_names.orig

# prepare tx2gene
ah <- AnnotationHub::AnnotationHub()
ahDb <- AnnotationHub::query(x = ah, pattern = c("ensDb", "Mus musculus"))

ah_latest_rec <- ahDb %>%
    mcols() %>%  # retrieve dataframe containing metadata columns
    rownames() %>%
    tail(n = 1)

ensDb <- ahDb[[ah_latest_rec]]

txData <- ensembldb::transcripts(x = ensDb, return.type = "DataFrame")
tx2gene <- txData[, c("tx_id", "gene_id")]

# create txi object
txi <- tximport::tximport(files = salmon_output_path, 
                          type = "salmon", 
                          tx2gene = tx2gene, 
                          ignoreTxVersion = T)

# prepare metadata
coldata <- data.frame(Sample = sample_names.orig)
coldata <- coldata %>%
    mutate(Group = c(rep("IM_B6", 2), rep("B6_6", 2), rep("B6_20", 2),
                     rep("B6_27", 2), rep("B6_28", 2), rep("B6_29", 2),
                     rep("B6_34", 2), rep("IM_129", 2), rep("129_3", 2),
                     rep("129_9", 2), rep("129_10", 2), rep("129_16", 2),
                     rep("129_22", 2), rep("129_23", 2)),
           MouseType = case_when(grepl("B6", Sample) ~ "B6",
                                 grepl("129", Sample) ~ "129")) %>%
    tibble::column_to_rownames(var = "Sample")

coldata$Group <- factor(coldata$Group, 
                        levels = c("IM_B6", "B6_6", "B6_20", "B6_27", "B6_28", "B6_29", "B6_34",
                                   "IM_129", "129_3", "129_9", "129_10", "129_16", "129_22", "129_23"))
coldata$MouseType <- factor(coldata$MouseType,
                            levels = c("B6", "129"))



#@@@
# 2. basic sample QC & EDA
#@@@

# create DESeq object
dds <- DESeq2::DESeqDataSetFromTximport(txi = txi, 
                                        colData = coldata, 
                                        design = ~Group)

dds <- dds[rowSums(BiocGenerics::counts(dds) >= 10) > 3, ]

ddsX <- DESeq2::DESeq(dds)

vst <- DESeq2::vst(ddsX, blind = F)

group_colors <- c(
    "IM_B6" = "#5e5e5e",
    "B6_6" = "#FFE4C4",
    "B6_20" = "#1E90FF",
    "B6_27" = "purple",
    "B6_28" = "magenta",
    "B6_29" = "#00CD00",
    "B6_34" = "#FF6347",
    "IM_129" = "#5e5e5e",
    "129_3" = "#FFE4C4",
    "129_9" = "purple",
    "129_10" = "#FF6347",
    "129_16" = "magenta",
    "129_22" = "#00CD00",
    "129_23" = "#1E90FF"
)

group_colors2 <- c("B6" = "#FFA07A", "129" = "#A4D3EE")

group_colors.transparent <- sapply(group_colors,
                                   function(x) grDevices::adjustcolor(x, alpha.f = .4))

group_shape <- c("B6" = 21, "129" = 24)


# check RLE plot
# Not working yet
#RLEplot(BiocGenerics::counts(ddsX, normalized = F), 
#        raw.mat = T,
#        group_colors = group_colors.transparent, 
#        box.only = T)
#
#RLEplot(vst, group_colors = group_colors.transparent, box.only = T)


vst.pca <- PCAplot(vst, 
                   grp_col = group_colors,
                   ext_grp = "MouseType",
                   grp_shape = group_shape,
                   pt_sz = 8, lab_sz = 4, glob_txt_sz = 16)


# heatmap annotation for mouse type
vst.col_annot <- ComplexHeatmap::HeatmapAnnotation(df = coldata[, "MouseType", drop = F],
                                                   col = list(MouseType = group_colors2),
                                                   border = T,
                                                   show_annotation_name = F,
                                                   show_legend = T)

vst.row_annot <- ComplexHeatmap::HeatmapAnnotation(df = coldata[, "MouseType", drop = F],
                                                   col = list(MouseType = group_colors2),
                                                   which = "row",
                                                   border = T,
                                                   show_annotation_name = F,
                                                   show_legend = F)

vst.cor <- CorHeatmap(vst, 
                      nr_col_split = 2,
                      nr_row_split = 2,
                      top_annot = vst.col_annot,
                      right_annot = vst.row_annot)

ComplexHeatmap::draw(vst.cor)


# heatmap annotation for groups within each mouse type
vst.col_annot2 <- ComplexHeatmap::HeatmapAnnotation(df = coldata[, "Group", drop = F],
                                                    col = list(Group = group_colors),
                                                    border = T,
                                                    show_annotation_name = F,
                                                    show_legend = T)

vst.row_annot2 <- ComplexHeatmap::HeatmapAnnotation(df = coldata[, "Group", drop = F],
                                                    col = list(Group = group_colors),
                                                    which = "row",
                                                    border = T,
                                                    show_annotation_name = F,
                                                    show_legend = F)

vst.cor.B6 <- CorHeatmap(vst[, coldata$MouseType == "B6"],
                         nr_col_split = 7,
                         nr_row_split = 7,
                         top_annot = vst.col_annot2[coldata$MouseType == "B6", ],
                         right_annot = vst.row_annot2[coldata$MouseType == "B6", ])

ComplexHeatmap::draw(vst.cor.B6)


vst.cor.129 <- CorHeatmap(vst[, coldata$MouseType == "129"],
                          nr_col_split = 3,
                          nr_row_split = 3,
                          top_annot = vst.col_annot2[coldata$MouseType == "129", ],
                          right_annot = vst.row_annot2[coldata$MouseType == "129", ])

ComplexHeatmap::draw(vst.cor.129)


# adipogenesis marker genes (very noisy)
adipogenesis_markers <- msigdbr::msigdbr(species = "mouse", category = "H") %>% 
    dplyr::filter(gs_name == "HALLMARK_ADIPOGENESIS") %>% 
    dplyr::select(gene_symbol, entrez_gene, ensembl_gene)

adipogenesis_markers.exprs_heatmap <- MarkersExprsHeatmap(vst, 
                                                          adipogenesis_markers$ensembl_gene,
                                                          col_annot = vst.col_annot)


adipogenesis_markers.exprs_heatmap.B6 <- MarkersExprsHeatmap(vst[, coldata$MouseType == "B6"], 
                                                             adipogenesis_markers$ensembl_gene,
                                                             col_annot = vst.col_annot2[coldata$MouseType == "B6", ])

ComplexHeatmap::draw(adipogenesis_markers.exprs_heatmap.B6)

adipogenesis_markers.exprs_heatmap.129 <- MarkersExprsHeatmap(vst[, coldata$MouseType == "129"], 
                                                              adipogenesis_markers$ensembl_gene,
                                                              col_annot = vst.col_annot2[coldata$MouseType == "129", ])

ComplexHeatmap::draw(adipogenesis_markers.exprs_heatmap.129)


# adipocyte marker genes
adipocyte_markers.safe <- biomart.cache %>% 
    dplyr::filter(external_gene_name %in% adipocyte_markers) %>% 
    dplyr::select(ensembl_gene_id, external_gene_name) %>% 
    dplyr::distinct(ensembl_gene_id, .keep_all = T) %>%
    dplyr::filter(ensembl_gene_id %in% rownames(vst))

vst.adipocyte_markers.annot <- ComplexHeatmap::HeatmapAnnotation(df = adipocyte_markers.df %>% 
                                                                     dplyr::filter(Gene %in% adipocyte_markers.safe$external_gene_name) %>%
                                                                     tibble::column_to_rownames("Gene"),
                                                                 col = list(Class = adipocyte_markers.colors),
                                                                 which = "row",
                                                                 border = T,
                                                                 show_annotation_name = F,
                                                                 show_legend = T)

adipocyte_markers.exprs_heatmap <- MarkersExprsHeatmap(vst, 
                                                       adipocyte_markers.safe$ensembl_gene_id,
                                                       scale = F,
                                                       legend_name = "Expression",
                                                       col_annot = vst.col_annot,
                                                       row_annot = vst.adipocyte_markers.annot,
                                                       show_row_names = T,
                                                       cluster_row = F,
                                                       nr_col_split = 5)

ComplexHeatmap::draw(adipocyte_markers.exprs_heatmap, padding = unit(c(2, 12, 2, 4), "mm"))

# let's separate B6 and 129 group
adipocyte_markers.B6.exprs_heatmap <- MarkersExprsHeatmap(vst[, colnames(vst)[grepl(x = colnames(vst), pattern = "B6")]], 
                                                          scale = F,
                                                          legend_name = "Norm.Exprs",
                                                          custom_col = test.col,
                                                       adipocyte_markers.safe$ensembl_gene_id,
                                                       col_annot = ComplexHeatmap::HeatmapAnnotation(df = coldata[grepl(x = rownames(coldata), pattern = "B6"), "Group", drop = F],
                                                                                                     col = list(Group = group_colors),
                                                                                                     border = T,
                                                                                                     show_annotation_name = F,
                                                                                                     show_legend = T),
                                                       row_annot = vst.adipocyte_markers.annot,
                                                       show_row_names = T,
                                                       cluster_row = F,
                                                       nr_col_split = 3)

ComplexHeatmap::draw(adipocyte_markers.B6.exprs_heatmap, padding = unit(c(2, 10, 2, 4), "mm"))

adipocyte_markers.129.exprs_heatmap <- MarkersExprsHeatmap(vst[, colnames(vst)[grepl(x = colnames(vst), pattern = "129")]],
                                                           scale = F,
                                                           legend_name = "Norm.Exprs",
                                                           custom_col = test.col,
                                                           adipocyte_markers.safe$ensembl_gene_id,
                                                           col_annot = ComplexHeatmap::HeatmapAnnotation(df = coldata[grepl(x = rownames(coldata), pattern = "B6"), "Group", drop = F],
                                                                                                     col = list(Group = group_colors),
                                                                                                     border = T,
                                                                                                     show_annotation_name = F,
                                                                                                     show_legend = T),
                                                       row_annot = vst.adipocyte_markers.annot,
                                                       show_row_names = T,
                                                       cluster_row = F,
                                                       nr_col_split = 2)

ComplexHeatmap::draw(adipocyte_markers.129.exprs_heatmap, padding = unit(c(2, 10, 2, 4), "mm"))



# TF inference using CollecTRI database
vst.ctri <- runCollecTRI(vst)
vst.ctri.sig <- vst.ctri %>% dplyr::filter(p_value < 0.05)

vst.ctri.sig.mat <- vst.ctri.sig %>%
    tidyr::pivot_wider(id_cols = "condition",
                       names_from = "source",
                       values_from = "score") %>%
    tibble::column_to_rownames(var = "condition") %>%
    as.matrix()

vst.ctr.sig.adipo.mat <- vst.ctri.sig.mat[, adipocyte_differentiation_tfs] 
vst.ctr.sig.adipo.mat.na_index <- is.na(vst.ctr.sig.adipo.mat)
vst.ctr.sig.adipo.mat[vst.ctr.sig.adipo.mat.na_index] <- 0
vst.ctr.sig.adipo.mat <- vst.ctr.sig.adipo.mat %>%
    scale() %>% 
    t()

vst.ctr.sig.adipo.mat[vst.ctr.sig.adipo.mat.na_index %>% t()] <- NA

vst.adipo_tfs.heatmap <- ComplexHeatmap::Heatmap(vst.ctr.sig.adipo.mat[, colnames(vst)],
                                                 rect_gp = gpar(col = "white", lwd = 0.5),
                                                 cluster_columns = T,
                                                 cluster_rows = F,
                                                 na_col = "#818181",
                                                 bottom_annotation = ComplexHeatmap::HeatmapAnnotation(df = coldata[colnames(vst), "MouseType", drop = F],
                                                                                                       col = list("MouseType" = group_colors2),
                                                                                                       border = T,
                                                                                                       show_annotation_name = F,
                                                                                                       show_legend = T),
                                                 right_annotation = ComplexHeatmap::HeatmapAnnotation(df = adipocyte_differentiation_tfs.df,
                                                                                                      which = "row",
                                                                                                      col = list("Sign" = adipocyte_differentiation_tfs.colors),
                                                                                                      border = T,
                                                                                                      show_annotation_name = F,
                                                                                                      show_legend = T),
                                                 row_split = NULL,
                                                 column_split = 7,
                                                 column_title = NULL,
                                                 column_names_rot = 45,
                                                 name = "Inference Score")

ComplexHeatmap::draw(vst.adipo_tfs.heatmap, padding = unit(c(2, 12, 2, 2), "mm"))

# let's separate B6 and 129 group
vst.ctr.sig.adipo.mat.B6 <- vst.ctri.sig.mat[grepl(x = rownames(vst.ctri.sig.mat), pattern = "B6"), 
                                             adipocyte_differentiation_tfs]
vst.ctr.sig.adipo.mat.B6.na_index <- is.na(vst.ctr.sig.adipo.mat.B6)
vst.ctr.sig.adipo.mat.B6[vst.ctr.sig.adipo.mat.B6.na_index] <- 0
vst.ctr.sig.adipo.mat.B6 <- vst.ctr.sig.adipo.mat.B6 %>%
    scale() %>% 
    t()

vst.ctr.sig.adipo.mat.B6[vst.ctr.sig.adipo.mat.B6.na_index %>% t()] <- NA

vst.adipo_tfs.B6.heatmap <- ComplexHeatmap::Heatmap(vst.ctr.sig.adipo.mat.B6[, colnames(vst)[grepl(x = colnames(vst), pattern = "B6")]],
                                                    rect_gp = gpar(col = "white", lwd = 0.5),
                                                    cluster_columns = T,
                                                    cluster_rows = F,
                                                    na_col = "#818181",
                                                    bottom_annotation = ComplexHeatmap::HeatmapAnnotation(df = coldata[grepl(x = rownames(coldata), pattern = "B6"), "Group", drop = F],
                                                                                                          col = list(Group = group_colors),
                                                                                                          border = T,
                                                                                                          show_annotation_name = F,
                                                                                                          show_legend = T),
                                                    right_annotation = ComplexHeatmap::HeatmapAnnotation(df = adipocyte_differentiation_tfs.df,
                                                                                                         which = "row",
                                                                                                         col = list("Sign" = adipocyte_differentiation_tfs.colors),
                                                                                                         border = T,
                                                                                                         show_annotation_name = F,
                                                                                                         show_legend = T),
                                                    row_split = NULL,
                                                    column_split = 4,
                                                    column_title = NULL,
                                                    column_names_rot = 45,
                                                    name = "Inference Score")

ComplexHeatmap::draw(vst.adipo_tfs.B6.heatmap, padding = unit(c(2, 12, 2, 2), "mm"))

vst.ctr.sig.adipo.mat.129 <- vst.ctri.sig.mat[grepl(x = rownames(vst.ctri.sig.mat), pattern = "129"), 
                                              adipocyte_differentiation_tfs]
vst.ctr.sig.adipo.mat.129.na_index <- is.na(vst.ctr.sig.adipo.mat.129)
vst.ctr.sig.adipo.mat.129[vst.ctr.sig.adipo.mat.129.na_index] <- 0
vst.ctr.sig.adipo.mat.129 <- vst.ctr.sig.adipo.mat.129 %>%  # dim : (14, 20) ; sample x gene
    scale() %>%  # scale columns --> gene-wise scaling for sample-wise comparison
    t()  # dim : (20, 14) ; gene x sample
    

vst.ctr.sig.adipo.mat.129[vst.ctr.sig.adipo.mat.129.na_index %>% t()] <- NA

vst.adipo_tfs.129.heatmap <- ComplexHeatmap::Heatmap(vst.ctr.sig.adipo.mat.129[, colnames(vst)[grepl(x = colnames(vst), pattern = "129")]],
                                                    rect_gp = gpar(col = "white", lwd = 0.5),
                                                    cluster_columns = T,
                                                    cluster_rows = F,
                                                    na_col = "#818181",
                                                    bottom_annotation = ComplexHeatmap::HeatmapAnnotation(df = coldata[grepl(x = rownames(coldata), pattern = "129"), "Group", drop = F],
                                                                                                          col = list(Group = group_colors),
                                                                                                          border = T,
                                                                                                          show_annotation_name = F,
                                                                                                          show_legend = T),
                                                    right_annotation = ComplexHeatmap::HeatmapAnnotation(df = adipocyte_differentiation_tfs.df,
                                                                                                         which = "row",
                                                                                                         col = list("Sign" = adipocyte_differentiation_tfs.colors),
                                                                                                         border = T,
                                                                                                         show_annotation_name = F,
                                                                                                         show_legend = T),
                                                    row_split = NULL,
                                                    column_split = 3,
                                                    column_title = NULL,
                                                    column_names_rot = 45,
                                                    name = "Inference Score")

ComplexHeatmap::draw(vst.adipo_tfs.129.heatmap, padding = unit(c(2, 12, 2, 2), "mm"))



# GSVA

target_gs <- c("HALLMARK_ADIPOGENESIS", 
               "WP_ADIPOGENESIS", 
               "WP_TRANSCRIPTIONAL_CASCADE_REGULATING_ADIPOGENESIS",
               "WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS",
               "WP_DIFFERENTIATION_OF_WHITE_AND_BROWN_ADIPOCYTE",
               "GOBP_ADIPOSE_TISSUE_DEVELOPMENT",
               "GOBP_REGULATION_OF_ADIPOSE_TISSUE_DEVELOPMENT",
               "GOBP_POSITIVE_REGULATION_OF_ADIPOSE_TISSUE_DEVELOPMENT",
               "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_WHITE_ADIPOCYTE_DIFFERENTIATION",
               "URS_ADIPOCYTE_DIFFERENTIATION_UP",
               "URS_ADIPOCYTE_DIFFERENTIATION_DN",
               "LI_ADIPOGENESIS_BY_ACTIVATED_PPARG",
               "NAKAMURA_ADIPOGENESIS_LATE_UP", 
               "NAKAMURA_ADIPOGENESIS_LATE_DN")

vst.gsva <- runGSVA(vst, target_gs)  # dim : (13, 28) ; geneset x sample

vst.gsva.heatmap <- ComplexHeatmap::Heatmap(vst.gsva %>% t() %>% scale() %>% t(),
                                            rect_gp = gpar(col = "white", lwd = 0.1),
                                            row_names_gp = gpar(fontsize = 8),
                                            column_names_gp = gpar(fontsize = 8),
                                            cluster_rows = T,
                                            cluster_columns = T,
                                            #row_split = 5,
                                            column_split = 4,
                                            row_title = NULL,
                                            column_title = NULL,
                                            column_names_rot = 45,
                                            bottom_annotation = ComplexHeatmap::HeatmapAnnotation(df = coldata[colnames(vst), "MouseType", drop = F],
                                                                                                  col = list("MouseType" = group_colors2),
                                                                                                  border = T,
                                                                                                  show_annotation_name = F,
                                                                                                  show_legend = F),
                                            show_heatmap_legend = F,
                                            #show_row_names = F,
                                            name = "Signature Score")

ComplexHeatmap::draw(vst.gsva.heatmap, padding = unit(c(2, 5, 2, 70), "mm"))

# let's separate B6 and 129 group
vst.gsva.B6.heatmap <- ComplexHeatmap::Heatmap(vst.gsva[, coldata$MouseType == "B6"] %>% t() %>% scale() %>% t(),
                                              rect_gp = gpar(col = "white", lwd = 0.1),
                                              row_names_gp = gpar(fontsize = 8),
                                              column_names_gp = gpar(fontsize = 8),
                                              cluster_rows = T,
                                              cluster_columns = T,
                                              #row_split = 5,
                                              column_split = 3,
                                              row_title = NULL,
                                              column_title = NULL,
                                              column_names_rot = 45,
                                              bottom_annotation = ComplexHeatmap::HeatmapAnnotation(df = coldata[grepl(x = rownames(coldata), pattern = "B6"), "Group", drop = F],
                                                                                                    col = list(Group = group_colors),
                                                                                                    border = T,
                                                                                                    show_annotation_name = F,
                                                                                                    show_legend = F),
                                              show_heatmap_legend = F,
                                              #show_row_names = F,
                                              name = "Signature Score")

ComplexHeatmap::draw(vst.gsva.B6.heatmap, padding = unit(c(2, 5, 2, 70), "mm"))

vst.gsva.129.heatmap <- ComplexHeatmap::Heatmap(vst.gsva[, coldata$MouseType == "129"] %>% t() %>% scale() %>% t(),
                                               rect_gp = gpar(col = "white", lwd = 0.1),
                                               row_names_gp = gpar(fontsize = 8),
                                               column_names_gp = gpar(fontsize = 8),
                                               cluster_rows = T,
                                               cluster_columns = T,
                                               #row_split = 5,
                                               column_split = 4,
                                               row_title = NULL,
                                               column_title = NULL,
                                               column_names_rot = 45,
                                               bottom_annotation = ComplexHeatmap::HeatmapAnnotation(df = coldata[grepl(x = rownames(coldata), pattern = "129"), "Group", drop = F],
                                                                                                     col = list(Group = group_colors),
                                                                                                     border = T,
                                                                                                     show_annotation_name = F,
                                                                                                     show_legend = T),
                                               show_heatmap_legend = T,
                                               #show_row_names = F,
                                               name = "Signature Score")

ComplexHeatmap::draw(vst.gsva.129.heatmap, padding = unit(c(2, 5, 2, 70), "mm"))


# B6 group differential expression analysis
dds.B6 <- dds[, grepl("B6", colnames(dds))]
dds.B6$Group <- droplevels(dds.B6$Group)
# (https://support.bioconductor.org/p/118090/) : To compare one vs all the rest, you don’t want an intercept. Use ~0 + condition.
design(dds.B6) <- ~0 + Group

ddsX.B6 <- DESeq2::DESeq(dds.B6)

vst.B6 <- DESeq2::vst(ddsX.B6, blind = F)

vst.B6.pca <- PCAplot(vst.B6, 
                      grp_col = group_colors,
                      pt_sz = 8, lab_sz = 4, glob_txt_sz = 16)


# (https://support.bioconductor.org/p/86347/) : You could also specify this with a numeric contrast, e.g. contrast=c(0,1,-1/5,-1/5,-1/5,-1/5,-1/5).
B6.Im.res <- DESeq2::lfcShrink(ddsX.B6, contrast = c(1, -1/6, -1/6, -1/6, -1/6, -1/6, -1/6), type = "ashr") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ENSG") %>%
    Ens2Sym(., "ENSG")

B6.6.res <- DESeq2::lfcShrink(ddsX.B6, contrast = c(-1/6, 1, -1/6, -1/6, -1/6, -1/6, -1/6), type = "ashr") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ENSG") %>%
    Ens2Sym(., "ENSG")

B6.20.res <- DESeq2::lfcShrink(ddsX.B6, contrast = c(-1/6, -1/6, 1, -1/6, -1/6, -1/6, -1/6), type = "ashr") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ENSG") %>%
    Ens2Sym(., "ENSG")

B6.27.res <- DESeq2::lfcShrink(ddsX.B6, contrast = c(-1/6, -1/6, -1/6, 1, -1/6, -1/6, -1/6), type = "ashr") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ENSG") %>%
    Ens2Sym(., "ENSG")

B6.28.res <- DESeq2::lfcShrink(ddsX.B6, contrast = c(-1/6, -1/6, -1/6, -1/6, 1, -1/6, -1/6), type = "ashr") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ENSG") %>%
    Ens2Sym(., "ENSG")

B6.29.res <- DESeq2::lfcShrink(ddsX.B6, contrast = c(-1/6, -1/6, -1/6, -1/6, -1/6, 1, -1/6), type = "ashr") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ENSG") %>%
    Ens2Sym(., "ENSG")

B6.34.res <- DESeq2::lfcShrink(ddsX.B6, contrast = c(-1/6, -1/6, -1/6, -1/6, -1/6, -1/6, 1), type = "ashr") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ENSG") %>%
    Ens2Sym(., "ENSG")

B6.Im.sig_res <- B6.Im.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)  # 2598

B6.6.sig_res <- B6.6.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)  # 2908

B6.20.sig_res <- B6.20.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)  # 3446

B6.27.sig_res <- B6.27.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)  # 2426

B6.28.sig_res <- B6.28.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)  # 1889

B6.29.sig_res <- B6.29.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)  # 2011

B6.34.sig_res <- B6.34.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)  # 1949


# 129 group differential expression analysis
dds.129 <- dds[, grepl("129", colnames(dds))]
dds.129$Group <- droplevels(dds.129$Group)
design(dds.129) <- ~0 + Group

ddsX.129 <- DESeq2::DESeq(dds.129)

vst.129 <- DESeq2::vst(ddsX.129, blind = F)

vst.129.pca <- PCAplot(vst.129, 
                       grp_col = group_colors,
                       pt_sz = 8, lab_sz = 4, glob_txt_sz = 16)





# IM-129 vs. each sample (consider IM-129 as a reference level)
levels(dds$Group) <- c("IM-129", 
                       levels(dds$Group)[levels(dds$Group) != "IM-129"])

ddsX <- DESeq2::DESeq(dds)

IM129_vs_B6.Im.res <- DESeq2::lfcShrink(ddsX, coef = "Group_IM.B6_vs_IM.129", type = "apeglm")
IM129_vs_B6.Im.res <- DDSAddSymbols(IM129_vs_B6.Im.res)

IM129_vs_B6.6.res <- DESeq2::lfcShrink(ddsX, coef = "Group_B6.6_vs_IM.129", type = "apeglm")
IM129_vs_B6.6.res <- DDSAddSymbols(IM129_vs_B6.6.res)

IM129_vs_B6.20.res <- DESeq2::lfcShrink(ddsX, coef = "Group_B6.20_vs_IM.129", type = "apeglm")
IM129_vs_B6.20.res <- DDSAddSymbols(IM129_vs_B6.20.res)

IM129_vs_B6.27.res <- DESeq2::lfcShrink(ddsX, coef = "Group_B6.27_vs_IM.129", type = "apeglm")
IM129_vs_B6.27.res <- DDSAddSymbols(IM129_vs_B6.27.res)

IM129_vs_B6.28.res <- DESeq2::lfcShrink(ddsX, coef = "Group_B6.28_vs_IM.129", type = "apeglm")
IM129_vs_B6.28.res <- DDSAddSymbols(IM129_vs_B6.28.res)

IM129_vs_B6.29.res <- DESeq2::lfcShrink(ddsX, coef = "Group_B6.29_vs_IM.129", type = "apeglm")
IM129_vs_B6.29.res <- DDSAddSymbols(IM129_vs_B6.29.res)

IM129_vs_B6.34.res <- DESeq2::lfcShrink(ddsX, coef = "Group_B6.34_vs_IM.129", type = "apeglm")
IM129_vs_B6.34.res <- DDSAddSymbols(IM129_vs_B6.34.res)

IM129_vs_B6.Im.sig_res <- IM129_vs_B6.Im.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

IM129_vs_B6.6.sig_res <- IM129_vs_B6.6.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

IM129_vs_B6.20.sig_res <- IM129_vs_B6.20.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

IM129_vs_B6.27.sig_res <- IM129_vs_B6.27.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

IM129_vs_B6.28.sig_res <- IM129_vs_B6.28.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

IM129_vs_B6.29.sig_res <- IM129_vs_B6.29.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

IM129_vs_B6.34.sig_res <- IM129_vs_B6.34.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

writexl::write_xlsx(IM129_vs_B6.Im.sig_res, path = "results/IM129_vs_B6.Im.sig_res.xlsx")
writexl::write_xlsx(IM129_vs_B6.6.sig_res, path = "results/IM129_vs_B6.6.sig_res.xlsx")
writexl::write_xlsx(IM129_vs_B6.20.sig_res, path = "results/IM129_vs_B6.20.sig_res.xlsx")
writexl::write_xlsx(IM129_vs_B6.27.sig_res, path = "results/IM129_vs_B6.27.sig_res.xlsx")
writexl::write_xlsx(IM129_vs_B6.28.sig_res, path = "results/IM129_vs_B6.28.sig_res.xlsx")
writexl::write_xlsx(IM129_vs_B6.29.sig_res, path = "results/IM129_vs_B6.29.sig_res.xlsx")
writexl::write_xlsx(IM129_vs_B6.34.sig_res, path = "results/IM129_vs_B6.34.sig_res.xlsx")

IM129_vs_B6.Im.enr <- DDSEnrich(IM129_vs_B6.Im.res)
IM129_vs_B6.6.enr <- DDSEnrich(IM129_vs_B6.6.res)
IM129_vs_B6.20.enr <- DDSEnrich(IM129_vs_B6.20.res)
IM129_vs_B6.27.enr <- DDSEnrich(IM129_vs_B6.27.res)
IM129_vs_B6.28.enr <- DDSEnrich(IM129_vs_B6.28.res)
IM129_vs_B6.29.enr <- DDSEnrich(IM129_vs_B6.29.res)
IM129_vs_B6.34.enr <- DDSEnrich(IM129_vs_B6.34.res)


IM129_vs_129.3.res <- DESeq2::lfcShrink(ddsX, coef = "Group_129.3_vs_IM.129", type = "apeglm")
IM129_vs_129.3.res <- DDSAddSymbols(IM129_vs_129.3.res)

IM129_vs_129.9.res <- DESeq2::lfcShrink(ddsX, coef = "Group_129.9_vs_IM.129", type = "apeglm")
IM129_vs_129.9.res <- DDSAddSymbols(IM129_vs_129.9.res)

IM129_vs_129.10.res <- DESeq2::lfcShrink(ddsX, coef = "Group_129.10_vs_IM.129", type = "apeglm")
IM129_vs_129.10.res <- DDSAddSymbols(IM129_vs_129.10.res)

IM129_vs_129.16.res <- DESeq2::lfcShrink(ddsX, coef = "Group_129.16_vs_IM.129", type = "apeglm")
IM129_vs_129.16.res <- DDSAddSymbols(IM129_vs_129.16.res)

IM129_vs_129.22.res <- DESeq2::lfcShrink(ddsX, coef = "Group_129.22_vs_IM.129", type = "apeglm")
IM129_vs_129.22.res <- DDSAddSymbols(IM129_vs_129.22.res)

IM129_vs_129.23.res <- DESeq2::lfcShrink(ddsX, coef = "Group_129.23_vs_IM.129", type = "apeglm")
IM129_vs_129.23.res <- DDSAddSymbols(IM129_vs_129.23.res)


IM129_vs_129.3.sig_res <- IM129_vs_129.3.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

IM129_vs_129.9.sig_res <- IM129_vs_129.9.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

IM129_vs_129.10.sig_res <- IM129_vs_129.10.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

IM129_vs_129.16.sig_res <- IM129_vs_129.16.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

IM129_vs_129.22.sig_res <- IM129_vs_129.22.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

IM129_vs_129.23.sig_res <- IM129_vs_129.23.res %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

writexl::write_xlsx(IM129_vs_129.3.sig_res, path = "results/IM129_vs_129.3.sig_res.xlsx")
writexl::write_xlsx(IM129_vs_129.9.sig_res, path = "results/IM129_vs_129.9.sig_res.xlsx")
writexl::write_xlsx(IM129_vs_129.10.sig_res, path = "results/IM129_vs_129.10.sig_res.xlsx")
writexl::write_xlsx(IM129_vs_129.16.sig_res, path = "results/IM129_vs_129.16.sig_res.xlsx")
writexl::write_xlsx(IM129_vs_129.22.sig_res, path = "results/IM129_vs_129.22.sig_res.xlsx")
writexl::write_xlsx(IM129_vs_129.23.sig_res, path = "results/IM129_vs_129.23.sig_res.xlsx")


IM129_vs_B6.Im.enr <- DDSEnrich(IM129_vs_B6.Im.res)
IM129_vs_B6.6.enr <- DDSEnrich(IM129_vs_B6.6.res)
IM129_vs_B6.20.enr <- DDSEnrich(IM129_vs_B6.20.res)
IM129_vs_B6.27.enr <- DDSEnrich(IM129_vs_B6.27.res)
IM129_vs_B6.28.enr <- DDSEnrich(IM129_vs_B6.28.res)
IM129_vs_B6.29.enr <- DDSEnrich(IM129_vs_B6.29.res)



#======================#
#   get regulon info   #
#======================#========================================================
regulons <- decoupleR::get_collectri(organism = "mouse")



#====================#
#   B6-27 vs B6-28   #
#====================#
# Use Contrast option in dds.B6 rather than subsetting only target samples to consider dispersion within B6 group
#===============================================================================
MTB6_27_vs_MTB6_28.vst <- vst.B6[, c("B6-27_1", "B6-27_2", "B6-28_1", "B6-28_2")]

MTB6_27_vs_MTB6_28.res <- DESeq2::lfcShrink(dds = ddsX.B6, 
                                            contrast = c("Group", "B6_28", "B6_27"), 
                                            type = "ashr") # genes with logFC > 0 are overexpressed in B6-28
MTB6_27_vs_MTB6_28.res <- DDSAddSymbols(MTB6_27_vs_MTB6_28.res)

MTB6_27_vs_MTB6_28.sig_res <- MTB6_27_vs_MTB6_28.res %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    dplyr::filter(!grepl("^Gm[0-9]", external_gene_name) & !grepl("[0-9]Rik", external_gene_name))

writexl::write_xlsx(MTB6_27_vs_MTB6_28.sig_res, path = "results/DEG_B6_27_vs_28//MTB6_27_vs_MTB6_28.sig_res.xlsx")

MTB6_27_vs_MTB6_28.res$padj_pc <- MTB6_27_vs_MTB6_28.res$padj + .Machine$double.xmin  # add pseudo-count to avoid 0 pvalues
MTB6_27_vs_MTB6_28.top_genes <- c(MTB6_27_vs_MTB6_28.sig_res %>% dplyr::top_n(n = 10, log2FoldChange) %>% dplyr::pull(external_gene_name),
                                  MTB6_27_vs_MTB6_28.sig_res %>% dplyr::top_n(n = -10, log2FoldChange) %>% dplyr::pull(external_gene_name))

Volcanoplot(res_df = MTB6_27_vs_MTB6_28.res,
            y = "padj_pc",
            target_genes = MTB6_27_vs_MTB6_28.top_genes,
            label_size = 5)

MTB6_27_vs_MTB6_28.sig_res %>% 
    dplyr::filter(external_gene_name %in% intersect(MTB6_27_vs_MTB6_28.sig_res$external_gene_name, regulons$source)) %>%
    View()

MarkersExprsHeatmap(X_vst = MTB6_27_vs_MTB6_28.vst, 
                    ENSG_Markers = biomart.cache %>% 
                        dplyr::filter(external_gene_name %in% intersect(MTB6_27_vs_MTB6_28.sig_res$external_gene_name, regulons$source)) %>% 
                        dplyr::pull(ensembl_gene_id) %>% 
                        unique())



#======================#
#   129-22 vs 129-23   #
#======================#========================================================
MT129_22_vs_MT129_23.vst <- vst.129[, c("129-22_1", "129-22_2", "129-23_1", "129-23_2")]

MT129_22_vs_MT129_23.res <- DESeq2::lfcShrink(ddsX.129, contrast = c("Group", "129_23", "129_22"), type = "ashr")  # genes with logFC > 0 are overexpressed in 129_23
MT129_22_vs_MT129_23.res <- DDSAddSymbols(MT129_22_vs_MT129_23.res)

MT129_22_vs_MT129_23.sig_res <- MT129_22_vs_MT129_23.res %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    dplyr::filter(!grepl("^Gm[0-9]", external_gene_name) & !grepl("[0-9]Rik", external_gene_name))

writexl::write_xlsx(MT129_22_vs_MT129_23.sig_res, path = "results/DEG_129_22_vs_23/MT129_22_vs_MT129_23.sig_res.xlsx")

MT129_22_vs_MT129_23.res$padj_pc <- MT129_22_vs_MT129_23.res$padj + .Machine$double.xmin  # add pseudo-count to avoid 0 pvalues
MT129_22_vs_MT129_23.top_genes <- c(MT129_22_vs_MT129_23.sig_res %>% dplyr::top_n(n = 10, log2FoldChange) %>% dplyr::pull(external_gene_name),
                                    MT129_22_vs_MT129_23.sig_res %>% dplyr::top_n(n = -10, log2FoldChange) %>% dplyr::pull(external_gene_name))

Volcanoplot(res_df = MT129_22_vs_MT129_23.res,
            y = "padj_pc",
            target_genes = MT129_22_vs_MT129_23.top_genes,
            label_size = 5)

MT129_22_vs_MT129_23.sig_res %>% 
    dplyr::filter(external_gene_name %in% intersect(MT129_22_vs_MT129_23.sig_res$external_gene_name, regulons$source)) %>%
    View()

MarkersExprsHeatmap(X_vst = MT129_22_vs_MT129_23.vst, 
                    ENSG_Markers = biomart.cache %>% 
                        dplyr::filter(external_gene_name %in% intersect(MT129_22_vs_MT129_23.sig_res$external_gene_name, regulons$source)) %>% 
                        dplyr::pull(ensembl_gene_id) %>% 
                        unique())



#======================================#
#   B6-27 & 129-22 vs B6-28 & 129-23   #
#======================================#
# B6-27 & 129-22 : Differentiation (++) ; Browning (++)
# B6-28 & 129-23 : Differentiation (++) ; Browning (--)
#===============================================================================
# add Differentiation and Browning information to metadata
coldata$AD_CL <- c("AD_N.CL_N", "AD_N.CL_N",  # IM-B6 --> AD(-), CL(-)
                   "AD_N.CL_N", "AD_N.CL_N",  # B6-6 --> AD(-), CL(-)
                   "AD_N.CL_N", "AD_N.CL_N",  # B6-20 --> AD(-), CL(-)
                   "AD_P.CL_P", "AD_P.CL_P",  # B6-27 --> AD(+), CL(+)
                   "AD_P.CL_N", "AD_P.CL_N",  # B6-28 --> AD(+), CL(-)
                   "AD_N.CL_N", "AD_N.CL_N",  # B6-29 --> AD(-), CL(-)
                   "AD_N.CL_N", "AD_N.CL_N",  # B6-34 --> AD(-), CL(-)
                   "AD_N.CL_N", "AD_N.CL_N",  # IM-129 --> AD(-), CL(-)
                   "AD_N.CL_N", "AD_N.CL_N",  # 129-3 --> AD(-), CL(-)
                   "AD_N.CL_N", "AD_N.CL_N",  # 129-9 --> AD(-), CL(-)
                   "AD_N.CL_N", "AD_N.CL_N",  # 129-10 --> AD(-), CL(-)
                   "AD_N.CL_N", "AD_N.CL_N",  # 129-16 --> AD(-), CL(-)
                   "AD_P.CL_P", "AD_P.CL_P",  # 129-22 --> AD(+), CL(+)
                   "AD_P.CL_N", "AD_P.CL_N")  # 129-23 --> AD(+), CL(-)

coldata$AD_CL <- factor(coldata$AD_CL, levels = c("AD_N.CL_N", "AD_P.CL_N", "AD_P.CL_P"))

dds.AD_CL <- DESeq2::DESeqDataSetFromTximport(txi = txi, 
                                              colData = coldata, 
                                              design = ~AD_CL)

dds.AD_CL <- dds.AD_CL[rowSums(BiocGenerics::counts(dds.AD_CL) >= 10) > 3, ]

ddsX.AD_CL <- DESeq2::DESeq(dds.AD_CL)

vst.AD_CL <- DESeq2::vst(ddsX.AD_CL, blind = F)

AD_CL.res <- DESeq2::lfcShrink(ddsX.AD_CL, contrast = c("AD_CL", "AD_P.CL_P", "AD_P.CL_N"), type = "ashr")  # genes with logFC > 0 are overexpressed in AD_P.CL_P
AD_CL.res <- DDSAddSymbols(AD_CL.res)

AD_CL.sig_res <- AD_CL.res %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    dplyr::filter(!grepl("^Gm[0-9]", external_gene_name) & !grepl("[0-9]Rik", external_gene_name))

writexl::write_xlsx(AD_CL.sig_res, path = "results/DEG_129_22_vs_23/AD_CL.sig_res.xlsx")

AD_CL.res$padj_pc <- AD_CL.res$padj + .Machine$double.xmin  # add pseudo-count to avoid 0 pvalues
AD_CL.top_genes <- c(AD_CL.sig_res %>% dplyr::top_n(n = 10, log2FoldChange) %>% dplyr::pull(external_gene_name),
                     AD_CL.sig_res %>% dplyr::top_n(n = -10, log2FoldChange) %>% dplyr::pull(external_gene_name))

Volcanoplot(res_df = AD_CL.res,
            y = "padj_pc",
            target_genes = AD_CL.sig_res$external_gene_name,
            label_size = 5,
            ylim = c(0, 8))

MT129_22_vs_MT129_23.sig_res %>% 
    dplyr::filter(external_gene_name %in% intersect(MT129_22_vs_MT129_23.sig_res$external_gene_name, regulons$source)) %>%
    View()

MarkersExprsHeatmap(X_vst = MT129_22_vs_MT129_23.vst, 
                    ENSG_Markers = biomart.cache %>% 
                        dplyr::filter(external_gene_name %in% intersect(MT129_22_vs_MT129_23.sig_res$external_gene_name, regulons$source)) %>% 
                        dplyr::pull(ensembl_gene_id) %>% 
                        unique())













