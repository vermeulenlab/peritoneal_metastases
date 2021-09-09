# Code used in the study (Figure 1 A): "Molecular characterization of colorectal cancer related peritoneal metastatic disease", Kristiaan J. Lenos et al.
# Please direct any questions or comments related to this code to Leandro Moreno: l.ferreiramoreno@amsterdamumc.nl

#CMS classification
library("tximport")
library("DESeq2")

samples <- read.table(file.path("/dir", "files_sorted.txt"), header = TRUE)
files <- file.path("count_RSEM", paste0(samples$ID))
samples$condition <- factor(rep(c("Round1","Round2"),each=90))
names(files) <- paste0("sample", 1:180)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

ddsTxi <- DESeqDataSetFromTximport(txi.rsem,
                                   colData = samples,
                                   design = ~ condition)

rlog.ddsTxiCOL = rlog(ddsTxiCOL)
assay_peri <- assay(rlog.ddsTxiCOL)
assay_peri_df <- as.data.frame(assay_peri)

#Convert SYMBOL to ENTREZID
library(org.Hs.eg.db)

assay_peri_df$ENTREID <- mapIds(org.Hs.eg.db,keys=rownames(assay_peri_df),column="ENTREZID",keytype="SYMBOL",multiVals="first")
assay_peri_df = assay_peri_df[!duplicated(assay_peri_df$ENTREID),]
assay_peri_df <- assay_peri_df[!is.na(assay_peri_df$ENTREID),]
row.names(assay_peri_df) <- assay_peri_df$ENTREID
assay_peri_df$ENTREID <- NULL

library(CMSclassifier)
Rfcms_SSP <- CMSclassifier::classifyCMS(assay_peri_df, method="SSP")[[3]]

#Generate Heatmap based on CMS classification
library(ComplexHeatmap)

metadata<- read.table("/Files/metadata.txt", sep="\t",header = TRUE,row.names = 1,check.names=FALSE)
metadata_t<-t(metadata)
metadata_t<-as.data.frame(metadata_t)
metadata_t<- metadata_t[,1, drop=FALSE]
metadata_t$CMS_SSP_nearest <- as.factor(metadata_t$CMS_SSP_nearest)

#Select top 500 variable genes 

assay_peri_SSP <- assay_peri_df[(row.names(assay_peri_df) %in% CMSclassifier::listModelGenes("SSP")), ]
z_assay_peri_SSP <- apply(assay_peri_SSP,1,scale)
rownames(z_assay_peri_SSP) <- colnames(assay_peri_SSP)
z_assay_peri_SSPmads=apply(t(z_assay_peri_SSP),1,mad)
z_assay_peri_SSP_500=t(z_assay_peri_SSP)[rev(order(z_assay_peri_SSPmads))[1:500],]

prob_CMS<- read.table("/Files/annotation_PM_prob2.txt", sep="\t",header = TRUE,row.names = 1,check.names=FALSE)
prob_CMS_subset <- prob_CMS[(row.names(prob_CMS) %in% row.names(metadata_t)), ]
prob_CMS_subset_sorted <- prob_CMS_subset[match(rownames(z_assay_peri_SSP_500), rownames(prob_CMS_subset)), ]
bottom_cor = HeatmapAnnotation(Correlation = anno_lines(prob_CMS_subset_sorted, which = c("column"),
                                                      height = unit(2, "cm"),gp = gpar(col = c("#e79e25", "#0f71ae", "#c47da1", "#169c73"), lwd = 1.3), 
                                                      add_points = TRUE, pt_gp = gpar(col = c("#e79e25", "#0f71ae", "#c47da1", "#169c73")), size = unit(1, "mm")
))

top_bar= HeatmapAnnotation(df = annotation_t,
                         col = list(CMS_SSP_nearest = c("cms1" = "#e79e25", "cms2" = "#0f71ae", "cms3" = "#c47da1", "cms4" = "#169c73")))

col_fun = colorRamp2(c(-2, 0, 2), c("#5798c9", "white", "#cd262f"))

Heatmap(z_assay_peri_SSP_500, clustering_distance_rows = "pearson",
        column_split=annotation_t, cluster_column_slices = TRUE, 
        col = col_fun,
        top_annotation = top_bar, border = TRUE,  show_row_names = FALSE,
        show_column_names = FALSE,
        show_row_dend = FALSE,
        show_column_dend = TRUE,
        bottom_annotation = bottom_cor,
        column_title = "SSP - 500 genes; 82 samples",
        heatmap_legend_param = list(
          at = c(-2, 0, 2),
          title = "z-score",
          legend_height = unit(4, "cm"),
          title_position = "lefttop-rot"
        ))



