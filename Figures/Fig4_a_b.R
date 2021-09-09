# Code used in the study (Figure 4 A and B): "Molecular characterization of colorectal cancer related peritoneal metastatic disease", Kristiaan J. Lenos et al.
# Please direct any questions or comments related to this code to Leandro Moreno: l.ferreiramoreno@amsterdamumc.nl

#Figure 4 A

cond <- read.table("Files/prim_vs_pm.txt", sep="\t",header = TRUE,row.names = 1,check.names=FALSE) 
cond$cond <- factor(cond$cond)
cond$pat <- factor(cond$pat)


prim_vs_pm <- subset(assay_peri_df, select=rownames(cond))
matrix_prim_vs_pm <- prim_vs_pm[, rownames(cond)]
all(rownames(cond) == colnames(prim_vs_pm))

dds <- DESeqDataSetFromMatrix(countData = prim_vs_pm,
                              colData = cond,
                              design = ~ pat + cond)

vsd_pm <- vst(dds, blind=FALSE)
sp2 <- plotPCA(vsd_pm, intgroup=c("cond")) + theme_bw() + geom_text_repel(aes(label=name),vjust=2, size = 2) 
sp2 <- sp2 + scale_color_manual(values=c("#cd262f", "#5798c9"))

#Figure 4 B

#Liver primary and metastasis data were downloaded from the Recount2 database; ID: SRP029880
load.Rdata( filename="C:/amc.intra/users/L/lferreiramoreno/home/rdata/rse_54study.Rdata", "dat.54study" ) 

patterns <- c("primary", "metastasized")
sample_info <- data.frame(colData(dat.54study)[grepl(paste(patterns, collapse="|"), colData(dat.54study)$title), ])
sample_info$patient=gsub("-.*","",gsub(".*_","pat",sample_info$title))
sample_info$title=gsub("c.*","",sample_info$title)
sample_info$connectors <- sample_info$patient

colData(dat.54study.subset)$group <- sample_info$title
colData(dat.54study.subset)$cmsclass <- classification2$SSP.nearestCMS
colData(dat.54study.subset)$patient <- sample_info$patient
colData(dat.54study.subset)$connectors <- sample_info$connectors

dat.54study.subset <- subset(dat.54study, select = colData(dat.54study)$sample %in% sample_info$sample)
dds <- DESeqDataSet(dat.54study.subset, ~ group)
vsd <- vst(dds, blind=FALSE)
vsd_54 <- assay(vsd)

sp54 <- plotPCA(vsd_54, intgroup=c("group")) + theme_bw() + geom_text_repel(aes(label=name), size = 2, vjust=2) #vsd_54 was produced by the script 54_study
sp54 <- sp54 + scale_color_manual(values=c("#cd262f", "#5798c9"))

#CMS Classification of Liver samples
dat.54study.subset$ENTREID <- mapIds(org.Hs.eg.db,keys=rownamesdat.54study.subset,column="ENTREZID",keytype="SYMBOL",multiVals="first")
dat.54study.subset = data54_study[!duplicated(dat.54study.subset$ENTREID),]
dat.54study.subset <- dat.54study.subset[!is.na(dat.54study.subset$ENTREID),]
row.names(dat.54study.subset) <- dat.54study.subset$ENTREID
dat.54study.subset$ENTREID <- NULL

Rssp_54study_SSP <- CMSclassifier::classifyCMS(dat.54study.subset, method="SSP")[[3]]
