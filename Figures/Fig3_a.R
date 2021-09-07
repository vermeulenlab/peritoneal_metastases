# Code used in the study (Figure 3 A): "Molecular characterization of colorectal cancer related peritoneal metastatic disease", Kristiaan J. Lenos et al.
# Please direct any questions or comments related to this code to Leandro Moreno: l.ferreiramoreno@amsterdamumc.nl


#Subgroups of CMS4 samples were compared for detecting DE genes using the ANOVA test available on the platform R2 (Genomics Analysis and Visualization Platform)
samples_CMS4 <- read.table("/Files/ANOVA_CMS4_45samples.txt", sep="\t",header = TRUE,row.names = 1,check.names=FALSE)

metadata<- read.table("/Files/metadata.txt", sep="\t",header = TRUE,row.names = 1,check.names=FALSE)
metadata<-t(metadata)
metadata<-as.data.frame(metadata)
metadata <- subset(annotation_groups, select=c("groups"))
metadata_45 <- subset(metadata, rownames(metadata) %in% colnames(samples_CMS4))

metadata_45_sorted <- as.data.frame(metadata_45[colnames(samples_CMS4),])
colnames(metadata_45_sorted) <- "groups"

ha_k_45 = HeatmapAnnotation(df = metadata_45_sorted,
                            col = list(groups = c("k1"="#e7f0c3", "k2"="#a4d4ae", "k3"="#32afa9")))


Heatmap(samples_CMS4, clustering_distance_rows = "pearson",
        column_split=metadata_45_sorted, cluster_column_slices = TRUE, 
        col = col_fun,heatmap_height = unit(12, "cm"), heatmap_width = unit(7, "cm"),
        top_annotation = ha_k_45, border = TRUE,  show_row_names = FALSE,
        show_column_names = TRUE,
        show_row_dend = FALSE,
        show_column_dend = TRUE,
        column_title = "1355 DE genes; CMS4 samples",
        heatmap_legend_param = list(
          at = c(-2, 0, 2),
          title = "z-score",
          legend_height = unit(4, "cm"),
          title_position = "lefttop-rot"
        ))
