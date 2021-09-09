# Code used in the study (Figure 5 C and D): "Molecular characterization of colorectal cancer related peritoneal metastatic disease", Kristiaan J. Lenos et al.
# Please direct any questions or comments related to this code to Leandro Moreno: l.ferreiramoreno@amsterdamumc.nl

peri_alluvial <- read.table("/Files/alluvial_peri.txt", sep="\t",header = TRUE,check.names=FALSE)
liver_alluvial <- read.table("/Files/alluvial_liver.txt", sep="\t",header = TRUE,check.names=FALSE)


peri <- ggplot(data = peri_alluvial,
               aes(axis1 = colon, axis2 = peritoneal,y = freq)) +
  scale_x_discrete(limits = c("Colon", "peritoneal"), expand = c(.2, .05)) +
  xlab("Subtype dynamic change in Peritoneal") +
  geom_alluvium(aes(fill = color), width = 1/12,  alpha = 1) +
  geom_stratum(width = 1/6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values=c("#c47da1", '#169c73'), breaks=c("cms3","cms4"), labels=c("cms3", "cms4"))

liver <- ggplot(data = liver_alluvial,
                aes(axis1 = colon, axis2 = liver,y = freq)) +
  scale_x_discrete(limits = c("Colon", "liver"), expand = c(.2, .05)) +
  xlab("Subtype dynamic change in Liver") +
  geom_alluvium(aes(fill = color), width = 1/12,  alpha = 1) +
  geom_stratum(width = 1/6) +
  scale_fill_manual(values=c("#e79e25", '#0f71ae', "#c47da1", "#169c73" ), breaks=c("cms1", "cms2", "cms3", "cms4"), labels=c("cms1", "cms2", "cms3", "cms4")) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)))

pdf("alluvial_plot.pdf", width=16, height=6)
ggarrange(peri, liver , ncol = 2, nrow = 1)
dev.off()
