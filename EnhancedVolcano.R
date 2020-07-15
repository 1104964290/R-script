# This script is used to generate volcano plot for RNA-seq data.
# The input file is the output of edgeR.
library(EnhancedVolcano)
setwd("~/ShaoXieXiangRNAseq/EnhancedVolcano/")
inputdata <- read.table("B_vs_E_all.out",header=TRUE)
# Plot volcano plot.
pdf("VolcanoPlot.pdf")
EnhancedVolcano(inputdata,lab = as.vector(inputdata$gene_name),x = 'log2FC',y = 'P_value',
selectLab = c("MYH3","MYH1","MYOD1","MYOG","PAX7","MYF5","CKM","MYMK","MYMX","MYL1",
                "MYH2","DES","MYF6","MYH4","TNNC2","DMD","MYH8","MYLPF","TTN",
                "COL1A1","COL3A1","SCX","MKX","TNC","TNMD","THBS4","COMP","FMOD","DCN","EGR1"),
xlim = c(-15, 10), ylim = c(0, 20),
pCutoff = 0.005, FCcutoff = 2, pointSize = 2.0, labSize = 3.0,
legendLabels=c('Not sig.','Log (base 2) FC','p-value','p-value & Log (base 2) FC'),
legendPosition = 'top',
legendLabSize = 10,
legendIconSize = 5.0,
drawConnectors = TRUE,
subtitle = "Differential expression",
title = 'B versus E')
dev.off()
# The End.
