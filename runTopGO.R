library(topGO)
setwd("~/topgo/")
# ---------------------------------Step 1 -----------------------------------------------------------
# The two input files: geneID2GO which is the mapping file, and the differential expressed gene list
# file: up regulated or down regulated genes
# geneID2GO file format:
# ENSG00000000003 GO:0005515,GO:0005887,GO:0039532
# ENSG00000000005 GO:0001886,GO:0001937,GO:0005515,GO:0005635,GO:0005737
# ENSG00000000419 GO:0004169,GO:0004582,GO:0005515,GO:0005634,GO:0005783,GO:0005789,GO:0006506
# ENSG00000000457 GO:0000139,GO:0005515,GO:0005524,GO:0005737,GO:0005794,GO:0006468,GO:0006954
# ENSG00000000938 GO:0001784,GO:0002768,GO:0004713,GO:0004715,GO:0005102,GO:0005515,GO:0005524
# The following are the up or down regulated gene file format:
# ENSG00000164161
# ENSG00000164283
# ENSG00000205426
# ENSG00000159167
# ---------------------------------Step 2 -----------------------------------------------------------
# Prepare the mapping file, the data type of geneID2GO is a list
geneID2GO <- readMappings("gene2go.map")
# Prepare the differential expressed gene list, myInterestingGenes is a named factor
differentialExpressedGenes <- read.table("sample19L095_and_sample19L094.up")
myInterestingGenes <- differentialExpressedGenes[[1]]
names(myInterestingGenes) <- myInterestingGenes
# Prepare the gene universe, geneList is a named factor.
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
# ---------------------------------Step 3 -----------------------------------------------------------
ontologyTerm <- c("BP","MF","CC");
test_GO_enrichment <- function(ontologyTerm){
        # Build the topGOdata object
        GOdata<- new("topGOdata", nodeSize = 5, ontology = ontologyTerm ,allGenes = geneList,
        annot = annFUN.gene2GO,gene2GO = geneID2GO)
        # Run the test
        result <- runTest(GOdata,algorithum="weight", statistic = "fisher")
        # The ‘GenTable’ function generates a summary of the results of the enrichment analysis.
        GOTerm_enrich <- GenTable(GOdata,pva_wei=result,numChar=10000,topNodes=20)
        return (GOTerm_enrich)
}
# GOTerm_top is a list containing BP, MF, CC and the related P values
GOTerm_top <- lapply(ontologyTerm,test_GO_enrichment)
# ---------------------------------Step 4 -----------------------------------------------------------
BP_color <- rep("red",20) # the topNodes=20
BP_Term <- GOTerm_top[[1]]$Term
BP_Pvalue<- gsub("<","",GOTerm_top[[1]]$pva_wei)
BP_Pvalue <- -log10(as.numeric(BP_Pvalue))
names(BP_Pvalue) <- BP_Term
BP_Pvalue <- sort(BP_Pvalue)

MF_color <- rep("green",20)
MF_Term <- GOTerm_top[[2]]$Term
MF_Pvalue<- gsub("<","",GOTerm_top[[2]]$pva_wei)
MF_Pvalue <- -log10(as.numeric(MF_Pvalue))
names(MF_Pvalue) <- MF_Term
MF_Pvalue <- sort(MF_Pvalue)

CC_color <- rep("blue",20)
CC_Term <- GOTerm_top[[3]]$Term
CC_Pvalue<- gsub("<","",GOTerm_top[[3]]$pva_wei)
CC_Pvalue <- -log10(as.numeric(CC_Pvalue))
names(CC_Pvalue) <- CC_Term
CC_Pvalue <- sort(CC_Pvalue)
data_plot <- c(BP_Pvalue,MF_Pvalue,CC_Pvalue)

# -------------------------------- Plot the bar plot -------------------------------------------------
pdf("95_94_up.pdf")
par(mar=c(1,22,3,1))
barplot(data_plot,col=c(BP_color,MF_color,CC_color),horiz=TRUE,border=NA,
        names.arg=names(data_plot),las=2,bg="white",axes=FALSE,space =2,cex.names=.5)
axis(3,cex.axis=0.5,mgp=c(2,0.5,0),tck=-0.02,tcl=0.2)
mtext(substitute(paste(-Log[10],italic(P),sep=' ')),side=3,line=1,cex=0.6)
title(main="95 vs 94 up",cex.main=0.7,font.main= 1,line=2)
legend("bottomright",legend=ontologyTerm,fill=c("red","green","blue"),bty="n",cex=.6,border=F)
dev.off()
# -------------------------------- The End -----------------------------------------------------------
