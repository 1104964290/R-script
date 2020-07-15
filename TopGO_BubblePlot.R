library(topGO)
setwd("~/ShaoXieXiangRNAseq/EditedGO/")
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
geneID2GO <- readMappings("gene2go_for_human.map")
# Prepare the differential expressed gene list, myInterestingGenes is a named factor
differentialExpressedGenes <- read.table("b_e.down")
myInterestingGenes <- differentialExpressedGenes[[1]]
names(myInterestingGenes) <- myInterestingGenes
# Prepare the gene universe, geneList is a named factor.
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
# ---------------------------------Step 3 -----------------------------------------------------------
# Set nodeSize=1
GOdata<- new("topGOdata", nodeSize = 1, ontology = "BP" ,allGenes = geneList,
        annot = annFUN.gene2GO,gene2GO = geneID2GO)
# Run the test
result <- runTest(GOdata,algorithum="weight", statistic = "fisher")
# The ‘GenTable’ function generates a summary of the results of the enrichment analysis.
GOTerm_enrich <- GenTable(GOdata,pvalue=result,numChar=10000,topNodes=35)
write.table(GOTerm_enrich,file="go_b_e_down.txt",sep="\t",quote=FALSE,row.names=FALSE)
# Generate the bubble plot
options(digits = 2)
input <- read.table("go_b_e_down.txt",header=TRUE,sep="\t")
input$Rich <- input$Significant/input$Annotated
library(ggplot2)
pdf("b_e_down.pdf")
ggplot(input, aes(x=-log10(pvalue),y=Term,size=Rich,color="red"))+ geom_point(alpha=1)+guides(col = FALSE)+
scale_size(range = c(.1, 10), name="Rich")+
theme(axis.text=element_text(color = "black"),
     plot.background = element_rect(fill = "white")
)+
labs(x="-log10(pvalue)",y="Term",face="bold")+xlim(c(1,30))
dev.off()
# The End.
