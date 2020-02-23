# using pheatmap package
# First create a matrix which contains the expression values, in this example, the expression values are the
# raw count generated from htseq-count, the following matrix artificial raw count generated in R
a <- sample(1:2000,1500,replace = TRUE)
express.matrix <- matrix(a,nrow = 100,ncol = 15,byrow = TRUE)
# add row name and column name to the matrix
rname <- c()
for (i in 1:100) {
        rname[i] <- paste("gene",i,sep = "")
}
rownames(express.matrix) <- rname

lname <- c()
for (i in 1:15) {
        lname[i] <- paste("sampe",i,sep = "")
}

colnames(express.matrix) <- lname
# transform the raw count number using log2
express.matrix <- log2(express.matrix+1)
# set colors
heatmap.color <- colorRampPalette(c("green","black","red"))(100)
# when using the kmeans_k parameter in pheatmap function, we should set the seed for reproducible result
set.seed(123)
res <- pheatmap(express.matrix,kmeans_k = 15,color = heatmap.color,border_color = NA,angle_col=315)
# extract the genes in each cluster
cluster_index <- res$kmeans$cluster # cluster_index is a vector with genes are the element name
output <- cbind(sort(res$kmeans$cluster),express.matrix[names((sort(res$kmeans$cluster))),])
# write the result output to a file
write.table(output,file = "result.txt",row.names = TRUE, sep ="\t", quote = FALSE)
