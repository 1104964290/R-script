# This script is used to extract the GO ids and related description from the .OBO file.
# The .OBO file can be downloaded from the website:     http://purl.obolibrary.org/obo/go.obo
setwd("~/topgo")
library(ontologyIndex)
go_obo <- get_ontology("go.obo", propagate_relationships = "is_a",extract_tags = "everything")
goID_description <- as.data.frame(go_obo$name)
write.table(goID_description,file="goID_description.txt",quote=FALSE,sep="\t")
