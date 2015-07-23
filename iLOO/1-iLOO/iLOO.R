source("iLOO.script.download.R")

Pheno1table <- read.table("Gene_Matrix_RSEM-counts.Pheno1.txt", header=TRUE)
Pheno2table <- read.table("Gene_Matrix_RSEM-counts.Pheno2.txt", header=TRUE)

Pheno1 <- data.matrix(Pheno1table)
Pheno2 <- data.matrix(Pheno2table)

Pheno1.iLOO <- iLOO(Pheno1)
Pheno2.iLOO <- iLOO(Pheno2)

results.iLOO <- cbind(Pheno1.iLOO, Pheno2.iLOO)
