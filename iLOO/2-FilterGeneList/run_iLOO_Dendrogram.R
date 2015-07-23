DFtable <- read.table("Gene_Matrix_RSEM-counts.DF.txt", header=TRUE)
Rtable <- read.table("Gene_Matrix_RSEM-counts.R.txt", header=TRUE)
R <- data.matrix(Rtable)
DF <- data.matrix(DFtable)
DF.iLOO = iLOO(DF)
R.iLOO = iLOO(R)
results.iLOO = cbind(DF.iLOO,R.iLOO)
good.gene.inxs = which(is.nan(rowMeans(results.iLOO,na.rm=T)))

colnames(R) = as.vector(sapply(colnames(R),function(x) paste("R",sub(".Sample","",strsplit(x,"_")[[1]][2]))))
colnames(DF) = as.vector(sapply(colnames(DF),function(x) paste("DF",sub(".Sample","",strsplit(x,"_")[[1]][2]))))

all.data = cbind(R,DF)
#colnames(all.data) = c(rep("R",1,ncol(R)),rep("DF",1,ncol(DF)))
hc = hclust(dist(t(all.data[good.gene.inxs,])))
# very simple dendrogram
plot(hc,main="Sample Distance Before iLOO Filter", xlab="Sample", sub = "")


hc = hclust(dist(t(all.data)))
# very simple dendrogram
plot(hc,main="Sample Distance After iLOO Filter", xlab="Sample", sub = "")
