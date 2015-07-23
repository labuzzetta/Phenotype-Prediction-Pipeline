DFSamples <- c(<Sample Numbers>) #INPUT REQUIRED
RSamples <- c(<Sample Numbers>) #INPUT REQUIRED

for(i in 1:length(RSamples)){
  
  for(j in 1:length(RSamples)){
    
    dir.name <- paste(<Create Pattern for directores>) #INPUT REQUIRED
    
    dir.create(dir.name)
    
    setwd(dir.name)
    
    data.dir <- ("./")
    
    ctrl.pattern <- paste(<REGEX>) #INPUT REQUIRED
    
    case.pattern <- paste(<REGEX>) #INPUT REQUIRED
    
    case.files <- dir(data.dir, pattern=case.pattern, full.names = TRUE)
    
    print("Case Files")
    
    print(case.files)
    
    control.files <- dir(data.dir, pattern=ctrl.pattern, full.names = TRUE)
    
    print("Control Files")
    
    print(control.files)
    
    geneset.file <- ("../../msigdb.v5.0.symbols.gmt.txt") #CHANGE AS NEEDED
    
    geneID.type <- "gene.symbol"
    
    output.prefix <- "SeqGSEA.L2O"
    
    nCores <- 1
    
    perm.times <- 1000
    
    DEonly <- FALSE
    
    DEweight <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    
    integrationMethod <- "linear"
    
    library(doParallel)

    cl <- makeCluster(30)

    registerDoParallel(cl)
    
    RCS <- loadExonCountData(case.files, control.files)
    
    RCS <- exonTestability(RCS, cutoff=5)
    
    geneTestable <- geneTestability(RCS)
    
    RCS <- subsetByGenes(RCS, unique(geneID(RCS))[ geneTestable ])
    
    geneIDs <- unique(geneID(RCS))
    
    RCS <- estiExonNBstat(RCS)
    
    RCS <- estiGeneNBstat(RCS)
    
    print("Generated RCS")
    
    permuteMat <- genpermuteMat(RCS, times=perm.times)
    
    print("Permuted Matrix")
    
    RCS <- DSpermute4GSEA(RCS, permuteMat)
    
    print("DSpermute4GSEA Done!")
    
    geneCounts <- getGeneCount(RCS)
    
    label <- label(RCS)
    
    DEG <-runDESeq(geneCounts, label)
    
    DEGres <- DENBStat4GSEA(DEG)
    
    print("DENBStat4GSEA Done!")
    
    DEpermNBstat <- DENBStatPermut4GSEA(DEG, permuteMat)
    
    DEscore.normFac <- normFactor(DEpermNBstat)
    
    DEscore <- scoreNormalization(DEGres$NBstat, DEscore.normFac)
    
    DEscore.perm <- scoreNormalization(DEpermNBstat, DEscore.normFac)
    
    print("DE Score Perm Done!")
    
    DSscore.normFac <- normFactor(RCS@permute_NBstat_gene)
    
    DSscore <- scoreNormalization(RCS@featureData_gene$NBstat, DSscore.normFac)
    
    DSscore.perm <- scoreNormalization(RCS@permute_NBstat_gene, DSscore.normFac)
    
    print("DS Score Perm Done!")
    
    gene.score <- geneScore(DEscore, DSscore, DEweight=0.5)
    
    gene.score.perm <- genePermuteScore(DEscore.perm, DSscore.perm, DEweight=0.5)
    
    print("Gene Score Perm Done!")
    
    #gene.set <- loadGenesets(geneset.file, geneIDs, geneID.type="gene.symbol", genesetsize.min = 5, genesetsize.max = 1000)
    
    #gene.set <- GSEnrichAnalyze(gene.set, gene.score, gene.score.perm, weighted.type=1)
    
    #print("Gene Set Enriched!")
    
    #GSEAres <- GSEAresultTable(gene.set, TRUE)
    
    #write.table(GSEAres, paste(output.prefix,".GSEA.result.txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE)
    
    #runSeqGSEA(data.dir=data.dir, case.pattern=case.pattern, ctrl.pattern=ctrl.pattern,geneset.file=geneset.file, geneID.type=geneID.type, output.prefix=output.prefix,nCores=nCores, perm.times=perm.times, integrationMethod=integrationMethod,DEonly=DEonly, DEweight=DEweight)
    
    write.table(gene.score, paste(output.prefix,".geneScores.txt",sep=""), quote=FALSE, sep="\t", row.names=TRUE)
    
    write.table(DEscore, paste(output.prefix,".DEScores.txt",sep=""), quote=FALSE, sep="\t", row.names=TRUE)
    
    write.table(DSscore, paste(output.prefix,".DSScores.txt",sep=""), quote=FALSE, sep="\t", row.names=TRUE)
    
    write.table(geneIDs, paste(output.prefix, ".geneIDs.txt",sep=""), quote=FALSE, sep="\t", row.names=TRUE)
    
    print("Tables written")
    
    setwd("..")
    
  }
  
}
