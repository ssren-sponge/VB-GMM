rm(list=ls())


###preliminaries
adj<-1e-15
set.seed(3)
dataDir <- "./data"
fileDir <- "./data"
plotDir <- "./outputs"
outputDir <- "./outputs"
archiveSubdir<-"GSE161529_RAW"

colIdx<-c("chartreuse4","skyblue1","tomato3","goldenrod3","slateblue3","sienna1","thistle1","turquoise3","violet","violetred","red","bisque3","aquamarine3","darkgrey","darkolivegreen3","darkseagreen3","navajowhite3","rosybrown3","dimgrey","darkorchid3","deepskyblue2","skyblue3","paleturquoise1","saddlebrown","tan3");#defining good colours for plotting

### load packages and functions
CRANpackages<-c("Matrix","rARPACK","flexclust","mclust","PMA","uwot","parallel","BiocManager")
for(iP in 1:length(CRANpackages)){
  tempPackage<-CRANpackages[iP] 
  pckg = try(require(tempPackage,character.only=T))
  if(!pckg){
    print(paste0("Installing '",tempPackage,"' from CRAN"))
    install.packages(tempPackage,repos="http://cran.r-project.org")
    require(tempPackage,character.only=T)	
  }
}
numCores<-detectCores()-2
BIOCpackages<-c("biomaRt","sva","SingleCellExperiment","scDblFinder","BiocParallel")
for(iP in 1:length(BIOCpackages)){
  tempPackage<-BIOCpackages[iP]
  pckg = try(require(tempPackage,character.only=T))
  if(!pckg){
    print(paste0("Installing '",tempPackage,"' from bioconductor"))
    BiocManager::install(tempPackage)
    require(tempPackage,character.only=T)
  }
}

ensembl<-useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl");#useful database with info about all the genes
geneLkp<-as.matrix(getBM(attributes=c("ensembl_gene_id","external_gene_name","chromosome_name","protein_id"),mart=ensembl))

MancMarkers<-as.matrix(read.csv(paste0(dataDir,"/MancBrMarkers.csv"),as.is=T,stringsAsFactors=F))
MancMarkers[MancMarkers=="Luminal_HR_pos"]<-"Luminal mature"
MancMarkers[MancMarkers=="luminal_HR_neg"]<-"Luminal progenitor"
# MancMarkers<-MancMarkers[MancMarkers[,"Type"]%in%c("Luminal progenitor","Luminal mature","Basal"),]. #!!!
MancTypes<-unique(MancMarkers[,"Type"])
MancCols<-colIdx[1:length(MancTypes)]
names(MancCols)<-MancTypes

source(paste0(dataDir,"/scatterColourAssign.R"))

### load data
load(paste0(dataDir,"/PalEApdata.rd"))
allGenes<-read.delim(paste0(dataDir,"/features.tsv"),header=F,stringsAsFactors=F)[,"V2"];#GSE161529
fileNames<-names(which(table(gsub("-matrix.mtx.gz","",gsub("-barcodes.tsv.gz","",dir(paste0(dataDir,"/",archiveSubdir)),fixed=T),fixed=T))==2))
batchPhen<-phenData[(phenData[,"characteristics_ch1.4"]%in%c("cancer type: Normal","cancer type: BRCA1 pre-neoplastic"))&(phenData[,"characteristics_ch1.2"]!="menopause: Post")&(phenData[,"characteristics_ch1.3"]%in%c("parity: Nulliparous"))&(phenData[,"characteristics_ch1.5"]%in%c("cell population: Total")),]
GSMs<-rownames(batchPhen)
print(phenData[GSMs,c("characteristics_ch1.2","characteristics_ch1.3","characteristics_ch1.4","characteristics_ch1.5")])

exprList<-sapply(GSMs,function(tempGSM){
  print(paste0(match(tempGSM,GSMs),"/",length(GSMs)))
  tempMat<-readMM(file=paste0(dataDir,"/",archiveSubdir,"/",tail(strsplit(batchPhen[tempGSM,"supplementary_file_2"],split="/",fixed=T)[[1]],1)))
  colnames(tempMat)<-read.delim(paste0(dataDir,"/",archiveSubdir,"/",tail(strsplit(batchPhen[tempGSM,"supplementary_file_1"],split="/",fixed=T)[[1]],1)),header=F,stringsAsFactors=F)[,"V1"]
  rownames(tempMat)<-allGenes
  return(tempMat)
},simplify=F,USE.NAMES=T)
batchNames<-gsub("patient: ","",batchPhen[GSMs,"characteristics_ch1"],fixed=T)
batch<-do.call(c,lapply(1:length(exprList),function(tempI){return(rep(batchNames[[tempI]],ncol(exprList[[tempI]])))}))
exprData<-do.call(cbind,exprList)
rm(list="exprList")
batchVals<-unique(batch)
batchCols<-colIdx[1:length(batchVals)]
names(batchCols)<-batchVals

toKeepR<-rowSums(exprData>0L)>=3
toKeepC<-colSums(exprData>0L)>=200
batch<-batch[toKeepC]
exprData<-exprData[toKeepR,toKeepC];#Seurat settings
numUMI<-colSums(exprData)
numTranscr<-colSums(exprData>0L)
allGenes<-rownames(exprData)
allChr<-geneLkp[match(allGenes,geneLkp[,"external_gene_name"]),"chromosome_name"];#find the chromosome for each gene
fracMitoRds<-colSums(exprData[grepl("^MT-",allGenes),])/colSums(exprData);#if a large proportion of the reads for a sample (i.e. column) are from mitochondia (chromosome = MT) then this indicates that the mitochondria in the sample have burst during the experiment, and therefore that the data may be contaminated
fracRiboReads<-colSums(exprData[grepl("^RP[LS]",allGenes),])/colSums(exprData)
fracLargestCount<-apply(exprData[rownames(exprData)!="MALAT1",],2,max)/numUMI
rm(list=ls()[grepl("toKeep",ls(),fixed=T)])

##Doublets remove
#convert exprData to sce object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(exprData)))
set.seed(3)
bp <- MulticoreParam(3, RNGseed=1234)
set.seed(3) #
sce_Dbl_para <- scDblFinder(sce, dbr=0.006*(ncol(exprData)/1000), BPPARAM=bp)
tempDblsc<-sce_Dbl_para@colData@listData[["scDblFinder.score"]]
allID <- colnames(exprData)
Dblnames <- colnames(exprData)[tempDblsc<=0.999]
# save(list = "allID","Dblnames", file = "~/Desktop/Dbl_BC1_ID.rd")
Dblsc<-tempDblsc<=0.999

toKeepP<-!(grepl("^MT-",allGenes)|grepl("^RP[LS]",allGenes))
toKeepG<-!(allGenes%in%c("5S_rRNA","7SK","Metazoa_SRP","uc_338","Y_RNA",paste0("U",1:9)))
toKeepCC<-rowSums(exprData>0L)>=50;

toKeepM<-fracMitoRds<0.1
toKeepRB<-fracRiboReads<0.5
toKeepGN<-colSums(exprData>0L)>=750
toKeepGS<-fracLargestCount<0.1

toKeepR<-toKeepP&toKeepG&toKeepCC
toKeepC<-toKeepM&toKeepRB&toKeepGN&toKeepGS&Dblsc
allGenes<-allGenes[toKeepR]
allChr<-allChr[toKeepR]
batch<-batch[toKeepC]
numUMI<-numUMI[toKeepC]
numTranscr<-numTranscr[toKeepC]
fracMitoRds<-fracMitoRds[toKeepC]
fracRiboReads<-fracRiboReads[toKeepC]
fracLargestCount<-fracLargestCount[toKeepC]
exprData<-exprData[toKeepR,toKeepC]
rm(list=ls()[grepl("toKeep",ls(),fixed=T)])
rawCounts<-exprData
logData<-log10(exprData+1)
rm(list="exprData")
logData<-ComBat(logData,batch=batch);#batch correction, to remove differences due to different labs or different individuals who donate data
logData<-logData-min(logData)+adj
exprData<-(10^logData)-1
rm(list="logData")

dupGenes<-names(which(table(allGenes)>1))
if(length(dupGenes)>0){;#aggregate data for genes that appear multiple times
  dupData<-t(sapply(dupGenes,function(tempGene){
    return(apply(exprData[allGenes%in%tempGene,],2,max))
  }))
  toKeep<-!allGenes%in%dupGenes
  exprData<-exprData[toKeep,]
  exprData<-rbind(exprData,dupData)
  rm(list="dupData")
}
allGenes<-rownames(exprData)

###marker gene processing
logData<-t(scale(t(log10(as.matrix(exprData)+1)),center=F,scale=T))
typeMeans<-t(sapply(MancTypes,function(tempType){
  return(colMeans(logData[intersect(rownames(exprData),MancMarkers[MancMarkers[,"Type"]%in%tempType,"Gene"]),,drop=F]))
},simplify=T,USE.NAMES=T))
rm(list="logData")
plotCols2<-sapply(MancTypes,function(tempType){
  return(scatterColourAssign(typeMeans[tempType,],"black",MancCols[[tempType]],100))
},simplify=F,USE.NAMES=T)

##dimension reduction
tempTimeTS<-system.time({
  topGenes<-names(sort(log10(apply(exprData,1,var))/log10(apply(exprData,1,mean)),decreasing=T))[1:2000]
  tempData<-scale(t(log10(exprData[topGenes,]+1)),center=T,scale=F)
})[[3]]
print(paste0("Transform and scale: t=",signif(tempTimeTS,digits=3),"s"))
set.seed(3)
tempTimePC<-system.time({SPCAdata<-SPC(tempData,sumabsv=0.95*sqrt(ncol(tempData)),niter=50,K=25,trace=F,center=F,compute.pve=F)})[[3]]
print(paste0("SPCA: t=",signif(tempTimePC,digits=3),"s"))
rm(list="tempData")
tempVars<-SPCAdata[["d"]][-1]^2
tempVars<-tempVars/sum(tempVars)
SPCAdata<-SPCAdata[["u"]][,-1,drop=F]
sigDim<-max(which(tempVars>0.025))
SPCAdata<-SPCAdata[,1:sigDim,drop=F]

##project to 2D
set.seed(3)
tempTimeUM<-system.time({tempProjSPC<-uwot::umap(SPCAdata,n_neighbors=50,min_dist=0.5,metric="cosine",seed = 3,n_threads=numCores,n_sgd_threads = 1,verbose=F)})[[3]]
print(paste0("UMAP: t=",signif(tempTimeUM,digits=3),"s"))

##GMM-LE clustering
tempK<-8
adjMat<-log10(exprData[intersect(rownames(exprData),topGenes),]+1)
nR<-nrow(adjMat);
nC<-ncol(adjMat);
degDistR<-rowSums(adjMat);
degDistC<-colSums(adjMat);
DR_2<-((degDistR+max(1,median(degDistR)))^-(1/2));
DC_2<-((degDistC+max(1,median(degDistC)))^-(1/2));
L<-t(t(DR_2*adjMat)*DC_2)
rm(list=c("DR_2","DC_2"))
set.seed(3)
tempSVD<-svds(L,k=100,nu=0,nv=0)
Lvars<-tempSVD[["d"]][-1]^2
Lvars<-round(Lvars/sum(Lvars),digits=3)
maxDim<-max(which(Lvars>=0.03)) 
plot(Lvars,pch=20,col=c(rep("red",maxDim),rep("black",100-maxDim)))
tempSVD<-svds(L,k=(maxDim+1),nu=(maxDim+1),nv=(maxDim+1))
rm(list="L")
V<-tempSVD[["v"]][,2:(maxDim+1)]
leveragesC<-sqrt(rowSums(V^2));
toKMclusterC<-which(leveragesC>=2/sqrt(nC));
notKMclusterC<-setdiff(1:nC,toKMclusterC);
V<-V/sqrt(rowSums(V^2)+adj);

DCSBMclustersTempC<-rep(0,nC);
tempTimeC<-system.time(tempGMMc<-Mclust(V[toKMclusterC,,drop=F],tempK))[[3]]
DCSBMclustersTempC[toKMclusterC]<-apply(tempGMMc[["z"]],1,which.max)
DCSBMclustersTempC[notKMclusterC]<-apply(predict(tempGMMc,V[notKMclusterC,,drop=F])[["z"]],1,which.max)
range(DCSBMclustersTempC)

#Save V for VI in Python
write.csv(V, paste0(outputDir, '/V_BRCA1.csv'))


###LE projection to 2D
set.seed(3)
tempTimeUM<-system.time({tempProjLE<-uwot::umap(V,n_neighbors=50,min_dist=0.5,metric="cosine",seed = 3,n_threads=numCores,n_sgd_threads = 1,verbose=F)})[[3]]
print(paste0("UMAP LE: t=",signif(tempTimeUM,digits=3),"s"))

###GMM-LE clusters
png(file=paste0(plotDir,"/Figure1.png"),height=4,width=4.2,bg="white",res=300,units="in");{  
  
  par(mar=c(1,1,1,1),oma=c(0,0,0,2))
  
  tempCols<-rep("black",nrow(V))
  for(tempIC in 1:tempK){
    tempCols[DCSBMclustersTempC%in%tempIC]<-colIdx[tempIC]
  }	
  plot(tempProjLE[,1],tempProjLE[,2],pch=20,cex=0.8,xlab="",ylab="",main="",cex.axis=1.5,bty="n",xaxt="n",yaxt="n",col=tempCols)
  box(lwd=2)
  par(xpd=NA)
  legend(max(tempProjLE[,1]),min(tempProjLE[,2]),legend=paste(sort(unique(DCSBMclustersTempC))),pch=20,cex=0.7,col=colIdx[sort(unique(DCSBMclustersTempC))],title="Cluster",xjust=-0.25,yjust=0.05)	  
}
dev.off()


##marker plots in LE
png(file=paste0(plotDir,"/Figure3.png"),height=8,width=16,bg="white",res=300,units="in");{
  par(mar=c(3,3,3,3),mfrow=c(2,4))
  invisible(lapply(MancTypes,function(tempType){
    plot(tempProjLE[,1],tempProjLE[,2],pch=20,col=plotCols2[[tempType]],cex=0.8,xlab="",ylab="",main="",cex.axis=1.5,bty="n",xaxt="n",yaxt="n")
    box(lwd=2)
    mtext(tempType,side=3,line=0.5,cex=1.5)
  }))
}
dev.off()

#VB-GMM plot
Bay_probs <- read.csv(paste0(outputDir,"/probs_BRCA_dbl_k8.csv"))
Bay_probs <- t(Bay_probs)

rownames(Bay_probs) <- c(1:nrow(Bay_probs))  
typeCols <- colIdx[1:nrow(Bay_probs)] #select colours for use
names(typeCols) <- rownames(Bay_probs)
clusterTypes <- rownames(Bay_probs)
Bay_probs <- t(t(Bay_probs)/colSums(Bay_probs)) 

###generate colours #rownames(Bay_probs)
plotCols<-sapply(clusterTypes,function(tempType){
  return(scatterColourAssign(Bay_probs[tempType,],"black",typeCols[[tempType]], 100))
},simplify=F,USE.NAMES=T)

##plot VB-GMM clusters on LE projection
png(file=paste0(plotDir,"/Figure2.png"),height=8,width=16,bg="white",res=300,units="in");{
  par(mar=c(3,3,3,3),mfrow=c(2,4))
  invisible(lapply(clusterTypes,function(tempType){
    plot(tempProjLE[,1],tempProjLE[,2],pch=20,col=plotCols[[tempType]],cex=0.8,xlab="",ylab="",main="",cex.axis=1.5,bty="n",xaxt="n",yaxt="n")
    box(lwd=2)
    mtext(tempType,side=3,line=0.5,cex=1.5)
  }))
}
dev.off()




