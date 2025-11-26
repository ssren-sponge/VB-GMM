rm(list=ls())

###preliminaries
adj<-1e-15;#in case we need a 'smallest value'
set.seed(1);

dataDir <- "./data"
fileDir <- "./data"
plotDir <- "./outputs"
outputDir <- "./outputs"

colIdx<-c("chartreuse4","slateblue3","tomato3","goldenrod3","thistle1","sienna1","skyblue1","turquoise3","violet","violetred","red","bisque3","aquamarine3","darkgrey","darkolivegreen3","darkseagreen3","navajowhite3","rosybrown3","dimgrey","darkorchid3","deepskyblue2","skyblue3","paleturquoise1","saddlebrown","tan3");#defining good colours for plotting

##load packages
CRANpackages<-c("httpuv","impute","rARPACK","flexclust","mclust","Matrix","PMA","Seurat","BiocManager")
for(iP in 1:length(CRANpackages)){
  tempPackage<-CRANpackages[iP]
  pckg = try(require(tempPackage,character.only=T))
  if(!pckg){
    print(paste0("Installing '",tempPackage,"' from CRAN"))
    install.packages(tempPackage,repos="http://cran.r-project.org")
    require(tempPackage,character.only=T)
  }
}

numCores<-5
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

GOI<-as.vector(as.matrix(read.csv(paste0(dataDir,"/excitatoryNeuronGenes.csv"),as.is=T,header=F))[,1]);#genes of interest from literature
oRGgenes<-as.vector(as.matrix(read.table(paste0(dataDir,"/oRGgenes.csv"),sep=",")));#gene definitions from kriegstein lab for 'outer radial glial' stem cells: important cells that lead to cortical neurons and also probably lead to brain cancer (GBM, 'glioblastoma multiforme')
vRGgenes<-as.vector(as.matrix(read.table(paste0(dataDir,"/vRGgenes.csv"),sep=",")));#gene definitions from kriegstein lab for 'ventricular radial glial' stem cells
markerData<-read.csv(paste0(dataDir,"/Gene_CellType_Correlations.csv"),as.is=T);#correlations from Kriegstein lab for a wider range of marker genes
geneNames<-markerData[,1];#extract gene names
geneNames<-names(which(table(geneNames)<2))
markerData<-markerData[match(geneNames,markerData[,1]),];#match to marker data
rownames(markerData)<-geneNames
markerData<-markerData[,-1]
markerData<-as.matrix(markerData)
cellTypes<-setdiff(colnames(markerData),c("MGE_interneuron","Pericyte","endothelial","RG"));#removes uninteresting cell types, add or remove RG
markerGenes<-sapply(cellTypes,function(tempType){
  return(names(which(abs(markerData[,tempType])>=0.5)))
},simplify=F,USE.NAMES=T);#define marker genes as those with Pearson correlation >0 with defined genes

source(paste0(dataDir,"/scatterColourAssign.R"))

### load data
fileNames<-read.csv(paste0(dataDir,"/Kriegstein_file_summary_9_29_20.csv"),as.is=T,stringsAsFactors=F)[,"sample"]
keepTerms<-c("V1$")
suppressWarnings(fileWeeks<-as.numeric(substr(fileNames,3,4)))
weekRange<-c(17,22)
fileNames<-fileNames[which((fileWeeks>=weekRange[1])&(fileWeeks<=weekRange[2])&sapply(fileNames,function(tempName){return(any(sapply(keepTerms,grepl,tempName)))}))]

dir.create(paste0(dataDir,"/Kriegstein2021/temp"))
exprList<-sapply(fileNames,function(tempName){
  print(paste0(match(tempName,fileNames),"/",length(fileNames)))
  tempFile<-head(sort(dir(paste0(dataDir,"/Kriegstein2021"))[sapply(dir(paste0(dataDir,"/Kriegstein2021")),function(temp){return(grepl(tempName,temp,fixed=T))})]),1)
  tempFiles<-untar(paste0(dataDir,"/Kriegstein2021/",tempFile),list=T)
  tempDir<-head(sort(unique(sapply(strsplit(tempFiles,split="/",fixed=T),head,1))),1)
  tempLoc<-getwd()
  setwd(paste0(dataDir,"/Kriegstein2021/temp"))
  untar(paste0(dataDir,"/Kriegstein2021/",tempFile))
  setwd(tempLoc)
  tempMat<-readMM(file=paste0(dataDir,"/Kriegstein2021/temp/",tempDir,"/matrix.mtx"))
  colnames(tempMat)<-read.delim(paste0(dataDir,"/Kriegstein2021/temp/",tempDir,"/barcodes.tsv"),header=F,stringsAsFactors=F)[,"V1"]
  rownames(tempMat)<-read.delim(paste0(dataDir,"/Kriegstein2021/temp/",tempFiles[grepl("enes",tempFiles,fixed=T)|grepl("eatures",tempFiles,fixed=T)]),header=F,stringsAsFactors=F)[,"V2"]
  return(tempMat)
},simplify=F,USE.NAMES=T)
batch<-do.call(c,lapply(1:length(exprList),function(tempI){return(rep(tempI,ncol(exprList[[tempI]])))}))
exprData<-do.call(cbind,exprList)
rm(list="exprList")
unlink(paste0(dataDir,"/Kriegstein2021/temp"),recursive=T)
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

hist(log(exprData["GAPDH",]+1),breaks=100,freq=F);#GADPH is a 'housekeeping gene' that is expressed similarly whatever the cell type

plot(numUMI,exprData["GAPDH",],pch=20,log="xy")
points(numUMI[batch==1],exprData["GAPDH",batch==1],col="tomato3",pch=20)
points(numUMI[batch==2],exprData["GAPDH",batch==2],col="chartreuse3",pch=20)
points(numUMI[batch==3],exprData["GAPDH",batch==3],col="royalblue3",pch=20)
points(numUMI[batch==4],exprData["GAPDH",batch==4],col="goldenrod3",pch=20)
points(numUMI[batch==5],exprData["GAPDH",batch==5],col="slateblue3",pch=20)
points(numUMI[batch==6],exprData["GAPDH",batch==6],col="thistle1",pch=20)
points(numUMI[batch==7],exprData["GAPDH",batch==7],col="sienna1",pch=20)

##Doublets remove
sce <- SingleCellExperiment(assays = list(counts = as.matrix(exprData))) #convert exprData to sce object
set.seed(3)
bp <- MulticoreParam(3, RNGseed=1234)
set.seed(3) 
sce_Dbl_para <- scDblFinder(sce, dbr=0.006*(ncol(exprData)/1000), BPPARAM=bp)
tempDblsc<-sce_Dbl_para@colData@listData[["scDblFinder.score"]]
allID <- colnames(exprData)
Dblnames <- colnames(exprData)[tempDblsc<=0.999]
# save(list = "allID","Dblnames", file = paste0(outputDir,"/Dbl_KS_gw1722ID.rd"))
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

logData<-log10(exprData+1)
rm(list="exprData")
logData<-ComBat(logData,batch=batch);#batch correction, to remove differences due to different labs or different individuals who donate data
logData<-logData-min(logData)+adj
exprData<-(10^logData)-1
rm(list="logData")
dim(exprData)

dupGenes<-names(which(table(allGenes)>1))
if(length(dupGenes)>0){;#aggregate data for genes that appear multiple times
  dupData<-t(sapply(dupGenes,function(tempGene){
    return(apply(exprData[allGenes%in%tempGene,],2,mean))
  }))
  toKeep<-!allGenes%in%dupGenes
  exprData<-exprData[toKeep,]
  exprData<-rbind(exprData,dupData)
  rm(list="dupData")
}
allGenes<-rownames(exprData)

###marker gene processing
logData<-t(scale(t(log10(as.matrix(exprData)+1)),center=F,scale=T))
# G2MmarkerScore<-colMeans(logData[intersect(rownames(logData),G2Mphase),]);#define 'scores' for definining colours for plots. These scores are based on gene expression levels
# SmarkerScore<-colMeans(logData[intersect(rownames(logData),Sphase),]);#define 'scores' for definining colours for plots. These scores are based on gene expression levels
oRGmarkerScore<-colMeans(logData[intersect(rownames(logData),oRGgenes),]);#define 'scores' for definining colours for plots. These scores are based on gene expression levels
vRGmarkerScore<-colMeans(logData[intersect(rownames(logData),vRGgenes),]);#define 'scores' for definining colours for plots.
markerGenesScore<-sapply(markerGenes,function(tempMarkers){#define 'scores' for assigning colours to samples for plots.
  return(colMeans(logData[intersect(rownames(logData),tempMarkers),,drop=F]))
},simplify=F,USE.NAMES=T)
rm(list="logData")
markerGenesScore[["oRG"]]<-oRGmarkerScore;#combine all scores into one object 'markerGenesScores'
markerGenesScore[["vRG"]]<-vRGmarkerScore;#combine all scores into one object 'markerGenesScores'
# markerGenesScore[["G2M"]]<-G2MmarkerScore;
rm(list=c("oRGmarkerScore","vRGmarkerScore"))
names(markerGenesScore)<-gsub("CGE_","",names(markerGenesScore),fixed=T)
markerGenesScore<-markerGenesScore[c("oRG","vRG","astrocytes","Oligo","Microglia","IPC","neuron","interneuron")] #original
#markerGenesScore<-markerGenesScore[c("RG","oRG","vRG","astrocytes","Oligo","Microglia","IPC","neuron","interneuron")]
#markerGenesScore<-markerGenesScore[c("oRG","vRG","G2M","Oligo","Microglia","IPC","neuron","interneuron")]
cellTypes<-names(markerGenesScore)
typeMeans<-do.call(rbind,markerGenesScore);#arrange as matrix
typeMeans[is.nan(typeMeans)]<-0
typeMeans[is.na(typeMeans)]<-0
rm(list="markerGenesScore")

cellCols<-colIdx[1:length(cellTypes)];#select colours for use
names(cellCols)<-cellTypes

GOI<-intersect(GOI,rownames(exprData))
GOIcols<-rep(colIdx,ceiling(length(GOI)/length(colIdx)))[1:length(GOI)]
names(GOIcols)<-GOI

plotCols<-sapply(cellTypes,function(tempType){
  return(scatterColourAssign(typeMeans[tempType,],"black",cellCols[[tempType]],100))
},simplify=F,USE.NAMES=T)

pythonColKey<-sapply(cellTypes,function(tempType){
  return(lapply(colorRampPalette(c("black", cellCols[[tempType]]))(100),col2rgb))
},simplify=F,USE.NAMES=T)
dim(exprData)

##dimension reduction
set.seed(3)
tempTimeTS<-system.time({
  topGenes<-names(sort(log10(apply(exprData,1,var))/log10(apply(exprData,1,mean)),decreasing=T))[1:2000]
  tempData<-scale(t(log10(as.matrix(exprData[topGenes,])+1)),center=T,scale=F)
})[[3]]
print(paste0("Transform and scale: t=",signif(tempTimeTS,digits=3),"s"))
tempTimePC<-system.time({SPCAdata<-SPC(tempData,sumabsv=0.95*sqrt(ncol(tempData)),niter=50,K=25,trace=F,center=F,compute.pve=F)})[[3]]
print(paste0("SPCA: t=",signif(tempTimePC,digits=3),"s"))
rm(list="tempData")
tempVars<-SPCAdata[["d"]][-1]^2
tempVars<-tempVars/sum(tempVars)
SPCAdata<-SPCAdata[["u"]][,-1,drop=F]

sigDim<-max(which(tempVars>=0.02))
SPCAdata<-SPCAdata[,1:max(10,sigDim),drop=F]

##project to 2D
set.seed(3)
tempTimeUM<-system.time({tempProjSPC<-uwot::umap(SPCAdata,n_neighbors=15,min_dist=0.1,metric="cosine",n_threads=numCores,verbose=F)})[[3]]
print(paste0("UMAP: t=",signif(tempTimeUM,digits=3),"s"))

##GMM-LE clustering
#tempK<-length(cellTypes)
tempK<-14
adjMat<-log10(exprData[intersect(rownames(exprData),topGenes),]+1) ##Simulation Data
nR<-nrow(adjMat);
nC<-ncol(adjMat);
degDistR<-rowSums(adjMat);
degDistC<-colSums(adjMat);
DR_2<-((degDistR+max(1,median(degDistR)))^-(1/2));
DC_2<-((degDistC+max(1,median(degDistC)))^-(1/2));
L<-t(t(DR_2*adjMat)*DC_2)
rm(list=c("DR_2","DC_2"))
# DC<-((degDistC+max(1,median(degDistC)))^(-1));
# L<-t(t(adjMat)*DC)
# rm(list="DC")
set.seed(3)
tempSVD<-svds(L,k=100,nu=0,nv=0)
Lvars<-tempSVD[["d"]][-1]^2
Lvars<-round(Lvars/sum(Lvars),digits=3)

# maxDim<-max(max(which(Lvars>=0.01)),tempK-1)
maxDim<-max(which(Lvars>=0.009))
plot(Lvars,pch=20,col=c(rep("red",maxDim),rep("black",100-maxDim)))
text(14, 0.017, maxDim)
tempSVD<-svds(L,k=(maxDim+1),nu=(maxDim+1),nv=(maxDim+1))
# rm(list="L") !!!
V<-tempSVD[["v"]][,2:(maxDim+1)]
# V<-t(tempSVD[["d"]][2:(maxDim+1)]*t(tempSVD[["v"]][,2:(maxDim+1)]))
leveragesC<-sqrt(rowSums(V^2));
toKMclusterC<-which(leveragesC>=2/sqrt(nC));
# toKMclusterC<-1:nC
notKMclusterC<-setdiff(1:nC,toKMclusterC);
V<-V/sqrt(rowSums(V^2)+adj);

# Vsc<-t(tempSVD[["d"]][2:(maxDim+1)]*t(tempSVD[["v"]][,2:(maxDim+1)]))
# Vsc<-Vsc/sqrt(rowSums(Vsc^2)+adj)

DCSBMclustersTempC<-rep(0,nC);
tempTimeC<-system.time(tempGMMc<-Mclust(V[toKMclusterC,,drop=F],tempK))[[3]]
DCSBMclustersTempC[toKMclusterC]<-apply(tempGMMc[["z"]],1,which.max)
DCSBMclustersTempC[notKMclusterC]<-apply(predict(tempGMMc,V[notKMclusterC,,drop=F])[["z"]],1,which.max)
range(DCSBMclustersTempC)

#Save V for VB-GMM
write.csv(V,paste0(outputDir, '/V_GW1722_D16.csv'))

#Save label & markers for AUC
# write.csv(DCSBMclustersTempC, paste0(outputDir,'/gw1722_k', tempK, '_GMM.csv'), row.names = FALSE)
# save(list = 'typeMeans', file = paste0(outputDir,'/GW1722_typeMeans.rd'))

###LE projection to 2D#--------------------------------------------------------------------
set.seed(3)
tempTimeUM<-system.time({tempProjLE<-uwot::umap(V,n_neighbors=15,min_dist=0.5,metric="cosine",n_threads=numCores,verbose=F)})[[3]]
print(paste0("UMAP LE: t=",signif(tempTimeUM,digits=3),"s"))

save(list = 'tempProjLE', file = paste0(outputDir,'/GW1722_tempProjLE.rd'))

##plot marker genes LE
png(file=paste0(plotDir,"/Figure5.png"),height=16,width=8,bg="white",res=300,units="in");{
  par(mar=c(3,3,3,3),mfrow=c(4,2))
  invisible(lapply(cellTypes,function(tempType){
    plot(tempProjLE[,1],tempProjLE[,2],pch=20,col=plotCols[[tempType]],cex=0.8,xlab="",ylab="",main="",cex.axis=1.5,bty="n",xaxt="n",yaxt="n")
    box(lwd=2)
    mtext(tempType,side=3,line=0.5,cex=1.5)
  }))
}
dev.off()


#plot GMM-LE clusters in LE projection
png(file=paste0(plotDir,"/Figure6.png"),height=4,width=4.2,bg="white",res=300,units="in");{  
  
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

###Seurat algorithm
expD <- CreateAssayObject(exprData)
srt<-CreateSeuratObject(counts=expD,project="KSneuro10X");#use 'Seurat' package, for comparison
srt<-NormalizeData(srt);#standard pipeline we compare our results with: normalise data,
srt<-FindVariableFeatures(srt,selection.method="vst",nfeatures=2000);#then find the variable features
all.genes<-rownames(srt)
srt<-ScaleData(srt,features=all.genes);#scale the data to N(0,1)
srt<-RunPCA(srt,features=VariableFeatures(object=srt));#Find top PCA components
srt<-FindNeighbors(srt,dims=1:25);#Find neighbours for kNN clustering
srt<-FindClusters(srt,resolution=0.325)
srt<-RunUMAP(srt,dims=1:25)
tempClust<-as.numeric(srt@active.ident);#extract cluster markers
names(tempClust)<-colnames(exprData)
#plot by Louvain/srt
png(file=paste0(plotDir,"/Figure4.png"),height=4,width=4,bg="white",res=300,units="in");{
  par(mar=c(3,3,3,3))
  tempCols<-rep("black",length(tempClust))
  for(tempI in unique(tempClust)){
    tempCols[tempClust%in%tempI]<-colIdx[tempI]
  }
  plot(tempProjLE[,1],tempProjLE[,2],pch=20,col=tempCols,cex=0.8,xlab="",ylab="",main="",cex.axis=1.5,bty="n",xaxt="n",yaxt="n")
  box(lwd=2)
  par(xpd=NA)
  legend(max(tempProjLE[,1]),min(tempProjLE[,2]),legend=paste(sort(unique(tempClust))),pch=20,cex=0.7,col=colIdx[sort(unique(DCSBMclustersTempC))],title="Cluster",xjust=-0.25,yjust=0.05)	  
  
}
dev.off()


### VB-GMM
Bay_probs <- read.csv(paste0(outputDir,"/probs_KSneuro_k14.csv")) #probs_KSneuro_k14.csv is saved from KSneuro.ipynb
Bay_probs <- t(Bay_probs)

###Plot for VI, wVI
rownames(Bay_probs) <- c(1:nrow(Bay_probs))  #change as different clusters 
typeCols <- colIdx[1:nrow(Bay_probs)] #select colours for use
names(typeCols) <- rownames(Bay_probs)
clusterTypes <- rownames(Bay_probs)
Bay_probs <- t(t(Bay_probs)/colSums(Bay_probs)) 

###generate colours #rownames(Bay_probs)
plotCols<-sapply(clusterTypes,function(tempType){
  return(scatterColourAssign(Bay_probs[tempType,],"black",typeCols[[tempType]], 100))
},simplify=F,USE.NAMES=T)

#LE
png(file=paste0(plotDir,"/Figure7.png"),height=16,width=16,bg="white",res=300,units="in");{
  par(mar=c(3,3,3,3),mfrow=c(4,4))
  invisible(lapply(clusterTypes,function(tempType){
    plot(tempProjLE[,1],tempProjLE[,2],pch=20,col=plotCols[[tempType]],cex=0.8,xlab="",ylab="",main="",cex.axis=1.5,bty="n",xaxt="n",yaxt="n")
    box(lwd=2)
    mtext(tempType,side=3,line=0.5,cex=1.5)
  }))
}
dev.off()



