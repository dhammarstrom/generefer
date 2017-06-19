#' Analyze reference gene stability with geNorm
#'
#' @param data a data.frame() with sample id, target id and raw cq values. If duplicates, average cq values per sample.
#' @param sample column number of sample ids
#' @param detector column number of target ids
#' @param cq column number of cq values
#' @return A list from NormqPCR with genorm parameters
#' @import "qpcR"
#' @import "NormqPCR"
#' @import "dplyr"
#' @import "tidyr"
#' @export

genorm<-function(data, sample=1, detector=2, cq=3, ngenes=2){
  
  
  
  names(data)[sample]<-"Sample"
  names(data)[detector]<-"Detector"
  names(data)[cq]<-"Cq"
  
  data<-data%>%
    dplyr::select(Sample, Detector, Cq)%>%
    group_by(Sample, Detector)%>%
    summarise(Cq=mean(Cq, na.rm=T))%>%
    spread(Detector, Cq)%>%
    na.omit()%>%
    gather(Detector, Cq, -Sample)%>%
    data.frame()
  
  
  # transforms cq to linear scale
  data$relq<-rep(NA, nrow(data))
  
  
  geneList<-unique(data[,detector])
  
  
  for(i in 1:length(geneList)){
    data[data$Detector==geneList[i],]$relq<-(1/(2^data[data$Detector==geneList[i],]$Cq))/max((1/(2^data[data$Detector==geneList[i],]$Cq)))
  }
  
  data<-data[,c(1,2,4)]
  names(data)[3]<-"Cq"
  
  dir.create(file.path("temp"), showWarnings = FALSE)
  write.table(data, "./temp/temp.txt", sep = "\t")
  data.batch<-read.qPCR("./temp/temp.txt")
  
  geNorm.results<-selectHKs(data.batch, method="geNorm", Symbols=featureNames(data.batch), minNrHKs=ngenes, log=FALSE, trace=FALSE)
  
  return(geNorm.results)
  
}
