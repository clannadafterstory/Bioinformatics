#20190102 xinzeshi#
#TMB 修改#
mut699=read.delim("/PUMC/94\ ESCC/TMB/699_allmutation_input.txt",header = T)
mut699[is.na(mut699)] <- 0
rownames(mut699)=mut699[,1]
mut699=mut699[,-1]
mut699_t=t(mut699)
mut699_t=as.data.frame(mut699_t)
linlab=as.data.frame(mut699_t)
notlinlab=as.data.frame(mut704_t[grep("_WGS_",rownames(mut704_t),invert=TRUE),])
nottcga=as.data.frame(notlinlab[grep("TCGA",rownames(notlinlab),invert=TRUE),])

#根据基因是否突变分组，两组case number大于等于3时进行Wilcox rank sum test，否则输出NA
#计算p value，median，median difference，standard deviation
linlab_tmb<-data.frame(NULL)
for(i in 2:ncol(linlab))
{
  v1<-unlist(linlab[linlab[,i]==1,1])
  v2<-unlist(linlab[linlab[,i]==0,1])
  linlab_tmb[3,i-1]<- length(v1)
  linlab_tmb[4,i-1]<- length(v2)
  
  if (length(v1) >= 3 & length(v2) >=3){
    linlab_tmb[1,i-1]<- wilcox.test(v1, v2)$p.value
    linlab_tmb[2,i-1]<- median(v1)-median(v2)
    linlab_tmb[5,i-1]<- sd(v1)
    linlab_tmb[6,i-1]<- sd(v2)
    linlab_tmb[7,i-1]<- median(v1)
    linlab_tmb[8,i-1]<- median(v2)
    
  }else{
    linlab_tmb[1,i-1]<-"NA"
    linlab_tmb[2,i-1]<-"NA"
    linlab_tmb[5,i-1]<- "NA"
    linlab_tmb[6,i-1]<- "NA"
    linlab_tmb[7,i-1]<- "NA"
    linlab_tmb[8,i-1]<- "NA"
    
  }
}
p<-linlab_tmb[1,]
p<-p[-which(p=="NA")]
pfdr<-p.adjust( p, method = "fdr", n = length(p))
pbonferroni<-p.adjust(p,method="bonferroni",n=length(p))
pbonferroni
linlab_tmb[9,]<-pfdr
linlab_tmb<-rbind( linlab_tmb,pfdr)
rownames(linlab_tmb)<-c("mutCT.pvalue","median.diff","no.mut","no.wt","sd.mut","sd.wt","median.mut","median.wt","pfdr")
colnames(linlab_tmb)<-colnames(linlab)[2:ncol(linlab)]
linlab_tmb_t=as.data.frame(t(linlab_tmb))
write.table(notlinlab_snv_t,"/Users/chenyamei/Documents/Lin's lab_CYM/ESCC PANEL/Drug Prediction/610_snv_tmb.txt",col.names = T,row.names = T,sep = "\t",quote = F)
write.table(linlab_snv_t,"/Users/chenyamei/Documents/Lin's lab_CYM/ESCC PANEL/Drug Prediction/94_snv_tmb.txt",col.names = T,row.names = T,sep = "\t",quote = F)
write.table(linlab_tmb_t,"/PUMC/94\ ESCC/TMB/20181227_linlab_tmb_t2.txt",col.names = T,row.names = T,sep = "\t",quote = F)
