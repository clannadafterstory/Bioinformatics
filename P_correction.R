#20190102 xinzeshi#
#对P值去NA FDR & bonferroni 校正 并贴回原数据表#

p<-data.frame(NULL)
p=readxl::read_excel("/Users/shixinze/Downloads/p.xlsx",sheet = 1, col_names = T)
p<-as.data.frame(p)
row.names(p)=p[,1]

#去NA#
p_clean<- p[complete.cases(p),]
p_clean_list<-p_clean[,2]

#计算FDR & bonferroni#
for (i in 1:nrow(p_clean))
  if(isFALSE(p_clean[i,2]=="NA")){
    p_clean[i,3]<- p.adjust( p_clean_list, method = "BH", n = length(p_clean_list))[i]
    p_clean[i,4]<- p.adjust( p_clean_list, method = "bonferroni")[i]
    
  }else{
    p_clean[i,3]<- NA
    p_clean[i,4]<- NA
  }

#回帖#
p_out<- merge(p,p_clean,by="gene",all=TRUE)
colnames(p_out)<-c("Gene","Pvalue","Pvalue","FDR","Bonferroni")

write.table(p_out,"/Users/shixinze/Downloads/p_out.txt",col.names = T,row.names = T,sep = "\t",quote = F)