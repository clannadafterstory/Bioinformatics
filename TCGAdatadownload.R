#20190102 xinzeshi#
#从TCGA下载数据 from 简书 547可是贼帅的547#

# 安装R包
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")

# 加载R包
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(stringr)

#下面填入要下载的癌症种类
request_cancer=c("ESCA")
for (i in request_cancer) {
  cancer_type=paste("TCGA",i,sep="-")
  print(cancer_type)
  #下载临床数据
  clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")
  write.csv(clinical,file = paste(cancer_type,"clinical.csv",sep = "-"))
  
  #下载rna-seq的counts数据
  query <- GDCquery(project = cancer_type, 
                    data.category = "Transcriptome Profiling", 
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")  # 需注意“-”前后的空格
  
  GDCdownload(query, method = "api", files.per.chunk = 100)
  expdat <- GDCprepare(query = query)
  count_matrix=assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"Counts.csv",sep = "-"))
  
  #下载miRNA数据
  query <- GDCquery(project = cancer_type, 
                    data.category = "Transcriptome Profiling", 
                    data.type = "miRNA Expression Quantification", 
                    workflow.type = "BCGSC miRNA Profiling")
  
  GDCdownload(query, method = "api", files.per.chunk = 50)
  expdat <- GDCprepare(query = query)
  count_matrix=assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"miRNA.csv",sep = "-"))
  
  #下载Copy Number Variation数据
  query <- GDCquery(project = cancer_type, 
                    data.category = "Copy Number Variation", 
                    data.type = "Copy Number Segment")
  
  GDCdownload(query, method = "api", files.per.chunk = 50)
  expdat <- GDCprepare(query = query)
  count_matrix=assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"Copy-Number-Variation.csv",sep = "-"))
  
  #下载甲基化数据
  query.met <- GDCquery(project =cancer_type,
                        legacy = TRUE,
                        data.category = "DNA methylation")
  GDCdownload(query.met, method = "api", files.per.chunk = 300)
  expdat <- GDCprepare(query = query)
  count_matrix=assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"methylation.csv",sep = "-"))
}
#下载rna-seq的counts数据
projectid <- "TCGA-ESCA"
query.count <- GDCquery(project= projectid,
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - Counts") 


#附录 从TCGA下载数据 from 丁香园stringson#

# 下载数据

GDCdownload(query.count)

# 获得表达矩阵

dataAssay = GDCprepare(query.count, summarizedExperiment = F)

rownames(dataAssay) = as.character(dataAssay[,1])

# dataAssay就是矩阵了，它此时在R的环境变量里、也就是在计算机内存中。你可以在R中对它进行进一步的分析。

# 也可以用write.table或write.csv命令把它从R里保存出来到硬盘，并保存为csv的格式，就可以用excel打开了。

write.csv(dataAssay, "TCGA-matrix.csv")  # 此时，保存的文件名为“TCGA-matrix.csv”

