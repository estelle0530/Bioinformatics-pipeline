library(edgeR)
library(biomaRt)
library(pheatmap)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(dplyr)

#####Download preprocessed HTseq WGBS consistent RNA seq and other RNA seq 
#####Data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106838

###Get gene names mapped to mouse ID

ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
r1<-read.delim("/Users/yaoestelle/Downloads/GSM2855060_SC4-t0_not_stranded_htseq_output.txt",sep="\t",header=FALSE)
gene<-getBM(attributes = c( "ensembl_gene_id", "external_gene_name"), 
      filters = "ensembl_gene_id", 
      values = r1$V1, 
      mart = ensembl)
# which("RNA18S1"%in%gene$external_gene_name)
# View(gene$external_gene_name)

#####WGBS consistent samples
setwd("/Users/yaoestelle/Downloads/bab1")
file_list <- list.files(full.names = TRUE)

All <- lapply(file_list,function(i){read.delim(i, header=FALSE, sep='\t')})
rnaseq<-data.frame(All[1])
for (i in 2:length(file_list)){
  rnaseq=cbind(rnaseq, All[[i]]$V2)
}
colnames(rnaseq)[1]<-"Gene"
colnames<-c(rep("MEF",2) , 
            rep("D3-inefficient",2),
            rep("D3-efficient",2), 
            rep("D6-inefficient",2),
            rep("D6-efficient",2),
            rep("iPSC",2))
colnames(rnaseq)[2:13]<-colnames
rnaseq<-merge(gene,rnaseq, by.y = "Gene", by.x="ensembl_gene_id")
rnaseq$external_gene_name<-toupper(rnaseq$external_gene_name)

###selected age related genes
gene_list<-c("DNMT3A",'DNMT3B',"DNMT1",'TET1','TET2','TET3',
             'HOXA1','HOXA2','HOXA3','HOXA4',
             'HOXA5','HOXA6','HOXA7','HOXA9','HOXA11','HOXA13','LMNB1','CHAF1B','CHAF1A','STAG1','STAG2')
rnaseq_sub <-rnaseq[which(rnaseq$external_gene_name%in%gene_list),]
rownames(rnaseq_sub)<-rnaseq_sub$external_gene_name
rnaseq_sub<-as.data.frame(t(rnaseq_sub[,-(1:2)]))
rnaseq_sub$sample<-colnames

###CPM normalized by library size 
y <- DGEList(counts=rnaseq[,-(1:2)], group=colnames)
y <- calcNormFactors(y)
y$samples
logcpm <- cpm(y, log=TRUE)

wgbs_seq<-logcpm

locpm_subset <- logcpm[which(rnaseq$external_gene_name%in%gene_list),]
rna_subset<-(rnaseq[which(rnaseq$external_gene_name%in%gene_list),])
colnames(locpm_subset)=colnames
rownames(locpm_subset)<-rna_subset$external_gene_name

###heatmap for subset of the genes 
pheatmap(locpm_subset, show_colnames=T, cluster_cols = F)

sample<-t(locpm_subset)
sample<-as.data.frame(sample)
sample$sample<-colnames
sample$sample<-factor(sample$sample, levels=c("MEF", "D3-inefficient", 'D3-efficient', 'D6-inefficient','D6-efficient','iPSC' ))

###box plots of distribution per sample for just WGBS consistent samples
for (i in 1:length(gene_list)){
  #i=1
  temp = sample[,c(i,ncol(sample))]
  colnames(temp) <- c("Gene",'Sample')
  
  g<-ggplot(temp, aes(y=Gene,
                     x=Sample, 
                     fill=Sample))+
    geom_boxplot()+
    theme( axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           text = element_text(size=15),
           axis.title.x = element_text(size=15,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
           axis.title.y = element_text(size=15,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
           axis.text.x = element_text(size=15),
           axis.text.y = element_text(size=15))+
    xlab("Time points and reprogramming efficiency")+
    ylab(paste0(colnames(sample)[i]," (Normalized log2(CPM+1))"))
  
  ggsave(paste0(colnames(sample)[i],".png"),g, width=10,height=6)
  
}

wgbssample<-sample

#####Get other RNA seq samples 
setwd("/Users/yaoestelle/Downloads/GSE106835_RAW")
file_list <- list.files(full.names = TRUE)
All <- lapply(file_list,function(i){read.delim(i, header=FALSE, sep='\t')})

rnaseq<-data.frame(All[1])
for (i in 2:length(file_list)){
  rnaseq=cbind(rnaseq, All[[i]]$V2)
}
colnames(rnaseq)[1]<-"Gene"
name<- gsub("./GSM28550","",file_list)
name<- gsub("_not_stranded_htseq_output.txt","",name)
colnames(rnaseq)[2:34]<-name
rnaseq<-merge(gene,rnaseq, by.y = "Gene", by.x="ensembl_gene_id")
rnaseq$external_gene_name<-toupper(rnaseq$external_gene_name)

colnames <- c(rep("D0",3), 
              rep("D3-ineffi",3), rep("D6-ineffi",3), 
              rep("D9-inefficient",3), rep("D12-inefficient",3), 
              rep('iPSC-inefficient',3),
              rep("D3-effi",3), rep("D6-effi",3), 
              rep("D9-efficient",3), rep("D12-efficient",3),
              rep("iPSC-efficient",3))
sample_group<-c(rep("D0",3), 
                rep("D3",3), rep("D6",3), rep("D9",3), rep("D12",3),
                rep("iPSC",3),
                rep("D3",3), rep("D6",3), rep("D9",3), rep("D12",3), 
                rep('iPSC',3) )
  
###CPM normalized by library size 
y <- DGEList(counts=rnaseq[,-c(1:2)], group=colnames)
y <- calcNormFactors(y)
logcpm <- cpm(y, log=TRUE)

locpm_subset <- logcpm[which(rnaseq$external_gene_name%in%gene_list),]
colnames(locpm_subset)=colnames
rownames(locpm_subset)<-rna_subset$external_gene_name
pheatmap(locpm_subset, show_colnames=T, cluster_cols = F)

sample<-t(locpm_subset)
sample<-as.data.frame(sample)
sample$colname<-colnames
sample$sample<-sample_group

###Combine both samples only keeping D9 and D12 for the additional RNAseq data
alter<-sample[which(sample$sample=="D9"| sample$sample=="D12"),-ncol(sample)] 
colnames(alter)<-colnames(wgbssample)
bind<-rbind(wgbssample,alter)

###combine entire dataframe
t1<-logcpm[,which(grepl('D9',colnames)|grepl('D12', colnames))]
tcomb<-cbind(wgbs_seq , t1)
colnames(tcomb)<-bind$sample
rownames(tcomb)<-gene$external_gene_name
tcomb<-as.data.frame(tcomb)
tcomb<-cbind(gene=gene$external_gene_name, tcomb )
write.csv(tcomb,"RNA-seq-log2CPM.csv" )


###
bind$sample<-factor(bind$sample, levels=c("MEF", "D3-inefficient", 'D3-efficient',
                                          'D6-inefficient','D6-efficient',
                                          'D9-inefficient','D9-efficient',
                                          'D12-inefficient','D12-efficient',
                                          'iPSC'))


bind$efficiency <- ifelse(bind$sample=="D3-inefficient"|
                            bind$sample=="D6-inefficient"|
                            bind$sample=="D9-inefficient"|
                            bind$sample=="D12-inefficient", "Inefficient", 'Efficient')
table(bind$efficiency)

for (i in 1:length(gene_list)){
  
  i= which(colnames(bind)=="TET3")
  temp = bind[,c(i,ncol(bind)-1,ncol(bind))]
  colnames(temp) <- c("Gene",'Sample','Efficiency')


  my_comparisons <- list( c("D3-inefficient", 'D3-efficient'), 
                          c('D6-inefficient','D6-efficient'), 
                          c('D9-inefficient','D9-efficient'),
                          c('D12-inefficient','D12-efficient')
                          )
  
  mean_group<- temp%>%group_by(Sample)%>% dplyr::summarize(Mean = median(Gene, na.rm=TRUE))
  
  mean_group$Efficiency<-ifelse(mean_group$Sample=="D3-inefficient"|
                                  mean_group$Sample=="D6-inefficient"|
                                  mean_group$Sample=="D9-inefficient"|
                                  mean_group$Sample=="D12-inefficient", "Inefficient", 'Efficient')
  mean_group$Days<-ifelse(grepl( 'D3',mean_group$Sample) , "D3",
                          ifelse(grepl('D6',mean_group$Sample), "D6",
                                 ifelse(grepl('D9', mean_group$Sample), "D9",
                                        ifelse(grepl( 'D12',mean_group$Sample), "D12", 
                                               ifelse( grepl('MEF', mean_group$Sample), "MEF", "iPSC" )))))
                       
  mean_group$Days<-factor(mean_group$Days, levels=c("MEF", 'D3','D6','D9','D12','iPSC'))
  
  g<-ggplot(mean_group, aes(y=Mean, x=Days, 
                            group=Efficiency, color=Efficiency
                            ))+
    geom_line(size=1.5)+
    scale_colour_manual(values=c(Efficient="#000066",Inefficient="#339999"))+
    theme( axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           text = element_text(size=15),
           axis.title.x = element_text(size=15,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
           axis.title.y = element_text(size=15,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
           axis.text.x = element_text(size=15, angle=60, hjust=1),
           axis.text.y = element_text(size=15))+
    xlab("Time points and reprogramming efficiency")+
    ylab(paste0(colnames(sample)[i]," Normalized log2(CPM+1)"))
  
  ggsave(paste0(colnames(sample)[i],"_line.png"),g, width=10,height=6)
  
  
  mean_list<-c()
  for (comp in my_comparisons){
    diff= as.numeric(mean_group[which(mean_group$Sample==comp[2]),"Mean"])-as.numeric(mean_group[which(mean_group$Sample==comp[1]),"Mean"])
    mean_list= append(mean_list, diff)
  }
  mean_df<-data.frame("Days"= c("D3",'D6','D9','D12'), "Difference"=mean_list)
  mean_df$Days<-factor(mean_df$Days, levels=c("D3",'D6','D9','D12'))
  

  g<-ggplot(mean_df, aes(x=Days, y=Difference, fill=Days))+
    geom_bar(stat="identity", show.legend = F)+
    theme( axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           text = element_text(size=15),
           axis.title.x = element_text(size=20,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
           axis.title.y = element_text(size=20,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
           axis.text.x = element_text(size=20),
           axis.text.y = element_text(size=20))+
    #labs(subtitle=paste0(colnames(sample)[i], " difference between efficient and inefficient samples"))+
    ylab(paste0('Differences in median expression') )
  
  ggsave(paste0(colnames(sample)[i],"_diff.png"),g, width=10,height=6)
  
  g<-ggplot(temp, aes(y=Gene,
                      x=Sample,
                      fill=Sample))+
    geom_boxplot()+
    theme( axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           text = element_text(size=15),
           axis.title.x = element_text(size=15,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
           axis.title.y = element_text(size=15,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
           axis.text.x = element_text(size=15, angle=60, hjust=1),
           axis.text.y = element_text(size=15))+
    xlab("Time points and reprogramming efficiency")+
    ylab(paste0(colnames(sample)[i]," (Normalized log2(CPM+1))"))

  ggsave(paste0(colnames(sample)[i],".png"),g, width=10,height=6)
  
}








###############Bismark M-bias samples processing

setwd("/Users/yaoestelle")
file_list <- list.files(pattern='.M-bias.txt',full.names = TRUE)

All <- lapply(file_list,function(i){read.delim(i, header=FALSE, sep='\t')})

###For CPG methylation 
cpg<-data.frame()
for (i in 1:6){
  #cpg=cbind(cpg, All[[i]][3:153,])
  
  d<-All[[i]][3:153,]
  colnames(d)=d[1,]
  colnames(d)[4]<-"methylation"
  d<-d[-1,]
  d$sample<-i
  cpg<-rbind(cpg, d)
}

cpg$position<-as.numeric(as.character(cpg$position))
cpg$methylation<-as.numeric(as.character(cpg$methylation))
cpg$Sample<-c(rep("MEF",150), rep("D3-inefficient",150), rep("D3-efficient",150), 
              rep("D6-inefficient",150), rep("D6-efficient",150), 
              rep("iPSC",150)
              )
g1<-ggplot( cpg,aes(x=position, y=methylation, 
                group=sample, color=Sample)) +
  geom_line() +
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=15),
         axis.title.x = element_text(size=15,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=15,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=15),
         axis.text.y = element_text(size=15))+
  xlab("Position in read [bp]")+
  ylab("% Methylation")+
  ggtitle("Bismark M-bias for CpG % Methylation")
  
ggsave('m_cpg.png',g1, width=10, height=6)

###For CHG  methylation 
chg<-data.frame()
for (i in 1:6){
  #cpg=cbind(cpg, All[[i]][3:153,])
  
  d<-All[[i]][157:308,]
  colnames(d)=d[1,]
  colnames(d)[4]<-"methylation"
  d<-d[-1,]
  d$sample<-i
  chg<-rbind(chg, d)
}
chg$position<-as.numeric(as.character(chg$position))
chg$methylation<-as.numeric(as.character(chg$methylation))
chg$Sample<-c(rep("MEF",151), rep("D3-inefficient",151), rep("D3-efficient",151), 
              rep("D6-inefficient",151), rep("D6-efficient",151), 
              rep("iPSC",151))

g1<-ggplot( chg,aes(x=position, y=methylation, 
                    group=sample, color=Sample)) +
  geom_line() +
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=15),
         axis.title.x = element_text(size=15,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=15,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=15),
         axis.text.y = element_text(size=15))+
  xlab("Position in read [bp]")+
  ylab("% Methylation")+
  ggtitle("Bismark M-bias for CHG % Methylation")

ggsave('m_chg.png',g1, width=10, height=6)

###For CHH  methylation
chh<-data.frame()
for (i in 1:6){
  #cpg=cbind(cpg, All[[i]][3:153,])
  
  d<-All[[i]][311:462,]
  colnames(d)=d[1,]
  colnames(d)[4]<-"methylation"
  d<-d[-1,]
  d$sample<-i
  chh<-rbind(chh, d)
}
chh$position<-as.numeric(as.character(chh$position))
chh$methylation<-as.numeric(as.character(chh$methylation))
chh$Sample<-c(rep("MEF",151), rep("D3-inefficient",151), rep("D3-efficient",151), 
              rep("D6-inefficient",151), rep("D6-efficient",151), 
              rep("iPSC",151))

g2<-ggplot( chh,aes(x=position, y=methylation, 
                    group=sample, color=Sample)) +
  geom_line() +
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=15),
         axis.title.x = element_text(size=15,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=15,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=15),
         axis.text.y = element_text(size=15))+
  xlab("Position in read [bp]")+
  ylab("% Methylation")+
  ggtitle("Bismark M-bias for CHH % Methylation")

ggsave('m_chh.png',g2, width=10, height=6)

combine<-rbind(data.frame(cpg,Methylation_type="CPG"), 
               data.frame(chg,Methylation_type="CHG"),
               data.frame(chh,Methylation_type="CHH")
               )
combine$combine_type<-paste0(combine$Sample,"-" ,combine$Methylation_type)

g<-ggplot( combine,aes(x=position, y=methylation, 
                group=combine_type, color=Sample)) +
  geom_line(aes(linetype=Methylation_type) )+
  scale_linetype_manual(values=c("twodash", "dotted",'solid'))+
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=15),
         axis.title.x = element_text(size=15,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=15,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=15),
         axis.text.y = element_text(size=15))+
  xlab("Position in read [bp]")+
  ylab("% Methylation")+
  ggtitle("Bismark M-bias for CHH % Methylation")

ggsave('m_all.png',g, width=10, height=6)

write.csv(combine, "Bismark_M-bias.csv")
