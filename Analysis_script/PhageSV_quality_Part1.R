#####################
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(forcats)
library(ggpubr)



color <- c(brewer.pal(12,'Set3'))
c(color[3],color[5],color[6],color[4])

sv_colors1=c("#339933", "#336699", "#CCCC33", "#CC6633")
sv_colors2=c("#CC6666", "#999933", "#339999", "#996699")
sv_colors3 = c("#FFCCCC", "#99CCFF", "#FFCC99", "#CCCCFF", "#99CCCC", "#FFCC99", "#CCFFCC")
sv_colors4=c("#9999CC","#CC9999","#99CC99")
sv_colors5 = c("#336699","#CCCC33")
sv_colors6 = c("#CC6666","#6699CC")

#############################################################################
#Benchmarking of SV calling approaches on simulated virome-enriched datasets
#############################################################################

# Fig. 1a Benchmarking of SV calling approaches on simulated virome-enriched datasets

# (1) Benchmarking on simulated virome-enriched datasets with "Even" condition
Precision=c(91.42,81.83,84.12,69.81,89.24)
Recall=c(63.80,60.24,60.13,61.20,51.06)
F1=c(77.82,69.40,70.13,65.22,64.96)
PreFilter_E2=data.frame(Precision=Precision,Recall=Recall,F1=F1,Tool=c("PSDP","PBSV","cuteSV","svim","Sniffles"),Group='Ref-filter')
Precision=c(91.64,83.44,85.60,74.81,89.88)
Recall=c(31.26,23.34,27.04,31.03,25.31)
F1=c(46.61,36.47,41.10,43.86,39.50)
NoFilter_E2=data.frame(Precision=Precision,Recall=Recall,F1=F1,Tool=c("PSDP","PBSV","cuteSV","svim","Sniffles"),Group='No-filter')
E2_performance = rbind(PreFilter_E2,NoFilter_E2)
E2_performance$Group = factor(E2_performance$Group,levels=c("No-filter","Ref-filter"))
E2_performance$Tool = factor(E2_performance$Tool,levels=c("svim","PBSV","cuteSV","Sniffles","PSDP"))


E2_performance_F1 = E2_performance[,c("Group","Tool","F1")]
colnames(E2_performance_F1) = c("Group","Tool","Metric")
E2_performance_Recall = E2_performance[,c("Group","Tool","Recall")]
colnames(E2_performance_Recall) = c("Group","Tool","Metric")
E2_performance_Precision = E2_performance[,c("Group","Tool","Precision")]
colnames(E2_performance_Precision) = c("Group","Tool","Metric")
E2_performance_F1$Metric_group='F1'
E2_performance_Recall$Metric_group='Recall'
E2_performance_Precision$Metric_group='Precision'
E2_performance_all = rbind(E2_performance_Recall,E2_performance_Precision)
E2_performance_all$Group = factor(E2_performance_all$Group,levels=c("No-filter","Ref-filter"))
E2_performance_all$Tool = factor(E2_performance_all$Tool,levels=c("svim","PBSV","cuteSV","Sniffles","PSDP"))
E2_performance_all$Metric_group = factor(E2_performance_all$Metric_group,levels=c("Precision","Recall"))

E2_performance_all$Group2 = paste(E2_performance_all$Metric_group,E2_performance_all$Group,sep='-')

pdf("Fig1a.pdf",width=6,height=5)
ggplot(E2_performance_all,aes(x=Tool,y=Metric)) +
    geom_point(aes(color=Metric_group,shape=Metric_group),size=3) +
    geom_line(aes(group=Group2,color=Metric_group,linetype=Group)) +
    theme_classic() +
        theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=60,size = 14, color = 'black',hjust=1),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.title = element_text(size=16,color="black"),
        legend.text = element_text(size=14,color="black"))+
        scale_color_manual(values = c(color[3],color[5],color[6],color[4])) +
    labs(x="Methods",y="Metric score",color=" ",shape=" ",linetype=" ") +
    guides(shape=FALSE)
dev.off()

# (2) Benchmarking on simulated virome-enriched datasets with "Uneven" condition
Precision=c(90.06,83.41,81.19,75.18,88.74)
Recall=c(30.67,28.60,17.47,39.45,27.56)
F1=c(45.76,28.76,51.75,42.06,44.73)
PreFilter_E1=data.frame(Precision=Precision,Recall=Recall,F1=F1,Tool=c("PSDP","PBSV","cuteSV","svim","Sniffles"),Group='Ref-filter')
Precision=c(91.90,84.60,84.09,78.71,90.37)
Recall=c(13.23,7.02,7.53,19.64,13.23)
F1=c(23.13,12.98,13.82,31.44,23.08)
NoFilter_E1=data.frame(Precision=Precision,Recall=Recall,F1=F1,Tool=c("PSDP","PBSV","cuteSV","svim","Sniffles"),Group='No-filter')

E1_performance = rbind(PreFilter_E1,NoFilter_E1)
E1_performance$Group = factor(E1_performance$Group,levels=c("No-filter","Ref-filter"))
E1_performance$Tool = factor(E1_performance$Tool,levels=c("svim","PBSV","cuteSV","Sniffles","PSDP"))


E1_performance_F1 = E1_performance[,c("Group","Tool","F1")]
colnames(E1_performance_F1) = c("Group","Tool","Metric")
E1_performance_Recall = E1_performance[,c("Group","Tool","Recall")]
colnames(E1_performance_Recall) = c("Group","Tool","Metric")
E1_performance_Precision = E1_performance[,c("Group","Tool","Precision")]
colnames(E1_performance_Precision) = c("Group","Tool","Metric")
E1_performance_F1$Metric_group='F1'
E1_performance_Recall$Metric_group='Recall'
E1_performance_Precision$Metric_group='Precision'
E1_performance_all = rbind(E1_performance_Recall,E1_performance_Precision)
E1_performance_all$Group = factor(E1_performance_all$Group,levels=c("No-filter","Ref-filter"))
E1_performance_all$Tool = factor(E1_performance_all$Tool,levels=c("svim","PBSV","cuteSV","Sniffles","PSDP"))
E1_performance_all$Metric_group = factor(E1_performance_all$Metric_group,levels=c("Precision","Recall"))

E1_performance_all$Group2 = paste(E1_performance_all$Metric_group,E1_performance_all$Group,sep='-')

pdf("FigS1a_uneven_simulated_benchmarking.pdf",width=6,height=5)
ggplot(E1_performance_all,aes(x=Tool,y=Metric)) +
    geom_point(aes(color=Metric_group,shape=Metric_group),size=3) +
    geom_line(aes(group=Group2,color=Metric_group,linetype=Group)) +
    theme_classic() +
        theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=60,size = 14, color = 'black',hjust=1),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.title = element_text(size=16,color="black"),
        legend.text = element_text(size=14,color="black"))+
        scale_color_manual(values = c(color[3],color[5],color[6],color[4])) +
    labs(x="Methods",y="Metric score",color=" ",shape=" ",linetype=" ") +
    guides(shape=FALSE)
dev.off()


######################################################=
# Pre-filtering VS. non-filtering on the real dataset
######################################################


data = read.csv("SV_stat.csv")
data = data[data$Group != 'Pre-filter(PB+NGS)',]
data$Group = factor(data$Group,levels=c("No-filter",'Pre-filter(PB)'),labels=c("No-filter",'Pre-filter'))

pdf("FigS1b_FilterVSnoFilter.pdf",width=4,height=5)
ggplot(data,aes(x=Group,y=Total_num)) +
    theme_bw() +
    geom_boxplot(aes(fill = Group), outlier.size=0, width = 0.6,alpha=0.6) +
    scale_fill_manual(values = c("#7472AE","#234D33","#C04D87")) +
    geom_point(size = 1,color='black') +
    geom_point(size = 1,shape=21) +
    geom_line(aes(group = sample), color = 'gray', lwd = 0.5) +
    theme(panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text.x = element_text(angle=60,size = 15, color = 'black',hjust=1),
        axis.text.y = element_text(size = 15, color = 'black'),
        axis.title = element_text(size = 15, color = 'black')) +
    labs(x=" ",y="SV number") +
    stat_compare_means(method = "t.test",paired = TRUE, comparisons=list(c("No-filter", "Pre-filter"))) +
    guides(fill=FALSE)
dev.off()

#######################################################
# Long-read SV calling VS. Short-read based SV calling
#######################################################

# Shell
# python3 scripts/PBvsNGS_stat.py --NGS_VCF_list NGS_VCF_path.txt --PacBio_VCF_list PB_VCF_path.txt --outputfile Results/NGS_vs_PB_sv_inf.csv

data = read.csv("Results/NGS_vs_PB_sv_inf.csv")

library(forcats)
library(ggpubr)

pdf("FigS1b_FilterVSnoFilter.pdf",width=3,height=5)
ggplot(data,aes(x=Group,y=Total)) +
    theme_bw() +
    geom_boxplot(aes(fill = Group), outlier.size=0, width = 0.6) +
    geom_point(size = 1,color='black') +
    geom_point(size = 1,shape=21) +
    geom_line(aes(group = sample), color = 'gray', lwd = 0.5) +
    theme(panel.grid = element_blank(),
    plot.title = element_text(size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 15, hjust = 0.5),
    axis.text.x = element_text(size = 15, color = 'black'),
    axis.text.y = element_text(size = 15, color = 'black'),
    panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
    axis.title = element_text(size = 15, color = 'black')) +
    labs(x=" ",y="SV number per sample") +
    scale_x_discrete(limits=c("Illumina","PacBio"),labels=c("NGS","PBS")) +
    stat_compare_means(method = "t.test",paired = TRUE, comparisons=list(c("PacBio", "Illumina"))) +
    scale_fill_manual(values = c(color[5],color[6])) +
    guides(fill=FALSE)
dev.off()


#######################################################
# Bar plots of all non-redundant viral SVs
#######################################################

Barplot_data = data.frame(Count = c(6690, 7014, 350, 284), SVTYPE=c("INS","DEL","DUP","INV"))
Barplot_data$SVTYPE = factor(Barplot_data$SVTYPE,levels=c("INS","DEL","DUP","INV"))

pdf("Fig1b.pdf",width=3,height=5)
ggplot(Barplot_data,aes(x=SVTYPE,y=Count)) +
    geom_bar(aes(fill=SVTYPE),stat='identity',color='black',width=0.8) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_rect(fill=NA,color="black",size=1,linetype="solid"),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title = element_text(size=16,color="black"),
    legend.text = element_text(size=14,color="black"),
    axis.ticks = element_blank(),
    strip.text = element_text(size=14,color="black")) +
    labs(x="SV types",y="SV number",fill=" ") +
    scale_fill_manual(values=sv_colors1) +
    guides(fill=FALSE)
dev.off()


#######################################################
# Length distribution of 14,338 non-redundant viral SVs
#######################################################

#python3 Scripts/VCFstat.py \
#    --VCF_file Results/Sample_SV_common_0.8_suppl.vcf \
#    --VCF_list Results/VCF_list.txt \
#    --outputfile Results/SV_bed_inf.tsv


ViralSV_inf = read.csv("Results/SV_bed_inf.tsv",sep="\t",row.names=1)
ViralSV_inf$SVTYPE = factor(ViralSV_inf$SVTYPE,levels=c("INS","DEL","DUP","INV"))

pdf("/Users/laisenying/Desktop/Fig1c.pdf",width=5,height=5)
ggplot(ViralSV_inf,aes(x=log2(abs(SVLEN)))) +
    facet_wrap(~SVTYPE,ncol=1)+
    geom_density(aes(fill=SVTYPE)) +
    theme_classic() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_text(color='black',size=16),
        axis.title.x = element_text(color='black',size=16),
        strip.text = element_text(size=16)
    ) +
    scale_fill_manual(values=sv_colors1) +
    labs(x="log2(SV length)") +
    guides(fill=FALSE)
dev.off()


#######################################################
# Correlation between SV and SNV density
#######################################################

# SNP all calcualted using inStrain

SNP_inf= read.csv("Results/all_genome_SNP_inf.tsv",sep="\t")
SNP_inf$SV_density = SNP_inf$SV_number/(SNP_inf$length/1000000)
cor.test(SNP_inf$SV_density,(SNP_inf$SNS_count+SNP_inf$SNV_count)/(SNP_inf$length*SNP_inf$breadth_minCov/1000000),na.rm=TRUE) # Cor = 0.61, pvalue

pdf("/Users/laisenying/Desktop/Part1_Fig1d.pdf",width=5,height=5)
ggplot(SNP_inf[SNP_inf$sample %in% unique(SNP_inf$sample)[1:10],],
    aes(x=SV_density,y=(SNS_count+SNV_count)/(length/1000000))) +
    geom_point(color=color[5]) +
    geom_smooth(method='lm',color=color[4],fill=color[4]) +
    theme_bw() +
    theme(
        panel.grid=element_blank(),
        axis.text=element_text(color='black',size=14),
        axis.title=element_text(color='black',size=16),
        panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
        plot.title=element_text(hjust=0.5)
    ) +
    labs(x="SV number per 1Mb",y="SNV number per 1Mb") +
    annotate("text",x=500,y=60000,label="Pearson's Cor=0.61, p=2.2e-16",size=5,color='black')
dev.off()


# Gene in SV region has significant higher nucleotide deversity

SNP_inf = read.csv("Results/all_gene_SNP_inf.tsv",sep="\t")
data_filter = SNP_inf[!is.na(SNP_inf$nucl_diversity),]
data_filter = data_filter[data_filter$breadth_minCov>0.4,]
data_filter = data_filter[data_filter$coverage>5,]
data_filter$SV.associated = factor(data_filter$SV.associated,
    levels=c("no_SV","noGE_SV","GE_SV"),labels=c("Gene in conserved region","Gene in SVs","Gene in SVs"))


pdf("Fig1f_1.pdf",width=5,height=5)
ggplot(data_filter,aes(x=log(coverage),y=nucl_diversity)) +
    geom_smooth(aes(group=SV.associated,color=SV.associated,fill=SV.associated)) +
    #geom_point(aes(color=SV.associated)) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=16,color="black"),
        panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
        legend.position = c(0.65,0.85),
        legend.text = element_text(size=14,color="black")
    ) +
    scale_color_manual(values=c(sv_colors1[4],sv_colors1[2])) +
    labs(x="log2(Coverage)",y="Nucleotide diversity",color=" ") +
    guides(fill=FALSE)
dev.off()

pdf("Fig1f_2.pdf",width=1,height=2)
ggplot(data_filter[sample(1:nrow(data_filter),10000),],aes(x=SV.associated,y=nucl_diversity)) +
    geom_boxplot(aes(fill=SV.associated,color=SV.associated),width=0.6,outlier.size=0) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
        legend.position = c(0.65,0.85),
        axis.ticks = element_blank(),
        legend.text = element_text(size=14,color="black")
    ) +
    scale_color_manual(values=c(sv_colors1[4],sv_colors1[2])) +
    scale_fill_manual(values=c(sv_colors1[4],sv_colors1[2])) +
    guides(fill=FALSE) +
    guides(color=FALSE) +
    geom_signif(comparisons=list(c("Gene in conserved region","Gene in SVs")))
dev.off()


#######################################################
# SV/SNV density across variable viral families
#######################################################

ViralSV_density_inf = read.csv("Results/ViralSV_density_inf.tsv",sep="\t",row.names=1)
ViralSV_density_inf$vOTU[ViralSV_density_inf$vOTU==''] = 'Others'
length(table(ViralSV_density_inf$vOTU))
ViralSV_density_inf$Family[ViralSV_density_inf$Family==''] = 'Others'
ViralSV_density_inf$Family[ViralSV_density_inf$Family=='Podoviridae '] = 'Podoviridae'
ViralSV_density_inf$Family[ViralSV_density_inf$Family=='Myoviridae '] = 'Myoviridae'
ViralSV_density_inf$Family[ViralSV_density_inf$Family=='crAss-like phage'] = 'crAssPhage'
ViralSV_density_inf$Family[ViralSV_density_inf$Family=='uc_o_Caudovirales'] = 'uc_Caudovirales'
ViralSV_density_inf$Family[ViralSV_density_inf$vOTU=='vOTU_41'] = 'uc_Caudovirales'
ViralSV_density_inf$Family[ViralSV_density_inf$vOTU=='vOTU_344'] = 'Gubaphage'
#ViralSV_density_inf$Family[ViralSV_density_inf$vOTU=='vOTU_344'] = 'vOTU_15'
ViralSV_density_inf$Family[ViralSV_density_inf$vOTU %in% names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='Myoviridae','vOTU']))[names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='Myoviridae','vOTU']))!='Others']] = 'Myoviridae'
ViralSV_density_inf$Family[ViralSV_density_inf$vOTU %in% names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='Podoviridae','vOTU']))[names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='Podoviridae','vOTU']))!='Others']] = 'Myoviridae'
ViralSV_density_inf$Family[ViralSV_density_inf$vOTU %in% names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='Podoviridae','vOTU']))[names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='Podoviridae','vOTU']))!='Others']] = 'Podoviridae'
ViralSV_density_inf$Family[ViralSV_density_inf$vOTU %in% names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='Myoviridae','vOTU']))[names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='Myoviridae','vOTU']))!='Others']] = 'Myoviridae'
ViralSV_density_inf$Family[ViralSV_density_inf$vOTU %in% names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='Siphoviridae','vOTU']))[names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='Siphoviridae','vOTU']))!='Others']] = 'Siphoviridae'
ViralSV_density_inf$Family[ViralSV_density_inf$vOTU %in% names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='uc_Caudovirales','vOTU']))[names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='uc_Caudovirales','vOTU']))!='Others']] = 'uc_Caudovirales'
ViralSV_density_inf$Family[ViralSV_density_inf$vOTU %in% names(table(ViralSV_density_inf[ViralSV_density_inf$Family=='crAssPhage','vOTU']))] = 'crAssPhage'

ViralSNP_density_inf = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/all_genome_SNP_inf.tsv",sep="\t",row.names=1)
ViralSNP_density_inf$SNV_density = ViralSNP_density_inf$SNV_count/(ViralSNP_density_inf$length/1000000)
SNP_inf = ViralSNP_density_inf[,c("Viral_contig","sample","nucl_diversity","SNV_density")]
merged_data = merge(ViralSV_density_inf,SNP_inf,by=c("Viral_contig","sample"))

sv = aggregate(SV_density~Family,ViralSV_density_inf,mean)
order_vOTUs = sv$Family[order(-sv$SV_density)]


color1 <- c(brewer.pal(9,'Set1'))
color2 <- c(brewer.pal(8,'Set2'))
color3 <- c(brewer.pal(12,'Set3'))

# boxplot showing the SV density across viral clades
pSNV=ggplot(merged_data,aes(x=Family,y=SNV_density)) +
    #geom_point(aes(color=Family))
    geom_boxplot(aes(fill=Family),outlier.size=0,width=0.7) +
    #geom_jitter(aes(color=Family),width=0.2,size=1,alpha=1) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14,color="black"),
        axis.title = element_text(size=16,color="black"),
        plot.title = element_text(hjust=0.5,size=18),
        panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=14,color="black"),
        axis.ticks = element_blank(),
        legend.position=c(0.78,0.76)
    ) +
    labs(x="Family",y="SNV number per 1Mb") +
    scale_fill_manual(values=c(color3[3:6],color2,color1)) +
    scale_color_manual(values=c(color3[3:6],color2,color1)) +
    #scale_y_continuous(breaks=seq(0, 0.002, 0.004))  +
    scale_x_discrete(limits=order_vOTUs[order_vOTUs != 'Microviridae']) +
    guides(fill=FALSE) +
    guides(color=FALSE) +
    coord_flip() +
    geom_signif(comparisons=list(c("Myoviridae","uc_Caudovirales")))

pSV=ggplot(ViralSV_density_inf,aes(x=Family,y=SV_density)) +
    #geom_point(aes(color=Family))
    geom_boxplot(aes(fill=Family),outlier.size=0,width=0.7) +
    #geom_jitter(aes(color=Family),width=0.2,size=1,alpha=1) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14,color="black"),
        axis.title = element_text(size=16,color="black"),
        plot.title = element_text(hjust=0.5,size=18),
        panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=14,color="black"),
        axis.ticks = element_blank(),
        legend.position=c(0.78,0.76)
    ) +
    labs(x="Family",y="SV number per 1Mb") +
    scale_fill_manual(values=c(color3[3:6],color2,color1)) +
    scale_color_manual(values=c(color3[3:6],color2,color1)) +
    #scale_y_continuous(breaks=seq(0, 0.002, 0.004))  +
    scale_x_discrete(limits=order_vOTUs[order_vOTUs != 'Microviridae']) +
    guides(fill=FALSE) +
    guides(color=FALSE) +
    coord_flip()

pdf("/Users/laisenying/Desktop/Fig1g.pdf",width=10,height=5)
plot_grid(pSV,pSNV,nrow=1,rel_widths=c(0.6,0.6))
dev.off()


### Phage genomes selected for each sample
import pandas as pd
import numpy as np
import os
from Bio import SeqIO

sample_list = pd.read_csv("sample.txt",sep="\t",header=None)

all_genome_data = pd.DataFrame({})
for sample in list(sample_list[0]):
    if os.path.exists("PBSV/"+sample+"/Mash_filter/screen.tab"):
        sample_data = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/PBSV/"+sample+"/Mash_filter/screen.tab",sep="\t",header=None)
        sample_data.columns = ['Identity','Shared-hashes','median-multiplicity','p-value','Genome_fasta','Genome']
        sample_data_filter = sample_data.loc[sample_data['Identity']>0.95,]
        sample_data_filter['sample'] = sample
        sample_data_filter = sample_data_filter.loc[:,['Genome','sample']]
        all_genome_data  = pd.concat([all_genome_data,sample_data_filter])
all_genome_data.to_csv("Results/Genomes_selected_for_sample.tsv",sep="\t") # 9,183



