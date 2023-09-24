###############################################################
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(forcats)
library(ggpubr)
library(ggsignif)
library(reshape2)
library(ggalluvial)

library(vegan)
library(ade4)
library(cluster)
library(vegan)
library(ggdendro)
library(sparcl)
library(factoextra)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(ggsci)
library(vegan)
library(cowplot)


color <- c(brewer.pal(12,'Set3'))
c(color[3],color[5],color[6],color[4])

sv_colors1=c("#339933", "#336699", "#CCCC33", "#CC6633")
sv_colors2=c("#CC6666", "#999933", "#339999", "#996699")
sv_colors3 = c("#FFCCCC", "#99CCFF", "#FFCC99", "#CCCCFF", "#99CCCC", "#FFCC99", "#CCFFCC")
sv_colors4=c("#9999CC","#CC9999","#99CC99")
sv_colors5 = c("#336699","#CCCC33")
sv_colors6 = c("#CC6666","#6699CC")


##########################################################################
####### The enrichment of SV transmission between phage and hosts ########
##########################################################################

# ! python3
Host_data = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/Host_identify/Bac2Viral_host_CRSPR.tsv",sep="\t",index_col=0)
BV_data = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_coevolution/Interactive_range/BV_all_interactive_data.tsv",sep="\t",index_col=0)
BV_data_filter = BV_data.loc[BV_data['Bacterial_Source'] == 'CH-Binning',]
BV_data_filter = BV_data_filter.loc[BV_data_filter['Viral_Source']=='HGV',]
BV_data_filter['interact_id'] = BV_data_filter['Viral_contig'] + "~" + BV_data_filter['Genus']
from collections import defaultdict
BV_inf = defaultdict(list)
for interact in set(BV_data_filter['interact_id']):
    subdata = pd.DataFrame(BV_data_filter.loc[BV_data_filter['interact_id']==interact,])
    BV_inf['Viral_contig'].append(list(subdata['Viral_contig'])[0])
    BV_inf['Taxonomy'].append(list(subdata['Taxonomy'])[0])
    BV_inf['Phylum'].append(list(subdata['Phylum'])[0])
    BV_inf['Class'].append(list(subdata['Class'])[0])
    BV_inf['Order'].append(list(subdata['Order'])[0])
    BV_inf['Family'].append(list(subdata['Family'])[0])
    BV_inf['Genus'].append(list(subdata['Genus'])[0])
    BV_inf['Species'].append(list(subdata['Species'])[0])
    BV_inf['B-to-V count'].append(len(set(subdata.loc[subdata['Direction']=='bacteria-to-viralSV','NR_SVID'])))
    BV_inf['V-to-B count'].append(len(set(subdata.loc[subdata['Direction']=='viral-to-bacSV','NR_SVID'])))

BV_inf_data = pd.DataFrame(BV_inf)
BV_inf_data['interact_id'] = BV_inf_data['Viral_contig'] + "~" + BV_inf_data['Genus']
BV_inf_data_filter = BV_inf_data.loc[BV_inf_data['B-to-V count']>0,]
Host_data['interact_id'] = Host_data['Viral_contig'] + "~" + Host_data["Bac_Genus"]


Host_data_filter = Host_data.loc[[x in set(BV_inf_data_filter['Viral_contig']) for x in Host_data['Viral_contig']],]


set(Host_data_filter['interact_id']).intersection(BV_inf_data_filter['interact_id']) # 31
len(set(Host_data_filter['interact_id'])) # 43

HGV_bac_inf = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/Bac_SV/Bacterial_genomes/bacterial_genome_inf.tsv",sep='\t',index_col=0)
HGV_bac_inf['Genus'] = [x.split(";")[5] for x in HGV_bac_inf['Taxonomy']]
HGV_genus = set(HGV_bac_inf['Genus'])


BV_inf_data_filter['Host_detected'] = [x in set(Host_data_filter['interact_id']) for x in BV_inf_data_filter['interact_id']]

all_possible_interact_ids = []
for x in set(BV_inf_data_filter['Viral_contig']):
    for y in HGV_genus:
        all_possible_interact_ids.append(x+"~"+y)


len(set(all_possible_interact_ids)) # 68460

# non Host-phage pairs: len(set(all_possible_interact_ids).difference(Host_data['interact_id'])) # 68417
# SV transmission between non Host-phage pairs: (set(all_possible_interact_ids).difference(Host_data['interact_id'])).intersection(BV_inf_data_filter['interact_id']) # 819


# plot
plot_data = data.frame(
    Group = c("PH pairs","non-PH pairs"),
    GE=c(31,819),
    nonGE = c(12,68417-819)
    )
colnames(plot_data) = c("Group","GE","nonGE")
plot_data[,2:3] = plot_data[,2:3]/rowSums(plot_data[,2:3])
df1 <- melt(plot_data,id.vars = 'Group',measure.vars = c("GE","nonGE"))
names(df1)[1:2] <- c("X","group")
df1$group = factor(df1$group,levels=c("nonGE","GE"))

color <- c(brewer.pal(12,'Set3'))
    category.colors <- c(
    "nonGE" = color[5],
    "GE" = color[4]
)

pdf("/Users/laisenying/Desktop/Fig4a.pdf",width=4.7,height=5.5)
ggplot(df1, aes( x = X,y=100 * value,fill = group,
    stratum = group, alluvium = group))+
    geom_stratum(width = 0.5, color='black')+
    geom_alluvium(alpha = 0.5,
    width = 0.5,
    curve_type = "linear") +
    theme_bw() +
    theme(panel.grid = element_blank(),
    panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
    axis.text.x = element_text(size = 14, color = 'black',angle=60,hjust=1),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title = element_text(size=16,color="black"),
    legend.text = element_text(size=14,color="black"),
    legend.title = element_text(size=14,color="black"),
    axis.ticks = element_blank(),
    strip.text = element_text(size=14,color="black"),
    plot.title = element_text(hjust=0.5,size=16)) +
    scale_fill_manual(values=as.vector(category.colors),limits=names(category.colors)) +
    labs(x=" ",y="Percentage",fill="SV transmission")
dev.off()

##########################################################################
########## Specificity of SV exchange between phage-host pairs ###########
##########################################################################


color <- c(brewer.pal(12,'Set3'))

color.categories = c()
i=1
for(x in c("Phylum","Class","Order","Family","Genus","Species","Kingdom")){
    color.categories[x] = color[i]
    i= i+ 1
}

color.categories = c(
    "Species"="#389E39",
    "Genus" = "#CCDD8D",
    "Family"="#2577AC",
    "Order" = "#B8D2E2",
    "Class" = "#CD2B20",
    "Phylum" = "#F2C173",
    "Kingdom"="#ECA1A5"
)

data_plot = data.frame(
    Count = c(52,55,21,9,8,1,3),
    Range = c("Species","Genus","Family","Order","Class","Phylum","Kingdom"),
    Direction = c("B-to-V","B-to-V","B-to-V","B-to-V","B-to-V","B-to-V","B-to-V")
)
data_plot$Range = factor(data_plot$Range,levels=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))


pdf("/Users/laisenying/Desktop/Fig4b.pdf",width=8,height=5.5)
ggplot(data_plot,aes(x=Direction,y=Count)) +
    theme_bw() +
    geom_bar(aes(fill=Range),stat='identity',position='fill',color='black',width=0.8) +
    theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size=14,color='black',angle=60,hjust=1),
    axis.text.y = element_text(size=14,color='black'),
    axis.title = element_text(size=16,color="black"),
    plot.title = element_text(hjust=0.5,size=18),
    panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
    legend.text = element_text(size=14,color="black"),
    legend.title = element_text(size=14,color="black"),
    axis.ticks = element_blank()
    ) +
    labs(x=" ",y="Proportion",fill="Interactive range") +
    scale_fill_manual(values=as.vector(color.categories),limits=names(color.categories)) +
    coord_polar(theta = 'y',start = pi )
dev.off()

##########################################################################
################ The influence of Host range on SV density ###############
##########################################################################

import pandas as pd
import numpy as np
from collections import defaultdict
from Bio import SeqIO

BV_data = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_coevolution/Interactive_range/BV_all_interactive_data.tsv",sep="\t",index_col=0)
BV_data_filter = BV_data.loc[BV_data['Direction']=='bacteria-to-viralSV',]
BV_data_filter = BV_data_filter.loc[BV_data_filter['Viral_Source']=='HGV',]

Viral_SV_inf = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/PacBio_SV/CAST_merge/Population/Result/SV_bed_inf.tsv",sep="\t",index_col=0)

Input = SeqIO.parse("/share/home1/Laisenying/Data-analysis/projects/PhageSV/ViralGenome/HGV.hq.genome.fa","fasta")
contig_len={}
for record in Input:
    contig_len[record.id] = len(record.seq)
        

Host_range_inf = defaultdict(list)
for contig in set(BV_data_filter['Viral_contig']):
    subdata = pd.DataFrame(BV_data_filter.loc[BV_data_filter['Viral_contig']==contig,])
    if subdata.shape[1] == 1:
        subdata = subdata.T
    if len(set(subdata['SV_name'])) < 2:
        continue
    subdata = subdata.drop_duplicates(['SV_name','Taxonomy'])
    subdata = subdata.iloc[np.argsort(-subdata['identity']),]
    subdata = subdata.drop_duplicates(['SV_name'])
    Host_range_inf['Viral_contig'].append(contig)
    Host_range_inf['Contig_length'].append(contig_len[contig])
    Host_range_inf['GESV_num'].append(len(set(subdata['SV_name'])))
    Host_range_inf['all_SV_count'].append(len(set(Viral_SV_inf.loc[Viral_SV_inf['contig']==contig,'SVID'])))
    Host_range_inf['Species_num'].append(len(set(subdata['Species'])))
    Host_range_inf['Species_purity'].append(pd.value_counts(subdata['Species'])[0]/len(set(subdata['SV_name'])))
    Host_range_inf['Genus_num'].append(len(set(subdata['Genus'])))
    Host_range_inf['Genus_purity'].append(pd.value_counts(subdata['Genus'])[0]/len(set(subdata['SV_name'])))
    Host_range_inf['Family_num'].append(len(set(subdata['Family'])))
    Host_range_inf['Family_purity'].append(pd.value_counts(subdata['Family'])[0]/len(set(subdata['SV_name'])))
    Host_range_inf['Order_num'].append(len(set(subdata['Order'])))
    Host_range_inf['Order_purity'].append(pd.value_counts(subdata['Order'])[0]/len(set(subdata['SV_name'])))
    Host_range_inf['Class_num'].append(len(set(subdata['Class'])))
    Host_range_inf['Class_purity'].append(pd.value_counts(subdata['Class'])[0]/len(set(subdata['SV_name'])))
    Host_range_inf['Phylum_num'].append(len(set(subdata['Phylum'])))
    Host_range_inf['Phylum_purity'].append(pd.value_counts(subdata['Phylum'])[0]/len(set(subdata['SV_name'])))
        

Host_range_inf = pd.DataFrame(Host_range_inf)
Host_range_inf['Host_range'] = 'Kingdom'
Host_range_inf.loc[Host_range_inf['Phylum_purity']==1,'Host_range'] = 'Phylum'
Host_range_inf.loc[Host_range_inf['Class_purity']==1,'Host_range'] = 'Class'
Host_range_inf.loc[Host_range_inf['Order_purity']==1,'Host_range'] = 'Order'
Host_range_inf.loc[Host_range_inf['Family_purity']==1,'Host_range'] = 'Family'
Host_range_inf.loc[Host_range_inf['Genus_purity']==1,'Host_range'] = 'Genus'
Host_range_inf.loc[Host_range_inf['Species_purity']==1,'Host_range'] = 'Species'
Host_range_inf.to_csv("/home1/Laisenying/Projects/PhageSV/Results/HGV_host_range_inf.tsv",sep="\t")


scp Laisenying@10.190.248.213:/home1/Laisenying/Projects/PhageSV/Results/HGV_host_range_inf.tsv '/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/HGV_host_range_inf.tsv'


Host_range_data = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/HGV_host_range_inf.tsv",sep="\t",row.names=1)
Host_range_data['SV_density'] = Host_range_data['all_SV_count']/(Host_range_data['Contig_length']/1000000)

color.categories = c(
    "Species"="#389E39",
    "Genus" = "#CCDD8D",
    "Family"="#2577AC",
    "Order" = "#B8D2E2",
    "Class" = "#CD2B20",
    "Phylum" = "#F2C173",
    "Kingdom"="#ECA1A5"
)

pdf("/Users/laisenying/Desktop/Fig4c.pdf",width=6,height=5.1)
ggplot(Host_range_data[Host_range_data['GESV_num']>=4,],aes(x=Host_range,y=SV_density)) +
    geom_violin(aes(fill=Host_range)) +
    geom_boxplot(width=0.15,aes(fill=Host_range),width=0.8,outlier.size=0) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size=14,color='black',angle=60,hjust=1),
        axis.text.y = element_text(size=14,color='black'),
        axis.title = element_text(size=16,color="black"),
        plot.title = element_text(hjust=0.5,size=18),
        panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=14,color="black"),
        axis.ticks = element_blank()
    ) +
    scale_fill_manual(values=as.vector(color.categories),limits=names(color.categories)) +
    labs(x="Host range",y="SV number per 1Mb") +
    scale_x_discrete(limits=c("Kingdom","Class","Order","Family","Genus","Species")) +
    guides(fill=FALSE) +
    geom_signif(comparisons=list(c("Kingdom","Family")),y_position=1250+150) +
    geom_signif(comparisons=list(c("Kingdom","Genus")),y_position=1150+150) +
    geom_signif(comparisons=list(c("Kingdom","Species")),y_position=1050+150) +
    geom_signif(comparisons=list(c("Class","Species")),y_position=950+150) +
    geom_signif(comparisons=list(c("Genus","Species")),y_position=850+150)
dev.off()


##########################################################################
################# Construction of BV interactive network #################
##########################################################################

import pandas as pd
import numpy as np

BV_interactive_data = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_coevolution/Interactive_range/BV_all_interactive_data.tsv",sep="\t",index_col=0)
BV_interactive_data_CHGV = BV_interactive_data.loc[BV_interactive_data['Viral_Source']=='HGV',]
BV_interactive_data_CHGV = BV_interactive_data_CHGV.loc[BV_interactive_data_CHGV['Bacterial_Source']=='CH-Binning',]
BV_interactive_data_CHGV = BV_interactive_data_CHGV.loc[BV_interactive_data_CHGV['Direction']=='bacteria-to-viralSV',]
BV_interactive_data_CHGV_filter = BV_interactive_data_CHGV.drop_duplicates(['NR_SVID','Rep_genome'])
BV_interactive_data_CHGV_filter = BV_interactive_data_CHGV_filter.iloc[np.argsort(-BV_interactive_data_CHGV_filter['identity']),]
BV_interactive_data_CHGV_filter = BV_interactive_data_CHGV_filter.drop_duplicates(['NR_SVID'])

BV_interactive_data_CHGV_filter['Interaction_id'] = BV_interactive_data_CHGV_filter['Genus'] +"~"+ BV_interactive_data_CHGV_filter['Viral_contig']
pd.value_counts(BV_interactive_data_CHGV_filter['Interaction_id'])
BV_interactive_data_CHGV_filter.index = BV_interactive_data_CHGV_filter['Interaction_id']
BV_interactive_data_CHGV_filter['Interact_num'] = pd.value_counts(BV_interactive_data_CHGV_filter['Interaction_id'])[BV_interactive_data_CHGV_filter['Interaction_id']]
BV_interactive_data_CHGV_filter['Count'] = 1
grouped = BV_interactive_data_CHGV_filter.groupby(['Viral_contig'])
BV_interactive_data_CHGV_filter.index = BV_interactive_data_CHGV_filter['Viral_contig']
BV_interactive_data_CHGV_filter['Viral_hit_num'] = (grouped['Count'].sum())[BV_interactive_data_CHGV_filter['Viral_contig']]
BV_interactive_data_CHGV_filter['Interactive_weight'] = BV_interactive_data_CHGV_filter['Interact_num']/BV_interactive_data_CHGV_filter['Viral_hit_num']
BV_interactive_data_CHGV_filter2 = BV_interactive_data_CHGV_filter.drop_duplicates(['Interaction_id'])
BV_interactive_data_CHGV_filter2.index = range(BV_interactive_data_CHGV_filter2.shape[0])

selected_genus = pd.value_counts(BV_interactive_data_CHGV_filter2['Genus'])[pd.value_counts(BV_interactive_data_CHGV_filter2['Genus'])>=10].index
BV_interactive_data_CHGV_filter3 = BV_interactive_data_CHGV_filter2.iloc[[x in selected_genus for x in BV_interactive_data_CHGV_filter2['Genus']],]


BV_interactive_data_CHGV_filter3.to_csv("/home1/Laisenying/Projects/PhageSV/Results/BV_Network_plot.tsv",sep="\t")

scp Laisenying@10.190.248.213:/home1/Laisenying/Projects/PhageSV/Results/BV_Network_plot.tsv '/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/BV_Network_plot.tsv'

#########################################################################################
################# The influence of lifestyle on phage diversification ###################
#########################################################################################

# 1. Temperate phage contains higher microdiversity

scp Laisenying@10.190.248.213:/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/Lifestyle_analysis/HGV_Viral_SV_density.csv /Users/laisenying/Desktop/HGV_Viral_SV_density.csv

data_HGV = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/HGV_Viral_SV_density.csv")

data_HGV$lifestyle[data_HGV$lifestyle=='uncertain temperate'] = 'Temperate'
data_HGV$lifestyle[data_HGV$lifestyle=='temperate'] = 'Temperate'
data_HGV$lifestyle[data_HGV$lifestyle=='virulent'] = 'Virulent'
data_HGV$lifestyle[data_HGV$lifestyle=='uncertain virulent'] = 'Virulent'

color <- c(brewer.pal(12,'Set3'))
category.colors <- c(
    "Virulent" = color[2],
    "Temperate" = color[5]
)

p1 = ggplot(data_HGV,aes(x=lifestyle,y=SV_density)) +
    theme_bw() +
    #geom_violin(aes(fill=lifestyle)) +
    geom_boxplot(aes(fill=lifestyle),outlier.size=0,width=0.6) +
    theme(panel.grid = element_blank(),
    panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
    axis.text.x = element_text(size = 14, color = 'black',angle=60,hjust=1),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title = element_text(size=16,color="black"),
    legend.text = element_text(size=14,color="black"),
    legend.title = element_text(size=14,color="black"),
    axis.ticks = element_blank(),
    strip.text = element_text(size=14,color="black"),
    plot.title = element_text(hjust=0.5,size=16))+
    labs(x=" ",y='SV number per 1Mb',title='CHGV') +
    scale_fill_manual(values=as.vector(category.colors),limits=names(category.colors)) +
    geom_signif(comparisons=list(c("Temperate","Virulent"))) +
    guides(fill=FALSE) +
    scale_x_discrete(limits=c("Temperate","Virulent"))

IMGVR_SV_inf = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/IMGVR_SV_summary.tsv",sep="\t",row.names=1)
IMGVR_SV_inf$lifestyle[IMGVR_SV_inf$lifestyle=='temperate'] = 'Temperate'
IMGVR_SV_inf$lifestyle[IMGVR_SV_inf$lifestyle=='virulent'] = 'Virulent'
#IMGVR_SV_inf$lifestyle[IMGVR_SV_inf$lifestyle=='uncertain virulent'] = 'Virulent'
IMGVR_SV_inf = IMGVR_SV_inf[IMGVR_SV_inf$Query_environment=='human-gut',]
IMGVR_SV_inf = IMGVR_SV_inf[IMGVR_SV_inf$Ref_environment=='human-gut',]

p2 = ggplot(IMGVR_SV_inf,aes(x=lifestyle,y=SV_density)) +
    theme_bw() +
    #geom_violin(aes(fill=lifestyle)) +
    geom_boxplot(aes(fill=lifestyle),outlier.size=0,width=0.6) +
    theme(panel.grid = element_blank(),
    panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
    axis.text.x = element_text(size = 14, color = 'black',angle=60,hjust=1),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title = element_text(size=16,color="black"),
    legend.text = element_text(size=14,color="black"),
    legend.title = element_text(size=14,color="black"),
    axis.ticks = element_blank(),
    strip.text = element_text(size=14,color="black"),
    plot.title = element_text(hjust=0.5,size=16))+
    labs(x=" ",y='SV number per 1Mb',title='IMG/VR') +
    scale_fill_manual(values=as.vector(category.colors),limits=names(category.colors)) +
    geom_signif(comparisons=list(c("Temperate","Virulent"))) +
    guides(fill=FALSE) +
    scale_x_discrete(limits=c("Temperate","Virulent"))


pdf("/Users/laisenying/Desktop/Fig12_microdiversity.pdf",width=6.5,height=5.5)
plot_grid(p1,p2,nrow=1)
dev.off()


# 2. More frequency of transmission of GE-like SVs
plot_data = data.frame(
Group = c("Background","B-to-V"),
    Temperate=c(4552+2611,262+268),
    Virulent = c(1144+1358,86+68)
)
fisher.test(rbind(c(530,154),c(7163,2602))) # OR = 1.25, p = 0.017

colnames(plot_data) = c("Group","Temperate","Virulent")
plot_data[,2:3] = plot_data[,2:3]/rowSums(plot_data[,2:3])
df1 <- melt(plot_data,id.vars = 'Group',measure.vars = c("Temperate","Virulent"))
names(df1)[1:2] <- c("X","group")
df1$group = factor(df1$group,levels=c("Virulent","Temperate"))

color <- c(brewer.pal(12,'Set3'))
category.colors <- c(
    "Virulent" = color[2],
    "Temperate" = color[5]
)

p1 = ggplot(df1, aes( x = X,y=100 * value,fill = group,
    stratum = group, alluvium = group))+
    geom_stratum(width = 0.5, color='black')+
    geom_alluvium(alpha = 0.5,
    width = 0.5,
    curve_type = "linear") +
    theme_bw() +
    theme(panel.grid = element_blank(),
    panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
    axis.text.x = element_text(size = 14, color = 'black',angle=60,hjust=1),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title = element_text(size=16,color="black"),
    legend.text = element_text(size=14,color="black"),
    legend.title = element_text(size=14,color="black"),
    axis.ticks = element_blank(),
    strip.text = element_text(size=14,color="black"),
    plot.title = element_text(hjust=0.5,size=16)) +
    scale_fill_manual(values=as.vector(category.colors),limits=names(category.colors)) +
    scale_x_discrete(limits=c("Background","B-to-V")) +
    labs(x=" ",y="Percentage",fill="Lifestyle",title="CHGV")



plot_data = data.frame(
Group = c("Background","B-to-V"),
    Temperate=c(15695+16367,268+262),
    Virulent = c(20078+12444,68+86)
)

fisher.test(rbind(c(268+262,68+86),c(15695+16367,20078+12444))) # OR = 3.49, p < 2.2e-16

colnames(plot_data) = c("Group","Temperate","Virulent")
plot_data[,2:3] = plot_data[,2:3]/rowSums(plot_data[,2:3])
df2 <- melt(plot_data,id.vars = 'Group',measure.vars = c("Temperate","Virulent"))
names(df2)[1:2] <- c("X","group")
df2$group = factor(df1$group,levels=c("Virulent","Temperate"))



p2 = ggplot(df2, aes( x = X,y=100 * value,fill = group,
    stratum = group, alluvium = group))+
    geom_stratum(width = 0.5, color='black')+
    geom_alluvium(alpha = 0.5,
    width = 0.5,
    curve_type = "linear") +
    theme_bw() +
    theme(panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
        axis.text.x = element_text(size = 14, color = 'black',angle=60,hjust=1),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.title = element_text(size=16,color="black"),
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=14,color="black"),
        axis.ticks = element_blank(),
        strip.text = element_text(size=14,color="black"),
        plot.title = element_text(hjust=0.5,size=16)) +
    scale_fill_manual(values=as.vector(category.colors),limits=names(category.colors)) +
    scale_x_discrete(limits=c("Background","B-to-V")) +
    labs(x=" ",y="Percentage",fill="Lifestyle",title="IMG/VR")

pdf("/Users/laisenying/Desktop/Fig11.pdf",width=8.5,height=5.5)
plot_grid(p1,p2,nrow=1)
dev.off()


#########################################################################################
############################# The influence of environment ###############################
#########################################################################################

Gut_num=19844
Oral_num = 281
other_num=1055
IMGVRSV_Gut_num=2920
IMGVRSV_Oral_num = 4
IMGVR_other_num = 92

plot_data = data.frame(
    Group = c("Background","Viral recipients"),
    Gut = c(Gut_num,IMGVRSV_Gut_num),
    Oral = c(Oral_num,IMGVRSV_Oral_num),
    Others = c(other_num,IMGVR_other_num)
)
colnames(plot_data) = c("Group","Human-gut","Human-oral","Others")
plot_data[,2:4] = plot_data[,2:4]/rowSums(plot_data[,2:4])
df1 <- melt(plot_data,id.vars = 'Group',measure.vars = c("Human-gut","Human-oral","Others"))


color <- c(brewer.pal(12,'Set3'))
category.colors <- c(
    "Others" = color[4],
    "Human-oral" = color[5],
    "Human-gut" = color[6]
)
names(df1)[1:2] <- c("X","group")
df1$group = factor(df1$group,levels=c("Others","Human-oral","Human-gut"))

fisher.test(rbind(c(IMGVRSV_Gut_num,IMGVRSV_Oral_num+IMGVR_other_num),c(Gut_num,Oral_num+other_num))) # p = 2.891e-13, OR = 2.047743

p1=ggplot(df1, aes( x = X,y=100 * value,fill = group,
    stratum = group, alluvium = group))+
    geom_stratum(width = 0.5, color='black')+
    geom_alluvium(alpha = 0.5,
    width = 0.5,
    curve_type = "linear") +
    theme_bw() +
    theme(panel.grid = element_blank(),
    panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
    axis.text.x = element_text(size = 14, color = 'black',angle=60,hjust=1),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title = element_text(size=16,color="black"),
    legend.text = element_text(size=14,color="black"),
    legend.title = element_text(size=14,color="black"),
    axis.ticks = element_blank(),
    strip.text = element_text(size=14,color="black"),
    plot.title = element_text(hjust=0.5,size=16)) +
    scale_fill_manual(values=as.vector(category.colors),limits=names(category.colors)) +
    scale_x_discrete(limits=c("Background","Viral recipients")) +
    labs(x=" ",y="Percentage",fill="Environment",title=" ")

pdf("/Users/laisenying/Desktop/Fig6_env.pdf",width=5,height=5.7)
plot_grid(p1,nrow=1)
dev.off()


##########################################################################
####### Replication results of IMG/VR ########
##########################################################################

IMGVR_SV_inf = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/IMGVR_Viral_SV_sequence_inf.tsv",sep="\t",row.names=1)


# Summary of SV count
Barplot_data = data.frame(Count = c(8810, 11712, 598, 60), SVTYPE=c("INS","DEL","DUP","INV"))
Barplot_data$SVTYPE = factor(Barplot_data$SVTYPE,levels=c("INS","DEL","DUP","INV"))


pdf("/Users/laisenying/Desktop/FigS5a.pdf",width=4,height=5)
ggplot(Barplot_data,aes(x=SVTYPE,y=Count)) +
    geom_bar(aes(fill=SVTYPE),stat='identity',color='black',width=0.7) +
    geom_text(aes(label=Count,y=Count+400),size=5) +
    theme_bw() +
    theme(panel.grid = element_blank(),
    axis.text.x = element_text(size=14,color='black'),
    panel.border = element_rect(fill=NA,color="black",size=1,linetype="solid"),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title = element_text(size=16,color="black"),
    legend.text = element_text(size=14,color="black"),
    axis.ticks = element_blank(),
    strip.text = element_text(size=14,color="black")) +
    labs(x=" ",y="SV number",fill=" ") +
    scale_fill_manual(values=sv_colors1) +
    guides(fill=FALSE)
dev.off()


# SV length distribution
IMGVR_SV_inf = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/IMGVR_Viral_SV_sequence_inf.tsv",sep="\t",row.names=1)
IMGVR_SV_inf$SV_type = factor(IMGVR_SV_inf$SV_type,levels=c("INS","DEL","DUP","INV"))

pdf("/Users/laisenying/Desktop/FigS5b.pdf",width=5,height=5)
ggplot(IMGVR_SV_inf,aes(x=log2(abs(SV_length)))) +
    facet_wrap(~SV_type,ncol=1,strip.position='left')+
    geom_density(aes(fill=SV_type)) +
    theme_classic() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size=14,color='black'),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color='black',size=16),
        axis.title.x = element_text(color='black',size=16),
        strip.text = element_text(size=16)
    ) +
    scale_fill_manual(values=sv_colors1) +
    labs(x="log2(SV length)",y=" ") +
    guides(fill=FALSE)
dev.off()

# SV density across different phylum
IMGVR_SV_inf = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/IMGVR_SV_summary.tsv",sep="\t",row.names=1)

env_mean = aggregate(SV_density~Query_environment,IMGVR_SV_inf,mean)
order_envs = env_mean$Query_environment[order(env_mean$SV_density)]
order_envs = order_envs[order_envs!='']
order_envs = order_envs[order_envs!='Unknown']
color1 <- c(brewer.pal(9,'Set1'))
color2 <- c(brewer.pal(8,'Set2'))

pdf("/Users/laisenying/Desktop/FigS5c.pdf",width=7,height=5.5)
ggplot(IMGVR_SV_inf,aes(x=Query_environment,y=SV_density)) +
    geom_boxplot(aes(fill=Query_environment),width=0.6,outlier.size=0) +
    #geom_jitter(width=0.5,aes(color=Query_environment)) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
        axis.text.x = element_text(size = 14, color = 'black',angle=60,hjust=1),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.title = element_text(size=16,color="black"),
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=14,color="black"),
        axis.ticks = element_blank(),
        strip.text = element_text(size=14,color="black"),
        plot.title = element_text(hjust=0.5,size=16)) +
    labs(x="Environment of phages",y="SV number per 1 Mb") +
    scale_x_discrete(limits=order_envs) +
    scale_fill_manual(values=c(color1,color2)) +
    guides(fill=FALSE) +
    guides(color=FALSE) +
    geom_signif(comparisons=list(c("Freshwater","human-gut"))) +
    geom_signif(comparisons=list(c("Marine","human-gut")),y_position=220)
dev.off()


#
IMGVR_SV_inf = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/IMGVR_SV_summary.tsv",sep="\t",row.names=1)
IMGSV_inf = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/IMGVR_all_Sequence_information_high_quality.tsv",sep="\t",row.names=1)
rownames(IMGSV_inf) = IMGSV_inf$X...UViG
IMGVR_SV_inf$phage_Taxonomy = IMGSV_inf[IMGVR_SV_inf$Viral_contig,'Taxonomic.classification']
IMGVR_SV_inf_filter = IMGVR_SV_inf[IMGVR_SV_inf$phage_Taxonomy!='',]
IMGVR_SV_inf_filter$Family = unlist(lapply(strsplit(IMGVR_SV_inf_filter$phage_Taxonomy,";"),function(x){x[5]}))
IMGVR_SV_inf_filter$order = unlist(lapply(strsplit(IMGVR_SV_inf_filter$phage_Taxonomy,";"),function(x){x[4]}))
IMGVR_SV_inf_filter$class = unlist(lapply(strsplit(IMGVR_SV_inf_filter$phage_Taxonomy,";"),function(x){x[3]}))
IMGVR_SV_inf_filter$phylum = unlist(lapply(strsplit(IMGVR_SV_inf_filter$phage_Taxonomy,";"),function(x){x[2]}))


IMGVR_SV_inf_filter[(IMGVR_SV_inf_filter$Family=='') & (IMGVR_SV_inf_filter$phylum!=''),'Family'] = paste('uc_',IMGVR_SV_inf_filter[(IMGVR_SV_inf_filter$Family=='') & (IMGVR_SV_inf_filter$phylum!=''),'phylum'],sep='')

color <- c(brewer.pal(12,'Set3'))

f_mean = aggregate(SV_density~Family,IMGVR_SV_inf_filter,mean)
order_fs = f_mean$Family[order(f_mean$SV_density)]


pdf("/Users/laisenying/Desktop/FigS5d.pdf",width=7,height=5.8)
ggplot(IMGVR_SV_inf_filter,aes(x=Family,y=SV_density)) +
    geom_boxplot(aes(fill=Family),width=0.6,outlier.size=0) +
    #geom_jitter(width=0.5,aes(color=Query_environment)) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
        axis.text.x = element_text(size = 14, color = 'black',angle=60,hjust=1),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.title = element_text(size=16,color="black"),
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=14,color="black"),
        axis.ticks = element_blank(),
        strip.text = element_text(size=14,color="black"),
        plot.title = element_text(hjust=0.5,size=16)) +
    labs(x="Family",y="SV number per 1 Mb") +
    scale_x_discrete(limits=order_fs) +
    scale_fill_manual(values=color) +
    guides(fill=FALSE) +
    guides(color=FALSE) +
    geom_signif(comparisons=list(c("Asfuvirales","Caudovirales"))) +
    geom_signif(comparisons=list(c("Imitervirales","Caudovirales")),y_position=220) +
    geom_signif(comparisons=list(c("Imitervirales","uc_Heunggongvirae")),y_position=240)
dev.off()


# Interactive range
color.categories = c(
    "Species"="#389E39",
    "Genus" = "#CCDD8D",
    "Family"="#2577AC",
    "Order" = "#B8D2E2",
    "Class" = "#CD2B20",
    "Phylum" = "#F2C173",
    "Kingdom"="#ECA1A5"
)

data_plot = data.frame(
    Count = c(88,10,28,5,10,0,31),
    Range = c("Species","Genus","Family","Order","Class","Phylum","Kingdom")
)
data_plot$Range = factor(data_plot$Range,levels=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
data_plot$Direction='x'

pdf("/Users/laisenying/Desktop/FigS5e.pdf",width=8,height=5.5)
ggplot(data_plot,aes(x=Direction,y=Count)) +
    theme_bw() +
    geom_bar(aes(fill=Range),stat='identity',position='fill',color='black',width=0.8) +
    theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size=14,color='black',angle=60,hjust=1),
    axis.text.y = element_text(size=14,color='black'),
    axis.title = element_text(size=16,color="black"),
    plot.title = element_text(hjust=0.5,size=18),
    panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
    legend.text = element_text(size=14,color="black"),
    legend.title = element_text(size=14,color="black"),
    axis.ticks = element_blank()
    ) +
    labs(x=" ",y="Proportion",fill="Interactive range") +
    scale_fill_manual(values=as.vector(color.categories),limits=names(color.categories)) +
    coord_polar(theta = 'y',start = pi )
dev.off()


# Identity
BV_data = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/BV_all_interactive_data.tsv",sep="\t")
BV_data = BV_data[BV_data$Viral_Source=='IMG/VR',]
BV_data = BV_data[BV_data$Bacterial_Source=='CH-Binning',]
BV_data = BV_data[order(-BV_data$identity),]
BV_data = BV_data[!duplicated(BV_data$NR_SVID),]

BV_data$Group='=100%'
BV_data$Group[BV_data$identity<100]='>99%'
BV_data$Group[BV_data$identity<100]='>99%'
BV_data$Group[BV_data$identity<95]='>90%'
BV_data$Group[BV_data$identity<90]='>85%'
BV_data$Group[BV_data$identity<85]='>80%'

plot_bar = data.frame(table(BV_data$Group))
plot_bar$Group = factor(plot_bar$Var1,levels=c("=100%",">99%",">95%",">90%",">85%",">80%"))

color <- c(brewer.pal(9,'Set1'))
p1=ggplot(plot_bar,aes(x=Group,y=Freq)) +
    geom_bar(aes(fill=Group),width=0.6,stat='identity',color='black') +
    geom_text(aes(label=Freq,y=Freq+30),size=5) +
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
        axis.ticks = element_blank()
    ) +
    labs(x="Sequence identity between phages and bacteria",y='Frequency',title='Hit to CHGB') +
    guides(fill=FALSE) +
    scale_fill_manual(values=color)


BV_data = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/BV_all_interactive_data.tsv",sep="\t")
BV_data = BV_data[BV_data$Viral_Source=='IMG/VR',]
BV_data = BV_data[BV_data$Bacterial_Source=='HumGut',]
BV_data = BV_data[order(-BV_data$identity),]
BV_data = BV_data[!duplicated(BV_data$NR_SVID),]

BV_data$Group='=100%'
BV_data$Group[BV_data$identity<100]='>99%'
BV_data$Group[BV_data$identity<100]='>99%'
BV_data$Group[BV_data$identity<95]='>90%'
BV_data$Group[BV_data$identity<90]='>85%'
BV_data$Group[BV_data$identity<85]='>80%'

plot_bar = data.frame(table(BV_data$Group))
plot_bar$Group = factor(plot_bar$Var1,levels=c("=100%",">99%",">95%",">90%",">85%",">80%"))

color <- c(brewer.pal(9,'Set1'))
p2=ggplot(plot_bar,aes(x=Group,y=Freq)) +
    geom_bar(aes(fill=Group),width=0.6,stat='identity',color='black') +
    geom_text(aes(label=Freq,y=Freq+30),size=5) +
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
        axis.ticks = element_blank()
    ) +
    labs(x="Sequence identity between phages and bacteria",y='Frequency',title='Hit to HumGut') +
    guides(fill=FALSE) +
    scale_fill_manual(values=color)


pdf("/Users/laisenying/Desktop/Fig. S6.pdf",width=12,height=6)
plot_grid(p1,p2,nrow=1)
dev.off()



##########################################################################
################# Construction of BV interactive network #################
##########################################################################

import pandas as pd
import numpy as np

BV_interactive_data = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_coevolution/Interactive_range/BV_all_interactive_data.tsv",sep="\t",index_col=0)
BV_interactive_data_CHGV = BV_interactive_data.loc[BV_interactive_data['Viral_Source']=='IMG/VR',]
BV_interactive_data_CHGV = BV_interactive_data_CHGV.loc[BV_interactive_data_CHGV['Direction']=='bacteria-to-viralSV',]
BV_interactive_data_CHGV_filter = BV_interactive_data_CHGV.drop_duplicates(['NR_SVID','Rep_genome'])
BV_interactive_data_CHGV_filter = BV_interactive_data_CHGV_filter.iloc[np.argsort(-BV_interactive_data_CHGV_filter['identity']),]
BV_interactive_data_CHGV_filter = BV_interactive_data_CHGV_filter.drop_duplicates(['NR_SVID'])

BV_interactive_data_CHGV_filter['Interaction_id'] = BV_interactive_data_CHGV_filter['Genus'] +"~"+ BV_interactive_data_CHGV_filter['Viral_contig']
pd.value_counts(BV_interactive_data_CHGV_filter['Interaction_id'])
BV_interactive_data_CHGV_filter.index = BV_interactive_data_CHGV_filter['Interaction_id']
BV_interactive_data_CHGV_filter['Interact_num'] = pd.value_counts(BV_interactive_data_CHGV_filter['Interaction_id'])[BV_interactive_data_CHGV_filter['Interaction_id']]
BV_interactive_data_CHGV_filter['Count'] = 1
grouped = BV_interactive_data_CHGV_filter.groupby(['Viral_contig'])
BV_interactive_data_CHGV_filter.index = BV_interactive_data_CHGV_filter['Viral_contig']
BV_interactive_data_CHGV_filter['Viral_hit_num'] = (grouped['Count'].sum())[BV_interactive_data_CHGV_filter['Viral_contig']]
BV_interactive_data_CHGV_filter['Interactive_weight'] = BV_interactive_data_CHGV_filter['Interact_num']/BV_interactive_data_CHGV_filter['Viral_hit_num']
BV_interactive_data_CHGV_filter2 = BV_interactive_data_CHGV_filter.drop_duplicates(['Interaction_id'])
BV_interactive_data_CHGV_filter2.index = range(BV_interactive_data_CHGV_filter2.shape[0])

selected_genus = pd.value_counts(BV_interactive_data_CHGV_filter2['Genus'])[pd.value_counts(BV_interactive_data_CHGV_filter2['Genus'])>=10].index
BV_interactive_data_CHGV_filter3 = BV_interactive_data_CHGV_filter2.iloc[[x in selected_genus for x in BV_interactive_data_CHGV_filter2['Genus']],]
BV_interactive_data_CHGV_filter3['Genus_name'] = BV_interactive_data_CHGV_filter3['Genus']

BV_interactive_data_CHGV_filter3.to_csv("/home1/Laisenying/Projects/PhageSV/Results/BV_Network_plot_IMGVR.tsv",sep="\t")

scp Laisenying@10.190.248.213:/home1/Laisenying/Projects/PhageSV/Results/BV_Network_plot_IMGVR.tsv '/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/BV_Network_plot_IMGVR.tsv'
