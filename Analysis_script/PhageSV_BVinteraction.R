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
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(forcats)
library(ggpubr)
library(factoextra)
library(cluster)
library(fpc)

color <- c(brewer.pal(12,'Set3'))
c(color[3],color[5],color[6],color[4])

sv_colors1=c("#339933", "#336699", "#CCCC33", "#CC6633")
sv_colors2=c("#CC6666", "#999933", "#339999", "#996699")
sv_colors3 = c("#FFCCCC", "#99CCFF", "#FFCC99", "#CCCCFF", "#99CCCC", "#FFCC99", "#CCFFCC")
sv_colors4=c("#9999CC","#CC9999","#99CC99")
sv_colors5 = c("#336699","#CCCC33")
sv_colors6 = c("#CC6666","#6699CC")


###############################################################################
#########1. PB interaction in phage-bacteria pair with GE exchange ############
###############################################################################
# 存在Genetic exchange的phage和bacteria之间其丰度也通常存在显著正相关
# Python3
from scipy.stats import pearsonr,spearmanr
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict

BV_all_interactive_data = pd.read_csv("/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_coevolution/Interactive_range/BV_all_interactive_data.tsv",sep="\t",index_col=0)
BV_all_interactive_data = BV_all_interactive_data.loc[BV_all_interactive_data['Viral_Source'] == 'HGV', ]

Bac_RPKM_profile = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/Bac_SV/Bac_abundance/RPKM_file/Bin_RPKM_inf.tsv",sep="\t",index_col=0)
Viral_RPKM_profile = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/ViralGenome/HGV_viral_RPKM_profile.tsv",sep="\t",index_col=0)
Bac_species_profile = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/ITS/data/NGS/BacSequencing/metaphlan_result/species_abundance_profile.csv",index_col=0)
shared_samples = set(Viral_RPKM_profile.columns).intersection(Bac_species_profile.columns)
Bac_species_profile_filter = Bac_species_profile.loc[:,shared_samples]
Viral_RPKM_profile_filter = Viral_RPKM_profile.loc[:,shared_samples]
Viral_RPKM_profile_filter = Viral_RPKM_profile_filter.fillna(0)

Bac_species_profile_filter['Genus'] = [x.split("|")[-2] for x in Bac_species_profile_filter.index]
Bac_species_profile_filter['Family'] = [x.split("|")[-3] for x in Bac_species_profile_filter.index]
BV_all_interactive_data.index = BV_all_interactive_data['Viral_contig']
BV_all_interactive_data = BV_all_interactive_data.loc[set(Viral_RPKM_profile.index).intersection(BV_all_interactive_data['Viral_contig']),]
BV_all_interactive_data['Genus'] = [x.split(';')[-2] for x in BV_all_interactive_data['Taxonomy']]
interaction_data = defaultdict(list)

for viral in set(BV_all_interactive_data['Viral_contig']):
    subdata = pd.DataFrame(BV_all_interactive_data.loc[viral,])
    if subdata.shape[1] == 1:
        subdata = subdata.T
    subdata = subdata.drop_duplicates(['Genus','NR_SVID','Direction'])
    for genus,count in zip(pd.value_counts(subdata['Genus']).index,pd.value_counts(subdata['Genus'])):
        if genus not in set(Bac_species_profile_filter['Genus']):
            genus = "_".join(genus.split("_")[:-1])
            if genus not in set(Bac_species_profile_filter['Genus']):
                continue
        bac_values = (Bac_species_profile_filter.loc[Bac_species_profile_filter['Genus']==genus,shared_samples]).sum()
        viral_values = Viral_RPKM_profile_filter.loc[viral,shared_samples]
        r,pvalue = spearmanr(bac_values,viral_values)
        interaction_data['Viral_contig'].append(viral)
        interaction_data['Bac_tax'].append(genus)
        interaction_data['pearson_r'].append(r)
        interaction_data['pvalue'].append(pvalue)
        interaction_data['Interact_count'].append(count)
        if (len(set(subdata['Genus']))==1):
            interaction_data['type'].append('Specialist')
        else:
            interaction_data['type'].append('Generalist')


interaction_inf = pd.DataFrame(interaction_data)
interaction_inf = interaction_inf.loc[interaction_inf['pearson_r']==interaction_inf['pearson_r'],]
interaction_inf = interaction_inf.iloc[np.argsort(-interaction_inf['Interact_count']),]
interaction_inf = interaction_inf.drop_duplicates(['Viral_contig','Bac_tax'])
interaction_inf.to_csv("/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_interaction/Phage_bac_interact_correlation.tsv",sep="\t")

interaction_inf = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_interaction/Phage_bac_interact_correlation.tsv",sep="\t")


# R plot
scp Laisenying@10.190.248.211:/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_interaction/Phage_bac_interact_correlation.tsv /Users/laisenying/Desktop/Phage_bac_interact_correlation.tsv


data = read.csv("/Users/laisenying/Desktop/Phage_bac_interact_correlation.tsv",sep="\t",row.names=1)
data[data$Interact_count>10,'type'] = 'Count > 10'
data[data$Interact_count<10 & data$Interact_count != 0,'type'] = '0 < count < 10'
data[data$Interact_count==0,'type'] = 'Count = 0'

data$type = factor(data$type,levels=c("Count > 10","0 < count < 10","Count = 0"),labels=c("Count > 10","0 < count < 10","Count = 0"))
data = data[!is.na(data$type),]

pdf("/Users/laisenying/Desktop/GV_correlation.pdf",width=11,height=4.5)
ggplot(data,aes(x=pearson_r)) +
    geom_density(aes(fill=type,group=type),alpha=0.4) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_text(size=14,color='black'),
        axis.title = element_text(size=16),
        legend.text = element_text(size=14,color='black'),
        legend.title = element_text(size=14,color='black'),
        axis.ticks = element_blank()
    ) +
    scale_fill_manual(values=c("#CC6666", "#999933", "#339999")) +
    labs(x='Phage-Bacteria correlation',y='Density',fill="GE-SV num") +
    geom_vline(xintercept=mean(data[data$type=='Count = 0','pearson_r'],na.rm=TRUE),color='#339999',linetype='dashed') +
    geom_vline(xintercept=mean(data[data$type=='0 < count < 10','pearson_r'],na.rm=TRUE),color='#999933',linetype='dashed') +
    geom_vline(xintercept=mean(data[data$type=='Count > 10','pearson_r'],na.rm=TRUE),color='#CC6666',linetype='dashed')
dev.off()

sum(data[data$type=='Count > 10','pvalue']<0.05)/sum(data$type=='Count > 10')
sum(data[data$type=='0 < count < 10','pvalue']<0.05)/sum(data$type=='0 < count < 10')
sum(data[data$type=='Count = 0','pvalue']<0.05)/sum(data$type=='Count = 0')

sum(data[data$type=='Count > 10' & data$pvalue<0.05,'pearson_r']>0)/sum(data[data$type=='Count > 10','pvalue']<0.05)
sum(data[data$type=='0 < count < 10' & data$pvalue<0.05,'pearson_r']>0)/sum(data[data$type=='0 < count < 10','pvalue']<0.05)
sum(data[data$type=='Count = 0' & data$pvalue<0.05,'pearson_r']>0)/sum(data[data$type=='Count = 0','pvalue']<0.05)


plot_bar = data.frame(
    Specialist=c(sum(data[data$type=='Count > 10','pvalue']<0.05),sum(data$type=='Count > 10')-sum(data[data$type=='Count > 10','pvalue']<0.05)),
    Generalist =  c(sum(data[data$type=='0 < count < 10','pvalue']<0.05),sum(data$type=='0 < count < 10')-sum(data[data$type=='0 < count < 10','pvalue']<0.05)),
    All = c(sum(data[data$type=='Count = 0','pvalue']<0.05),sum(data$type=='Count = 0')-sum(data[data$type=='Count = 0','pvalue']<0.05))
)
plot_bar$Sig = c("Yes","No")
library(reshape2)
plot_bar2 = melt(plot_bar)
plot_bar2$variable = factor(plot_bar2$variable,levels=c("Specialist","Generalist","All"),labels=c("Count > 10","0 < count < 10","Count = 0"))
plot_bar2$Sig = factor(plot_bar2$Sig,levels=c("No","Yes"),labels=c("p>0.05","p<0.05"))
pdf("/Users/laisenying/Desktop/GV_correlation_bar.pdf",width=3,height=3)
ggplot(plot_bar2,aes(x=variable,y=value)) +
    geom_bar(aes(fill=Sig),stat='identity',position='fill',color='black',width=0.6) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size=14,color='black',angle=60,hjust=1),
        axis.text.y = element_text(size=14,color='black'),
        axis.title = element_text(size=16),
        legend.text = element_text(size=14,color='black')
    ) +
    scale_fill_manual(values=c(color[5],color[4])) +
    labs(x=" ",y="Proportion",fill=" ") +
    scale_x_discrete(limits=c("Count = 0","0 < count < 10","Count > 10"),labels=c("=0","1~10",">10"))
dev.off()



###############################################################################
#######################2. Strain-level correlation ############################
###############################################################################


scp Laisenying@10.190.248.213:/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_interaction/Phage_bac_interact_correlation.tsv '/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/Phage_bac_interact_correlation.tsv'
scp Laisenying@10.190.248.213:/share/home1/Laisenying/Data-analysis/projects/ITS/data/NGS/BacSequencing/metaphlan_result/species_abundance_profile.csv '/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/species_abundance_profile.csv'
scp Laisenying@10.190.248.211:/share/home1/Laisenying/Data-analysis/projects/PhageSV/PacBio_SV/CAST_merge/Population/Merge/HGV_SV_profile.tsv '/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/CHGV_SV_profile.tsv'
scp Laisenying@10.190.248.211:/share/home1/Laisenying/Data-analysis/projects/PhageSV/ViralGenome/sample_inf.tsv '/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/sample_inf.tsv'
scp Laisenying@10.190.248.211:/share/home1/Laisenying/Data-analysis/projects/PhageSV/ViralGenome/HGV_viral_RPKM_profile.tsv '/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/HGV_viral_RPKM_profile.tsv'
scp Laisenying@10.190.248.211:/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_interaction/SV_bacdiversity/Viral_SV_bacdiversity.tsv '/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/Viral_SV_bacdiversity.tsv'



BV_data = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/BV_all_interactive_data.tsv",sep="\t")
BV_data = BV_data[BV_data$Viral_Source=='HGV',]
BV_data = BV_data[BV_data$Bacterial_Source=='CH-Binning',]
Species_abun_data = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/species_abundance_profile.csv",row.names=1)
interaction_inf = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/Phage_bac_interact_correlation.tsv",sep="\t")
Viral_geno_profile = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/CHGV_SV_profile.tsv",sep="\t",row.names=1) # Viral SV profile
sample_inf = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/sample_inf.tsv",sep="\t",row.names=1)
Viral_RPKM = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/HGV_viral_RPKM_profile.tsv",sep="\t",row.names=1)

SV_bacdiversity = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/Viral_SV_bacdiversity.tsv",sep='\t')

rownames(sample_inf) = sample_inf$Run
cor.test(sample_inf[colnames(Viral_RPKM),'Bac_Shannon'],as.numeric(Viral_RPKM['SK05_contig-120_1',]),na.rm=TRUE)



# 能够对Bacterial community产生显著影响的SV
Viral_geno_profile = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/CHGV_SV_profile.tsv",sep="\t",row.names=1) # Viral SV profile
Viral_geno_profile = Viral_geno_profile[(rowSums((Viral_geno_profile=='0/1') | (Viral_geno_profile=='1/1'))>=5)&(rowSums(Viral_geno_profile=='0/0')>=5),]
SV_bac_stat = data.frame(Viral_contig=c(),SVID=c(),bray_pvalue=c(),bray_R2=c(),jaccard_pvalue=c(),jaccard_R2=c())
for(SV in rownames(Viral_geno_profile)){
    SV_inf = data.frame(Viral_geno_profile[SV,])
    contig = SV_inf$Chr
    SV_inf  = SV_inf[,!(colnames(SV_inf) %in% c("Chr","Start","End","SVTYPE","SVLEN","SVID.1"))]
    SV_inf = data.frame(t(SV_inf))
    SV_inf$sample = rownames(SV_inf)
    SV_inf = data.frame(SV_inf[SV_inf[SV]!='Noexist',])
    colnames(SV_inf) = c('Genotype','sample')
    if(nrow(SV_inf)<10){next;}
    shared_samples = SV_inf$sample[SV_inf$sample %in% colnames(Bac_genus_profile)]
    SV_inf = SV_inf[shared_samples,]
    dist_bac_bray = vegdist(t(Bac_genus_profile[,shared_samples]),method='bray')
    dist_bac_bray[is.na(dist_bac_bray)] = 0
    if(length(names(table(SV_inf$Genotype)))==1){next;}
    SV_adnois_bray = adonis(dist_bac_bray~Genotype,data = SV_inf,permutations = 999)
    bray_pvalue = SV_adnois_bray$aov.tab['Genotype','Pr(>F)']
    bray_R2 = SV_adnois_bray$aov.tab['Genotype','R2']
    dist_bac_jaccard = vegdist(t(Bac_genus_profile[,shared_samples]),method='jaccard')
    dist_bac_jaccard[is.na(dist_bac_jaccard)] = 0
    SV_adnois_jaccard = adonis(dist_bac_jaccard~Genotype,data = SV_inf,permutations = 999)
    jaccard_pvalue = SV_adnois_jaccard$aov.tab['Genotype','Pr(>F)']
    jaccard_R2 = SV_adnois_jaccard$aov.tab['Genotype','R2']
    SV_bac_stat = rbind(SV_bac_stat,data.frame(Viral_contig=contig,SVID=SV,bray_pvalue=bray_pvalue,bray_R2=bray_R2,jaccard_pvalue=jaccard_pvalue,jaccard_R2=jaccard_R2))
    
}
for(contig in unique(SV_bac_stat$Viral_contig)){
    SV_bac_stat[SV_bac_stat$Viral_contig==contig,'bray_fdr'] = p.adjust(SV_bac_stat[SV_bac_stat$Viral_contig==contig,'bray_pvalue'],method='fdr')
    SV_bac_stat[SV_bac_stat$Viral_contig==contig,'jaccard_fdr'] = p.adjust(SV_bac_stat[SV_bac_stat$Viral_contig==contig,'jaccard_pvalue'],method='fdr')
}
write.table(SV_bac_stat,"/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/viralSV_BacComm_adonis.tsv",sep="\t")



#######################################------------------#######################################
### Summary
#BV_comnunity_cor_data = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/SVprofile_BacComm_cor.tsv",sep="\t",row.names=1) # 与Community存在显著关联viral SV profile
#Bac_ViralSV_stat = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/viralSV_BacTaxa_adonis.tsv",sep="\t",row.names=1) # 与Bac taxa存在显著关联的Viral SV profile
SV_bacCom_stat = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/viralSV_BacComm_adonis.tsv",sep="\t",row.names=1)
#SVprofile_diversity_adonis = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/SVprofile_bacdiversity_adonis.tsv",sep="\t")
#SV_bacDiversity_stat = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/SVID_correlate_BacDiversity.tsv",sep="\t")


SV_function = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/ViralSV_gene_functional_category.tsv",sep="\t")

BacCom_asso_SVs = SV_bacCom_stat[(SV_bacCom_stat$bray_fdr<0.05) | (SV_bacCom_stat$jaccard_fdr<0.05),'SVID']
SV_function[SV_function$SVID %in% BacCom_asso_SVs,]

# DEL009658SUR: Phage_integrase, HicA_toxin, Eco57I,ResIII,T5orf172  ,Eco57I
# DEL008746SUR: Resolvase, N terminal domain Recombinase impB/mucB/samB family Glycosyl transferases group 1
# DEL003101SUR Glycosyltransferase family 36
# DEL003193SUR  RecX family GO:regulation of DNA repair
# DEL0010961SUR zur, furB; Fur family transcriptional regulator, ferric uptake regulator
# INS008252SUR Integrase core domain
# DEL001355SUR stage V sporulation protein G
# DEL001704SUR Belongs to the glycosyl hydrolase 43 family
# DEL003208SUR robable diguanylate cyclase DgcQ diguanylate cyclase
# DEL0010853SUR Unknown
# DEL003087SUR Unknown
# INS003570SUR Unknown
# INS009420SUR Unknown_function
# INS003618SUR Branched-chain amino acid transport protein (AzlD)
# DEL003629SUR Flg_new

## Visualization
SV_bacCom_stat[SV_bacCom_stat$bray_fdr<0.05,]
Viral_geno_profile = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/CHGV_SV_profile.tsv",sep="\t",row.names=1) # Viral SV profile

#pc = list()
#i=1
#for(SVID in unique(SV_bacCom_stat[SV_bacCom_stat$bray_pvalue<0.05,'SVID'])){
SVID='DEL003193SUR'
Geno_inf = data.frame(t(Viral_geno_profile[SVID,!(colnames(Viral_geno_profile) %in% c("Chr","Start","End","SVTYPE","SVLEN","SVID.1"))]))
Geno_inf$sample = rownames(Geno_inf)
colnames(Geno_inf) = c("Genotype","sample")
Geno_inf = Geno_inf[Geno_inf$Genotype!='Noexist',]

bac.dist = vegdist(t(Bac_genus_profile),method='bray')
obs.pcoa = dudi.pco(bac.dist, scannf = F, nf = 3)
Geno_inf$Bac_X = obs.pcoa$li[Geno_inf$sample,'A1']
Geno_inf$Bac_Y = obs.pcoa$li[Geno_inf$sample,'A2']

geno.colors=c(
    "0/0"=color[3],
    "0/1"=color[4],
    "1/1"=color[5]
)

pdf("/Users/laisenying/Desktop/DEL008746SUR_baccomm.pdf",width=6,height=5)
ggplot(Geno_inf,aes(x=Bac_X,y=Bac_Y)) +
    geom_point(aes(color=Genotype),size=2) +
    #stat_ellipse(aes(group=Genotype,color=Genotype),level = 0.8,size=0.8,alpha=0.8) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size=14,color='black'),
        axis.text.y = element_text(size=14,color='black'),
        axis.title = element_text(size=16,color='black'),
        strip.text = element_text(size=14,color='black'),
        plot.title = element_text(size=16,color='black',hjust=0.5),
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=14,color="black")
    ) +
    scale_color_manual(values=geno.colors,limits=names(geno.colors)) +
    labs(x="PCo1",y="PCo2",title=SVID,color='Genotype')
dev.off()






## Association between viral SV profiles and bacterial community

all_viral_contigs = unique(Viral_geno_profile$Chr)
Viral_geno_profile = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/CHGV_SV_profile.tsv",sep="\t",row.names=1) # Viral SV profile
pic = list()
i=1
ViralSVprofile_bacCom = data.frame(Viral_contig=c(),pvalue=c(),R2=c())
Viral_RPKM = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/HGV_viral_RPKM_profile.tsv",sep="\t",row.names=1)

RPKM_pc = list()
j=1
for(contig in all_viral_contigs){
    subdata = Viral_geno_profile[Viral_geno_profile$Chr==contig,]
    selected_samples = colnames(subdata)[colSums(subdata=='Noexist')==0][!(colnames(subdata)[colSums(subdata=='Noexist')==0] %in% c("Chr","Start","End","SVTYPE","SVLEN","SVID.1"))]
    if(length(selected_samples)<10){next;}
    subdata = subdata[,colnames(subdata) %in% selected_samples]
    subdata[subdata == '0/0'] = 0
    subdata[subdata == '0/1'] = 1
    subdata[subdata == '1/1'] = 2
    if(nrow(subdata)<3){next;}
    subdata2 = as.data.frame(apply(subdata,2,as.numeric))
    rownames(subdata2) = rownames(subdata)
    subdata2 = data.frame(t(subdata2))
    sample_d <- dist(subdata2, method = "canberra")
    sample_d[is.na(sample_d)] = 0
    pamk.best = pamk(sample_d,krange=1:3)
    if(pamk.best$nc==1){
        print(1)
        next
    }
    optimal_k = pamk.best$nc
    fit_sample <- hclust(sample_d, method="ward.D")
    data.cluster = cutree(fit_sample, optimal_k)
    cluster_inf = data.frame(data.cluster)
    colnames(cluster_inf) = c('cluster')
    cluster_inf$sample = rownames(cluster_inf)
    bac.dist = vegdist(t(Bac_genus_profile),method='bray')
    obs.pcoa = dudi.pco(bac.dist, scannf = F, nf = 3)
    cluster_inf$Bac_X = obs.pcoa$li[cluster_inf$sample,'A1']
    cluster_inf$Bac_Y = obs.pcoa$li[cluster_inf$sample,'A2']
    cluster_inf$Strain = paste('Strain',cluster_inf$cluster,sep=' ')
    cluster_inf = cluster_inf[!is.na(cluster_inf$Bac_X),]
    if(length(table(cluster_inf$cluster))==1){
        next;
    }
    if(sum(cluster_inf$cluster!=1)<5){
        next;
    }
    if(sum(cluster_inf$cluster!=2)<5){
        next;
    }
    if(sum(cluster_inf$cluster!=3)<5){
        next;
    }
    if((sum(cluster_inf$cluster==1)<5) & (sum(cluster_inf$cluster==2)<5)){
        next;
    }
    if((sum(cluster_inf$cluster==1)<5) & (sum(cluster_inf$cluster==3)<5)){
        next;
    }
    if((sum(cluster_inf$cluster==2)<5) & (sum(cluster_inf$cluster==3)<5)){
        next;
    }
    cluster_inf$Viral_RPKM = as.numeric(Viral_RPKM[contig,cluster_inf$sample])
    cluster_inf$Viral_RPKM[is.na(cluster_inf$Viral_RPKM)] = 0
    kt = kruskal.test(cluster_inf$Viral_RPKM,cluster_inf$Strain,na.rm=TRUE)
    if(kt$p.value<0.05){
        RPKM_pc[[j]] = ggplot(cluster_inf,aes(x=Strain,y=Viral_RPKM)) +
            #geom_violin(aes(color=Strain),fill='white') +
            geom_boxplot(aes(fill=Strain),width=0.6) +
            theme_bw() +
            theme(
                panel.grid = element_blank(),
                axis.text = element_text(size=14,color='black'),
                axis.title = element_text(size=16,color='black'),
                plot.title = element_text(size=18,color='black',hjust=0.5),
                strip.text = element_text(size=14,color='black'),
                legend.text = element_text(size=14,color="black"),
                legend.title = element_text(size=14,color="black")
            ) +
            scale_fill_manual(values=c(color[4],color[5],color[6])) +
            scale_color_manual(values=c(color[4],color[5],color[6])) +
            guides(fill=FALSE) +
            guides(color=FALSE) +
            geom_signif(comparisons=list(c("Strain 1","Strain 2"))) +
            geom_signif(comparisons=list(c("Strain 2","Strain 3"))) +
            geom_signif(comparisons=list(c("Strain 1","Strain 3")),y_position=1.1*max(cluster_inf$Viral_RPKM)) +
            labs(x=" ",y="Viral RPKM",title=contig)
       j = j + 1
            
    }
    Bac_genus_profile_filter = Bac_genus_profile[,rownames(cluster_inf)]
    bac.dist = vegdist(t(Bac_genus_profile_filter),method='bray')
    ad_result = adonis(bac.dist~Strain,cluster_inf)
    pvalue = ad_result$aov.tab['Strain','Pr(>F)']
    R2 = ad_result$aov.tab['Strain','R2']
    table(cluster_inf$cluster)
    if(pvalue<0.05){
        pic[[i]] = ggplot(cluster_inf,aes(x=Bac_X,y=Bac_Y)) +
            geom_point(aes(color=Strain)) +
            stat_ellipse(aes(group=Strain,color=Strain),level = 0.8,size=0.8,alpha=0.8) +
            theme_bw() +
            theme(
                panel.grid = element_blank(),
                axis.text.x = element_text(size=14,color='black'),
                axis.text.y = element_text(size=14,color='black'),
                axis.title = element_text(size=16,color='black'),
                strip.text = element_text(size=14,color='black'),
                plot.title = element_text(size=16,color='black',hjust=0.5),
                legend.text = element_text(size=14,color="black"),
                legend.title = element_text(size=14,color="black")
            ) +
            scale_color_manual(values=c(color[4],color[5],color[6])) +
            labs(x="PCo1",y="PCo2",color=" ",title=contig)
            i = i + 1
    }
    ViralSVprofile_bacCom = rbind(ViralSVprofile_bacCom,data.frame(Viral_contig=contig,pvalue=pvalue,R2=R2,RPKM_pvalue=kt$p.value))
}
ViralSVprofile_bacCom$RPKM_fdr = p.adjust(ViralSVprofile_bacCom$RPKM_pvalue,method='fdr')

write.table(ViralSVprofile_bacCom[ViralSVprofile_bacCom$RPKM_fdr<0.05,],"/Users/laisenying/Desktop/Table S1.tsv",sep="\t")

pdf("/Users/laisenying/Desktop/contig_PCoA.pdf",width=6,height=5)
for(i in 1:12){
    plot(pic[[i]])
}
dev.off()


