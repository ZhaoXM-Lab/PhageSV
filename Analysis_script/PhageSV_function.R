#####################
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(forcats)
library(ggpubr)
library(ggalluvial)
color <- c(brewer.pal(12,'Set3'))
c(color[3],color[5],color[6],color[4])

sv_colors1=c("#339933", "#336699", "#CCCC33", "#CC6633")
sv_colors2=c("#CC6666", "#999933", "#339999", "#996699")
sv_colors3 = c("#FFCCCC", "#99CCFF", "#FFCC99", "#CCCCFF", "#99CCCC", "#FFCC99", "#CCFFCC")
sv_colors4=c("#9999CC","#CC9999","#99CC99")
sv_colors5 = c("#336699","#CCCC33")
sv_colors6 = c("#CC6666","#6699CC")




#############################################################################
################# Viral SV associated functional category ###################
#############################################################################

#scp Laisenying@10.190.248.211:/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/SV_functional_enrichment/ViralSV/ViralSV_Manual_Category_greater_enriched.tsv /Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/ViralSV_Manual_Category_greater_enriched.tsv

KO_enrich = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/ViralSV_Manual_Category_greater_enriched.tsv",sep="\t")
KO_enrich$INSDEL_fdr = p.adjust(KO_enrich$INSDEL_pvalue,method='fdr')
KO_enrich$DUP_fdr = p.adjust(KO_enrich$DUP_pvalue,method='fdr')
KO_enrich$INV_fdr = p.adjust(KO_enrich$INV_pvalue,method='fdr')
color <- c(brewer.pal(12,'Set3'))
Level1_colors=c()
i=1
for(x in names(table(KO_enrich$Level1_category))){
    Level1_colors[x] = color[i]
    i= i + 1
}

INSDEL_KOs_data = KO_enrich[KO_enrich$INSDEL_fdr<0.05,]
INSDEL_KOs_data = INSDEL_KOs_data[order(INSDEL_KOs_data$INSDEL_foldchange),]
INSDEL_KOs_data = INSDEL_KOs_data[,c("INSDELRatio","INSDEL_foldchange","INSDEL_fdr","Level1_category","Level2_category")]
colnames(INSDEL_KOs_data) = c("GeneRatio","foldchange","fdr","Level1","Level2")

DUP_KOs_data = KO_enrich[KO_enrich$DUP_fdr<0.05,]
DUP_KOs_data = DUP_KOs_data[order(DUP_KOs_data$DUP_foldchange),]
DUP_KOs_data = DUP_KOs_data[,c("DUPRatio","DUP_foldchange","DUP_fdr","Level1_category","Level2_category")]
colnames(DUP_KOs_data) = c("GeneRatio","foldchange","fdr","Level1","Level2")


INV_KOs_data = KO_enrich[KO_enrich$INV_fdr<0.05,]
INV_KOs_data = INV_KOs_data[order(INV_KOs_data$INV_foldchange),]
INV_KOs_data = INV_KOs_data[,c("INVRatio","INV_foldchange","INV_fdr","Level1_category","Level2_category")]
colnames(INV_KOs_data) = c("GeneRatio","foldchange","fdr","Level1","Level2")

INSDEL_KOs_data$KO_id = paste('INSDEL',INSDEL_KOs_data$Level2,sep='') # 42
DUP_KOs_data$KO_id = paste('DUP',DUP_KOs_data$Level2,sep='') # 1
INV_KOs_data$KO_id = paste('INV',INV_KOs_data$Level2,sep='') # 26
INSDEL_KOs_data$SVTYPE='INS & DEL'
DUP_KOs_data$SVTYPE='DUP'
INV_KOs_data$SVTYPE='INV'
KO_enrich_plot = rbind(INV_KOs_data,DUP_KOs_data,INSDEL_KOs_data)
KO_enrich_plot$SVTYPE = factor(KO_enrich_plot$SVTYPE,levels=c("INV","DUP","INS & DEL"))
KO_enrich_plot = KO_enrich_plot[order(KO_enrich_plot$SVTYPE,KO_enrich_plot$Level1,-KO_enrich_plot$foldchange),]
KO_enrich_plot$Level2[KO_enrich_plot$Level2=='Transpotion & Transposase'] = 'Transposition & Transposase'
Color.category = c(
    "INS & DEL" = sv_colors1[1],
    "DUP" = sv_colors1[2],
    "INV" = sv_colors1[3]
)

HGT_associated=c('Conjugative DNA transfer','DNA mediated transformation','Integrase','Transposition & Transposase','Recombination & Recombinase')
DNA_Methylation=c('Methyltransferase','Methylase')
CRSIPR=c("CRISPR-cas system")
Antibiotic_resistance=c('Antibiotic resistance')

KO_enrich_plot$Text_color='black'
KO_enrich_plot[KO_enrich_plot$Level2 %in% HGT_associated,'Text_color'] = sv_colors2[1]
KO_enrich_plot[KO_enrich_plot$Level2 %in% DNA_Methylation,'Text_color'] = sv_colors2[2]
KO_enrich_plot[KO_enrich_plot$Level2 %in% CRSIPR,'Text_color'] = Level1_colors['CRISPR-cas system']
KO_enrich_plot[KO_enrich_plot$Level2 %in% Antibiotic_resistance,'Text_color'] = sv_colors2[3]


p1=ggplot(KO_enrich_plot,aes(x=log2(foldchange),y=KO_id)) +
    geom_point(aes(color=SVTYPE,size=-log10(fdr))) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size=14,color='black'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=16,color='black'),
        plot.title = element_text(size=18,color='black',hjust=0.5),
        panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
        legend.text = element_text(size=14,color='black'),
        legend.title = element_text(size=14,color='black'),
        axis.ticks.y = element_blank()
    ) +
    labs(x='log2(Fold change)',y='KEGG Orthology(KOs)',color=" ") +
    scale_y_discrete(limits=KO_enrich_plot$KO_id,labels=KO_enrich_plot$Level2) +
    scale_color_manual(values=as.vector(Color.category),limits=names(Color.category)) +
    guides(size=FALSE) +
    geom_hline(yintercept=16.5,linetype='dashed') +
    geom_hline(yintercept=24.5,linetype='dashed')

p2 = ggplot(KO_enrich_plot,aes(x=" ",y=KO_id)) +
    geom_tile(aes(fill=Level1))+
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size=14,color='black'),
        axis.text.y = element_text(size=14,color=KO_enrich_plot$Text_color),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.text = element_text(size=14,color='black'),
        legend.title = element_text(size=14,color='black'),
        legend.position='left'
    )  +
    labs(x=' ',y=' ',fill=' ')+
    scale_y_discrete(limits=KO_enrich_plot$KO_id,labels=KO_enrich_plot$Level2) +
    scale_fill_manual(values=as.vector(Level1_colors),limits=names(Level1_colors)) +
    theme(plot.margin = unit(c(0.17,-0.1,0.76,0.1), 'cm'))

pdf("/Users/laisenying/Desktop/Category_enrich_viral.pdf",width=11,height=12)
plot_grid(p2,p1,nrow=1,rel_widths=c(0.4,0.4))
dev.off()

#############################################################################
################# Five major protein families of recombinase ################
#############################################################################

scp Laisenying@10.190.248.211:/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/Lifestyle_analysis/HGV_Viral_SV_density_with_proMGE_Recombinase.tsv /Users/laisenying/Desktop/HGV_Viral_SV_density_with_proMGE_Recombinase.tsv

Recombinase_data = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/SV_density_with_proMGE_Recombinase.tsv",sep="\t",row.names=1)
Recombinase_data = Recombinase_data[!duplicated(Recombinase_data$Viral_contig),]

Recombinase_data$HUH='N'
Recombinase_data$HUH[Recombinase_data$HUH.count>0]='Y'
Recombinase_data$DDE='N'
Recombinase_data$DDE[Recombinase_data$DDE.count>0]='Y'
Recombinase_data$Ser='N'
Recombinase_data$Ser[Recombinase_data$Ser.count>0]='Y'
Recombinase_data$Tyr='N'
Recombinase_data$Tyr[Recombinase_data$Tyr.count>0]='Y'
Recombinase_data$Cas='N'
Recombinase_data$Cas[Recombinase_data$Cas.count>0]='Y'
Recombinase_data$Recombinase='N'
Recombinase_data$Recombinase[Recombinase_data$Recombinase_count>0]='Yes'
Recombinase_data$all_SV_density = Recombinase_data$all_SV_count/(Recombinase_data$Genome_len/1000000)
plot_data = Recombinase_data[,c("HUH","DDE","Ser","Tyr","Cas","all_SV_density")]

library(reshape2)
plot_data2 = melt(plot_data,id.vars=c("all_SV_density"))

pdf("/Users/laisenying/Desktop/SV_density_Recombinase.pdf",width=6,height=5)
ggplot(plot_data2,aes(x=variable,y=log(all_SV_density))) +
    geom_boxplot(aes(fill=value),width=0.6,outlier.size=0) +
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
        strip.text = element_text(size=14,color='black')
    ) +
    scale_fill_manual(values=c(color[3],color[7])) +
    labs(x="Recombinase category",y="log(Viral SV number per 1Mb)",fill='Carry gene') +
    geom_signif(comparisons=list(c("Y","N")))
dev.off()

#############################################################################
############ Distribution of likely sources of the phage genes ##############
#############################################################################
plot_data = data.frame(
Group=c("Background","INS","DEL","DUP","INV"),
    Unk=c(381033,4768,5031,920,567),
    Unc=c(3415,40,219,32,29),
    bac=c(126263,2532,4076,495,511),
    viral=c(225920,1195,3313,575,287)
)

fisher.test(rbind(c(2532,4768+40+1195),c(126263,225920+3415+381033))) #  INS: OR = 2.04, p-value < 2.2e-16
fisher.test(rbind(c(4076,5031+32+3313),c(126263,225920+3415+381033))) # DEL: OR = 2.35, p-value < 2.2e-16
fisher.test(rbind(c(495,920+32+575),c(126263,225920+3415+381033))) # DUP OR = 1.57  p-value < 2.2e-16
fisher.test(rbind(c(511,567+29+287),c(126263,225920+3415+381033))) # INV OR = 2.79 p-value < 2.2e-16



colnames(plot_data) = c("Group","Unknown function","Uncertain source","Bacteria","Viral")

plot_data[,c("Unknown function","Uncertain source","Bacteria","Viral")] = plot_data[,c("Unknown function","Uncertain source","Bacteria","Viral")]/rowSums(plot_data[,c("Unknown function","Uncertain source","Bacteria","Viral")])
df1 <- melt(plot_data,id.vars = 'Group',measure.vars = c("Unknown function","Uncertain source","Bacteria","Viral"))
names(df1)[1:2] <- c("X","group")

color <- c(brewer.pal(12,'Set3'))
category.colors <- c(
    "Uncertain source" = color[3],
    "Unknown function" = color[4],
    "Viral" = color[5],
    "Bacteria" = color[6]
)


df1$group = factor(df1$group,levels=c("Viral","Unknown function","Uncertain source","Bacteria"))
df1$Category='Viral SVs'

plot_data = data.frame(
    Group=c("Background","INS","DEL","DUP","INV"),
    Unk=c(16910,9517,5003,104,87),
    Unc=c(75,1729,23,0,0),
    bac=c(81894,27257,9270,168,228),
    viral=c(7529,5888,2029,24,67)
)


colnames(plot_data) = c("Group","Unknown function","Uncertain source","Bacteria","Viral")

plot_data[,c("Unknown function","Uncertain source","Bacteria","Viral")] = plot_data[,c("Unknown function","Uncertain source","Bacteria","Viral")]/rowSums(plot_data[,c("Unknown function","Uncertain source","Bacteria","Viral")])
df2 <- melt(plot_data,id.vars = 'Group',measure.vars = c("Unknown function","Uncertain source","Bacteria","Viral"))
names(df2)[1:2] <- c("X","group")
df2$group = factor(df2$group,levels=c("Bacteria","Uncertain source","Unknown function","Viral"))
df2$Category = 'Bacterial SVs'
df_plot = rbind(df1,df2)
library(ggalluvial)

pdf("/Users/laisenying/Desktop/Fig1.pdf",width=6,height=6)
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
    scale_x_discrete(limits=c("Background","INS","DEL","DUP","INV")) +
    labs(x=" ",y="Percentage",fill="Gene source",title='Viral SVs')
dev.off()

pdf("/Users/laisenying/Desktop/Fig2.pdf",width=6,height=6)
ggplot(df2, aes( x = X,y=100 * value,fill = group,
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
    scale_x_discrete(limits=c("Background","INS","DEL","DUP","INV")) +
    labs(x=" ",y="Percentage",fill="Gene source",title='Bacterial SVs')
dev.off()

