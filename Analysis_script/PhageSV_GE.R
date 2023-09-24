#####################
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(forcats)
library(ggpubr)
library(ggsignif)
library(reshape2)
library(ggalluvial)

color <- c(brewer.pal(12,'Set3'))
c(color[3],color[5],color[6],color[4])

sv_colors1=c("#339933", "#336699", "#CCCC33", "#CC6633")
sv_colors2=c("#CC6666", "#999933", "#339999", "#996699")
sv_colors3 = c("#FFCCCC", "#99CCFF", "#FFCC99", "#CCCCFF", "#99CCCC", "#FFCC99", "#CCFFCC")
sv_colors4=c("#9999CC","#CC9999","#99CC99")
sv_colors5 = c("#336699","#CCCC33")
sv_colors6 = c("#CC6666","#6699CC")


###############################################################
#### Distributions of likely sources of phage SV sequences ####
###############################################################

# SV sequences
# /share/home1/Laisenying/Data-analysis/projects/PhageSV/PacBio_SV/CAST_merge/Population/Result/SV_sequence/SV_200_all_sequence.fasta
# /share/home1/Laisenying/Data-analysis/projects/PhageSV/Bac_SV/Bacterial_genomes/HumGut/HumGut_HQ_genome.fasta
# /share/home1/Laisenying/Data-analysis/projects/PhageSV/DataBase/IMGVR/IMGVR_all_nucleotides_HQ.fna
# /share/home1/Laisenying/Data-analysis/projects/PhageSV/Bac_SV/Bacterial_genomes/all_bacterial_genomes.fna
# /share/home1/Laisenying/Data-analysis/projects/PhageSV/DataBase/IMGVR/SV_result/SV_sequence/IMGVR_Viral_SV_sequence.fasta

# Viral SV hit HumGut
/home1/Laisenying/miniconda3/bin/makeblastdb \
    -dbtype nucl \
    -in /share/home1/Laisenying/Tools/data/HumGut/all_HumGut_genome.fna \
    -input_type fasta \
    -parse_seqids \
    -out /home1/Laisenying/Projects/PhageSV/HumGut/all_HumGut_genome

/home1/Laisenying/miniconda3/bin/makembindex -iformat blastdb -input /home1/Laisenying/Projects/PhageSV/HumGut/all_HumGut_genome

Viral_SV=/home1/Laisenying/Projects/PhageSV/SV_sequence/CHGV_SV_200_all_sequence.fasta
~/miniconda3/bin/blastn -task megablast -db /home1/Laisenying/Projects/PhageSV/HumGut/all_HumGut_genome \
    -query ${Viral_SV} \
    -evalue 1e-5 \
    -outfmt '6 std qlen slen' \
    -num_threads 4 \
    -use_index true \
    -out /home1/Laisenying/Projects/PhageSV/SV_sequence/ViralSV_blastn_HumGut.out

# Viral SV hit CHGB
Viral_SV=/home1/Laisenying/Projects/PhageSV/SV_sequence/CHGV_SV_200_all_sequence.fasta
Bac_genome=/home1/Laisenying/Data-analysis/projects/PhageSV/Bac_SV/Bacterial_genomes/all_bacterial_genomes.fna

~/miniconda3/bin/blastn -task megablast -db /share/home1/Laisenying/Data-analysis/projects/PhageSV/Bac_SV/Bacterial_genomes/all_bacterial_genomes \
    -query ${Viral_SV} \
    -evalue 1e-5 \
    -outfmt '6 std qlen slen' \
    -num_threads 4 \
    -use_index true \
    -out /home1/Laisenying/Projects/PhageSV/SV_sequence/ViralSV_blastn_CHGB.out


# IMGVR SV hit HumGut and CHGB
Viral_SV=/home1/Laisenying/Projects/PhageSV/SV_sequence/IMGVR_Viral_SV_sequence.fasta
~/miniconda3/bin/blastn -task megablast -db /home1/Laisenying/Projects/PhageSV/HumGut/all_HumGut_genome \
    -query ${Viral_SV} \
    -evalue 1e-5 \
    -outfmt '6 std qlen slen' \
    -num_threads 4 \
    -use_index true \
    -out /home1/Laisenying/Projects/PhageSV/SV_sequence/IMGVR_blastn_HumGut.out

Viral_SV=/home1/Laisenying/Projects/PhageSV/SV_sequence/IMGVR_Viral_SV_sequence.fasta
~/miniconda3/bin/blastn -task megablast -db /share/home1/Laisenying/Data-analysis/projects/PhageSV/Bac_SV/Bacterial_genomes/all_bacterial_genomes \
    -query ${Viral_SV} \
    -evalue 1e-5 \
    -outfmt '6 std qlen slen' \
    -num_threads 4 \
    -use_index true \
    -out /home1/Laisenying/Projects/PhageSV/SV_sequence/IMGVR_blastn_CHGB.out


###############################################################

BV_interactive_data = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_coevolution/Interactive_range/BV_all_interactive_data.tsv",sep="\t",index_col=0)

HGVSV_inf = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/PacBio_SV/CAST_merge/Population/Result/SV_bed_inf.tsv",sep="\t",index_col=0)
IMGVRSV_inf = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/DataBase/IMGVR/SV_result/SV_sequence/IMGVR_Viral_SV_sequence_inf.tsv",sep="\t",index_col=0)


# CHGV Viral SV
NoINS_number = 1406
NoDEL_number = 1532
NoDUP_number = 111
NoINV_number = 116

sharedhit_INS_number = 545
sharedhit_DEL_number = 679
sharedhit_DUP_number = 68
sharedhit_INV_number = 50

CHGBhit_INS_number = 377
CHGBhit_DEL_number = 444
CHGBhit_DUP_number = 26
CHGBhit_INV_number = 31

HumGuthit_INS_number = 380
HumGuthit_DEL_number = 448
HumGuthit_DUP_number = 40
HumGuthit_INV_number = 20


data_plot_viral = data.frame(
    SVTYPE = c("INS","DEL","DUP","INV","INS","DEL","DUP","INV","INS","DEL","DUP","INV","INS","DEL","DUP","INV"),
    Group = c("Both HumGut and CHGB","Both HumGut and CHGB","Both HumGut and CHGB","Both HumGut and CHGB","CHGB","CHGB","CHGB","CHGB",
    "HumGut","HumGut","HumGut","HumGut","No hit","No hit","No hit","No hit"),
    Number=c(sharedhit_INS_number,sharedhit_DEL_number,sharedhit_DUP_number,sharedhit_INV_number,
    CHGBhit_INS_number,CHGBhit_DEL_number,CHGBhit_DUP_number,CHGBhit_INV_number,HumGuthit_INS_number,HumGuthit_DEL_number,
    HumGuthit_DUP_number,HumGuthit_INV_number,NoINS_number,NoDEL_number,NoDUP_number,NoINV_number)
)

data_plot_viral$SVTYPE = factor(data_plot_viral$SVTYPE,levels=c("INS","DEL","DUP","INV"))
data_plot_viral$Group = factor(data_plot_viral$Group,levels=c("No hit","HumGut","Both HumGut and CHGB","CHGB"))

Color.categorys = c(
    "No hit" = color[3],
    "Both HumGut and CHGB" = color[1],
    "HumGut" = color[2],
    "CHGB" = color[6]
)

p1 = ggplot(data_plot_viral,aes(x=SVTYPE,y=Number)) +
    geom_bar(aes(fill=Group),stat='identity',color='black',position='fill',width=0.8) +
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
    )  +
    scale_fill_manual(values=as.vector(Color.categorys),limits=names(Color.categorys)) +
    labs(x=" ",y="Proportion",fill="Source of phage SVs",title="CHGV")

## IMG/VR viral SVs
NoINS_number = 6583
NoDEL_number = 7483
NoDUP_number = 91
NoINV_number = 6583

sharedhit_INS_number = 199
sharedhit_DEL_number = 456
sharedhit_DUP_number = 17
sharedhit_INV_number = 0

CHGBhit_INS_number = 306
CHGBhit_DEL_number = 429
CHGBhit_DUP_number = 6
CHGBhit_INV_number = 4

HumGuthit_INS_number = 492
HumGuthit_DEL_number = 1059
HumGuthit_DUP_number = 43
HumGuthit_INV_number = 5


data_plot_viral = data.frame(
    SVTYPE = c("INS","DEL","DUP","INV","INS","DEL","DUP","INV","INS","DEL","DUP","INV","INS","DEL","DUP","INV"),
    Group = c("Both HumGut and CHGB","Both HumGut and CHGB","Both HumGut and CHGB","Both HumGut and CHGB","CHGB","CHGB","CHGB","CHGB",
    "HumGut","HumGut","HumGut","HumGut","No hit","No hit","No hit","No hit"),
    Number=c(sharedhit_INS_number,sharedhit_DEL_number,sharedhit_DUP_number,sharedhit_INV_number,
    CHGBhit_INS_number,CHGBhit_DEL_number,CHGBhit_DUP_number,CHGBhit_INV_number,HumGuthit_INS_number,HumGuthit_DEL_number,
    HumGuthit_DUP_number,HumGuthit_INV_number,NoINS_number,NoDEL_number,NoDUP_number,NoINV_number)
)

data_plot_viral$SVTYPE = factor(data_plot_viral$SVTYPE,levels=c("INS","DEL","DUP","INV"))
data_plot_viral$Group = factor(data_plot_viral$Group,levels=c("No hit","HumGut","Both HumGut and CHGB","CHGB"))

Color.categorys = c(
    "No hit" = color[3],
    "Both HumGut and CHGB" = color[1],
    "HumGut" = color[2],
    "CHGB" = color[6]
)

p3 = ggplot(data_plot_viral,aes(x=SVTYPE,y=Number)) +
    geom_bar(aes(fill=Group),stat='identity',color='black',position='fill',width=0.8) +
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
    )  +
    scale_fill_manual(values=as.vector(Color.categorys),limits=names(Color.categorys)) +
    labs(x=" ",y="Proportion",fill="Source of phage SVs",title="IMG/VR")


pdf("/Users/laisenying/Desktop/IMGVR_SV_sources.pdf",width=6,height=5)
p3
dev.off()


##########################################################################
#### Genetic divergence of viral SV sequences against CHGB and HumGut ####
##########################################################################

BV_interactive_data = pd.read_csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/BV_all_interactive_data.tsv",sep="\t",index_col=0)
BV_interactive_data_HGV = BV_interactive_data.loc[BV_interactive_data['Viral_Source']=='HGV',]
BV_interactive_data_HGV.loc[BV_interactive_data_HGV['Bacterial_Source']=='CH-Binning',]
BV_interactive_data_HGV.loc[BV_interactive_data_HGV['Bacterial_Source']=='HumGut',]

BV_interactive_data = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/BV_all_interactive_data.tsv",sep="\t")

color <- c(brewer.pal(12,'Set3'))
BV_interactive_data_HGV = BV_interactive_data[BV_interactive_data$Viral_Source=='HGV',]

BV_interactive_data_HGV$Bacterial_Source = factor(BV_interactive_data_HGV$Bacterial_Source,levels=c("CH-Binning","HumGut"),labels=c("CHGB","HumGut"))

pdf("/Users/laisenying/Desktop/GeneDiver_SV.pdf",width=4,height=5.2)
ggplot(BV_interactive_data_HGV,aes(x=Bacterial_Source,y=100 - identity)) +
    geom_violin(aes(fill=Bacterial_Source)) +
    geom_boxplot(aes(fill=Bacterial_Source),color='white',outlier.size=0,width=0.08) +
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
    labs(x=" ",y="Genetic divergence(%)") +
    guides(fill=FALSE) +
    geom_signif(comparisons=list(c("CHGB","HumGut")),y_position=20) +
    scale_fill_manual(values=c(color[5],color[6]))
dev.off()


##########################################################################
#################### HT gene index of GE-like phage SVs ##################
##########################################################################


scp Laisenying@10.190.248.213:/share/home1/Laisenying/Data-analysis/projects/PhageSV/PacBio_SV/CAST_merge/Population/Result/BacSource_SV/ViralGene_allHGTscore_inf.tsv /Users/laisenying/Desktop/ViralGene_allHGTscore_inf.tsv



HGT_inf = read.csv("/Users/laisenying/Desktop/ViralGene_allHGTscore_inf.tsv",sep="\t")

color <- c(brewer.pal(12,'Set3'))
HGT_inf$Type = factor(HGT_inf$Type,levels=c("Conserved region","noGE-like SVs","GE-like SVs"))

pdf("/Users/laisenying/Desktop/HGT_gene_index_GESV.pdf",width=4.5,height=6.5)
ggplot(HGT_inf[HGT_inf$Viral_genome == "NC173_NODE_1",],aes(x=Type,y=HT_index)) +
    geom_violin(aes(fill=Type),width=0.6,outlier.size=0) +
    geom_boxplot(aes(fill=Type),color='white',width=0.1,outlier.size=0) +
    #geom_jitter(aes(color=Type),width=0.2,size=1,alpha=1) +
    theme_bw() +
    theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size=14,color="black",angle=60,hjust=1),
    axis.text.y = element_text(size=14,color="black"),
    axis.title = element_text(size=16,color="black"),
    plot.title = element_text(hjust=0.5,size=18),
    panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
    legend.text = element_text(size=14,color="black"),
    legend.title = element_text(size=14,color="black"),
    axis.ticks = element_blank()
    )  +
    labs(x=" ",y="HT score",fill=" ") +
    scale_fill_manual(values=c(color[3],color[4],color[5],color[6])) +
    scale_color_manual(values=c(color[3],color[4],color[5],color[6])) +
    scale_x_discrete(limits=c("Conserved region","noGE-like SVs","GE-like SVs")) +
    guides(fill=FALSE) +
    guides(color=FALSE) +
    geom_signif(comparisons=list(c("Conserved region","noGE-like SVs")),y_position=-3.3) +
    geom_signif(comparisons=list(c("GE-like SVs","noGE-like SVs")),y_position=-3.2) +
    geom_signif(comparisons=list(c("GE-like SVs","Conserved region")),y_position=-3.1)
dev.off()


#######################################################################################
#################### Identify HGT genes between viruses and bacteria ##################
#######################################################################################

HGV_protein=/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/BacViral_protein/allViral_protein.pep
output=/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/BacViral_protein
/home1/Laisenying/Tools/cd-hit/cd-hit -i ${HGV_protein} -o ${output}/allViral_cdhit95_protein.pep  -d 200 -c 0.95 -uL 0.1 -uS 0.1 -n 5 -T 8 -M 16000 > ${output}/Vial_cd-hit95.log

HGV_protein=/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/BacViral_protein/allBac_protein.pep
output=/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/BacViral_protein
/home1/Laisenying/Tools/cd-hit/cd-hit -i ${HGV_protein} -o ${output}/allBac_cdhit95_protein.fa  -d 200 -c 0.95 -uL 0.1 -uS 0.1 -n 5 -T 8 -M 16000 > ${output}/cd-hit95.log


Viral_rep_protein=/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/BacViral_protein/allViral_cdhit_protein.pep
Bac_rep_protein=/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/BacViral_protein/allBac_protein.pep
cat ${Viral_rep_protein} ${Bac_rep_protein} > /home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/BacViral_protein/BacViral_protein.pep


# 多序列比对建树
for PC in `cat /share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/Protein_family/all_PC_new.list`
#for PC in `cat /home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/Protein_family/GESV_PCs.list`
do
PC_dir=/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/Protein_family/PC_fasta_new
if [ ! -s "/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/Protein_family/MSA_file_new/${PC}_align.fa" ];
then
/home1/Laisenying/Tools/MAFFT/mafft-7.158-without-extensions/core/bin/mafft \
    --auto \
    ${PC_dir}/${PC}.fasta > /home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/Protein_family/MSA_file_new/${PC}_align.fa
fi

if [ ! -s "/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/IQtree_output/${PC}_align_filter.fa" ];
then
/home1/Laisenying/Tools/trimal-trimAl/source/trimal \
    -gt 0.2 \
    -in /home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/Protein_family/MSA_file_new/${PC}_align.fa \
    -out /home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/IQtree_output/${PC}_align_filter.fa
fi

if [ ! -s "/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/IQtree_output/${PC}_align_filter.fa.treefile" ];
then
/home1/Laisenying/Tools/iqtree-1.6.12-Linux/bin/iqtree \
    -s /home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/IQtree_output/${PC}_align_filter.fa \
    -m LG+F+R5 \
    -alrt 1000
fi

if [ ! -s "/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/IQtree_output/${PC}_align_filter.fa.treefile.rooted" ];
then
/home1/Laisenying/Tools/mad/mad -m  /home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/IQtree_output/${PC}_align_filter.fa.treefile
sed -i '/^$/d' /share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/IQtree_output/${PC}_align_filter.fa.treefile.rooted
sed -i '/^>>/d' /share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/IQtree_output/${PC}_align_filter.fa.treefile.rooted
sed -i '/^>/d' /share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/IQtree_output/${PC}_align_filter.fa.treefile.rooted
fi
done

protein_inf = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/BacViral_protein/BacViral_all_protein_inf.tsv",sep="\t",index_col=0)
protein_inf.index = ['_'.join(x.split("~")) for x in protein_inf['geneID']]
PC_inf = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/Protein_family/BacViral_protein_family_inf_filter_new.tsv",sep="\t",index_col=0)
HGT_data = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/HGT_result/NVOG_hgt_summary.2022_09_27.tab",sep="\t")
HGT_data['Re_tax'] = [protein_inf.loc[x.split(",")[0],'Tax'] for x in HGT_data['recipients']]
HGT_data.loc[(HGT_data['direction']=='unkn')*(HGT_data['Re_tax']=='Viral'),'direction'] = 'Uncertain PtoV'
HGT_data.loc[(HGT_data['direction']=='unkn')*(HGT_data['Re_tax']=='Prokaryotic'),'direction'] = 'Uncertain VtoP'


protein_inf.index = ['_'.join(x.split("~")) for x in protein_inf['geneID']]
VtoP_HGT_genes = []
for x in HGT_data.loc[(HGT_data['direction']=='VtoP')+(HGT_data['direction']=='Uncertain VtoP'),'donors']:
    VtoP_HGT_genes.extend(x.split(','))


VtoP_HGT_genes = list(set(VtoP_HGT_genes))
VtoP_reps = set(protein_inf.loc[set(VtoP_HGT_genes).intersection(protein_inf.index),'Rep'])
protein_inf.index = protein_inf['Rep']
VtoP_HGT_genes_all = set(protein_inf.loc[VtoP_reps,'geneID']) # 7,9033 genes

protein_inf.index = ['_'.join(x.split("~")) for x in protein_inf['geneID']]
PtoV_HGT_genes = []
for x in HGT_data.loc[(HGT_data['direction']=='PtoV')+(HGT_data['direction']=='Uncertain PtoV'),'recipients']:
    PtoV_HGT_genes.extend(x.split(','))
    
PtoV_reps = set(protein_inf.loc[PtoV_HGT_genes,'Rep'])
protein_inf.index = protein_inf['Rep']
PtoV_HGT_genes_all = set(protein_inf.loc[PtoV_reps,'geneID']) # 33,572 genes

ViralSV_gene_inf = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/PacBio_SV/CAST_merge/Population/Result/SV_protein_function/SV_gene_function_all_category_noKO_inf.tsv",sep="\t",index_col=0)
SVV_genes = set(ViralSV_gene_inf.loc[ViralSV_gene_inf['SVID']==ViralSV_gene_inf['SVID'],'geneID']) # 22,667

HGT_overlaped_SV_genes = PtoV_HGT_genes_all.intersection(SVV_genes) # 3,506 P-to-V genes located in SV regions
HGT_nooverlaped_SV_genes = PtoV_HGT_genes_all.difference(SVV_genes) # 30,066

Carry_HGT_gene_SVs = set(ViralSV_gene_inf.loc[HGT_overlaped_SV,'SVID']) # 2,147 SVs / 6,273 SVs

GE_SVs = set(ViralSV_gene_inf.loc[ViralSV_gene_inf['SV_Source']=='bacteria','SVID']) # 2587 SVs

BV_interactive_data = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/Analysis/BV_coevolution/Interactive_range/BV_all_interactive_data.tsv",sep="\t",index_col=0)

GEV_genes = set(ViralSV_gene_inf.loc[ViralSV_gene_inf['SV_Source']=='bacteria','geneID']) # 7395

PtoV_HGT_genes_all.intersection(GEV_genes)



# all_viral_genes: 768788
# all P-to-V genes: 30,066
# all genes in SVs: 22,667
# P-to-V genes located within SVs: 3,506
# GEV_genes: 7395
# GE-like SVs: 2587 SVs
# Number of GE-like SVs carrying HGT genes: 957 / 2587
# Number of noGE-like SVs carrying HGT genes: 1190 / 5856

# P-to-V genes located within GElike SVs: 1625 / (7395 genes in GE-like SVs)
# P-to-V genes located within noGE-like SVs: 1881 / (15272 genes in noGE-like SVs)
# P-to-V genes located within conserved regions: 30066 / (746121 genes in conserved regions)


data_plot = data.frame(Region=c("Conserved region","noGE-like SVs","GE-like SVs"),Proportion=c(30066/747121,1881/15272,1625/7395))
color <- c(brewer.pal(12,'Set3'))
data_plot$Region = factor(data_plot$Region,levels=c("Conserved region","noGE-like SVs","GE-like SVs"))
color.category = c(
    "Conserved region" = color[3],
    "noGE-like SVs" = color[4],
    "GE-like SVs" = color[5]
)

pdf("/Users/laisenying/Desktop/HGT_proportion.pdf",width=4,height=6.5)
ggplot(data_plot,aes(x=Region,y=Proportion)) +
    geom_bar(aes(fill=Region),stat='identity',width=0.7,color='black') +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size=14,color="black",angle=60,hjust=1),
        axis.text.y = element_text(size=14,color="black"),
        axis.title = element_text(size=16,color="black"),
        plot.title = element_text(hjust=0.5,size=18),
        panel.border = element_rect(fill=NA,color="black",size=0.8,linetype="solid"),
        legend.text = element_text(size=14,color="black"),
        legend.title = element_text(size=14,color="black"),
        axis.ticks = element_blank()
    )  +
    labs(x=" ",y="Proportion of B-to-P genes",fill=" ") +
    scale_fill_manual(values=color.category,labels=names(color.category)) +
    guides(fill=FALSE)
dev.off()



###########################################################################################################
#################### Functional difference of GE-like SVs and noGE-like SVs ###############################
###########################################################################################################

data = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/Functional_enrichment_GEvsnoGE_SVs.tsv",sep="\t")

data_INSDEL = data[,c("Functional_category","INSDELCount","GeneCount","INSDEL_pvalue","INSDEL_foldchange","Group")]
colnames(data_INSDEL) = c("Functional_category","SV_GeneCount","GeneCount","pvalue","Foldchange","Group")
data_INSDEL$SVTYPE = "INS & DEL"
data_INSDEL$fdr = p.adjust(data_INSDEL$pvalue,method='fdr')
INSDEL_functions = unique(data_INSDEL[data_INSDEL$fdr<0.05,'Functional_category'])

data_INV = data[,c("Functional_category","INVCount","GeneCount","INV_pvalue","INV_foldchange","Group")]
colnames(data_INV) = c("Functional_category","SV_GeneCount","GeneCount","pvalue","Foldchange","Group")
data_INV$SVTYPE = "INV"
data_INV$fdr = p.adjust(data_INV$pvalue,method='fdr')
INV_functions = unique(data_INV[data_INV$fdr<0.05,'Functional_category'])

data_DUP = data[,c("Functional_category","DUPCount","GeneCount","DUP_pvalue","DUP_foldchange","Group")]
colnames(data_DUP) = c("Functional_category","SV_GeneCount","GeneCount","pvalue","Foldchange","Group")
data_DUP$SVTYPE = "DUP"
data_DUP$fdr = p.adjust(data_DUP$pvalue,method='fdr')
DUP_functions = unique(data_DUP[data_DUP$fdr<0.05,'Functional_category'])

data_allSV = data[,c("Functional_category","SVCount","GeneCount","SV_pvalue","SV_foldchange","Group")]
colnames(data_allSV) = c("Functional_category","SV_GeneCount","GeneCount","pvalue","Foldchange","Group")
data_allSV$SVTYPE = 'allSVs'
data_allSV$fdr = p.adjust(data_allSV$pvalue,method='fdr')
all_functions = unique(data_allSV[data_allSV$fdr<0.05,'Functional_category'])


Category = c('HUH','DDE','Tyr','Ser','Integration & Integrase','Conjugative system protein',
'DNA transformation protein','Antibiotic resistance','Virulence-associated protein','Topoisomerase','Relaxase','Resolvase',
'CRISPR-cas system','DNA packaging protein','Methyltransferase','Methylase')

data_plot = data_allSV[data_allSV$Functional_category %in% Category,]
subdata = data_plot[data_plot$Group=='GE-like SVs',]
order_functions = subdata$Functional_category[order(subdata$Foldchange)]

data_plot$Functional_category = factor(data_plot$Functional_category,levels=order_functions)
data_plot = data_plot[order(data_plot$Functional_category),]

data_plot$Group = factor(data_plot$Group,levels=c("noGE-like SVs","GE-like SVs"),labels=c("noGE-like","GE-like"))

p_all = ggplot(data_plot,aes(x=Foldchange,y=Functional_category)) +
    geom_bar(aes(fill=Group),stat='identity',position='dodge',width=0.9,color='black') +
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
    )  +
    scale_y_discrete(limits=order_functions) +
    scale_fill_manual(values=c(color[5],color[4])) +
    labs(x="Fold change",y=" ",fill='SV category',title=' ') +
    geom_vline(xintercept=1,linetype='dashed')


data_plot = data_INSDEL[data_INSDEL$Functional_category %in% Category,]
subdata = data_plot[data_plot$Group=='GE-like SVs',]
order_functions = subdata$Functional_category[order(subdata$Foldchange)]

data_plot$Group = factor(data_plot$Group,levels=c("noGE-like SVs","GE-like SVs"),labels=c("noGE-like","GE-like"))

p_INSDEL = ggplot(data_plot,aes(x=Foldchange,y=Functional_category)) +
    geom_bar(aes(fill=Group),stat='identity',position='dodge',width=0.9,color='black') +
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
    )  +
    scale_y_discrete(limits=order_functions) +
    scale_fill_manual(values=c(color[5],color[4])) +
    labs(x="Fold change",y=" ",fill='SV category',title='INS & DEL') +
    geom_vline(xintercept=1,linetype='dashed')

data_plot = data_INV[data_INV$Functional_category %in% Category,]
subdata = data_plot[data_plot$Group=='GE-like SVs',]
order_functions = subdata$Functional_category[order(subdata$Foldchange)]

data_plot$Group = factor(data_plot$Group,levels=c("noGE-like SVs","GE-like SVs"),labels=c("noGE-like","GE-like"))

p_INV = ggplot(data_plot,aes(x=Foldchange,y=Functional_category)) +
    geom_bar(aes(fill=Group),stat='identity',position='dodge',width=0.9,color='black') +
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
    )  +
    scale_y_discrete(limits=order_functions) +
    scale_fill_manual(values=c(color[5],color[4])) +
    labs(x="Fold change",y=" ",fill='SV category',title='INV') +
    geom_vline(xintercept=1,linetype='dashed')


data_plot = data_DUP[data_DUP$Functional_category %in% Category,]
subdata = data_plot[data_plot$Group=='GE-like SVs',]
order_functions = subdata$Functional_category[order(subdata$Foldchange)]

data_plot$Group = factor(data_plot$Group,levels=c("noGE-like SVs","GE-like SVs"),labels=c("noGE-like","GE-like"))

p_DUP = ggplot(data_plot,aes(x=Foldchange,y=Functional_category)) +
    geom_bar(aes(fill=Group),stat='identity',position='dodge',width=0.9,color='black') +
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
    )  +
    scale_y_discrete(limits=order_functions) +
    scale_fill_manual(values=c(color[5],color[4])) +
    labs(x="Fold change",y=" ",fill='SV category',title='DUP') +
    geom_vline(xintercept=1,linetype='dashed')

pdf("/Users/laisenying/Desktop/Functional_enrich_allSV.pdf",width=7.5,height=5.5)
p_all
dev.off()

pdf("/Users/laisenying/Desktop/Functional_enrich_INSDEL.pdf",width=7.5,height=5.5)
p_INSDEL
dev.off()

pdf("/Users/laisenying/Desktop/Functional_enrich_INV.pdf",width=7.5,height=5.5)
p_INV
dev.off()

pdf("/Users/laisenying/Desktop/Functional_enrich_DUP.pdf",width=7.5,height=5.5)
p_DUP
dev.off()

##########################################################################
#### Distribution of identity ####
##########################################################################

BV_interactive_data = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/BV_all_interactive_data.tsv",sep="\t")
BV_interactive_data_HGV = BV_interactive_data[BV_interactive_data$Viral_Source=='HGV',]
BV_interactive_data_HGV = BV_interactive_data_HGV[BV_interactive_data_HGV$Direction=='bacteria-to-viralSV',]
BV_interactive_data_HGV = BV_interactive_data_HGV[order(-BV_interactive_data_HGV$identity),]

BV_interactive_data_HGV_CH = BV_interactive_data_HGV[BV_interactive_data_HGV$Bacterial_Source=='CH-Binning',]
BV_interactive_data_HGV_filter = BV_interactive_data_HGV_CH[!duplicated(BV_interactive_data_HGV_CH$NR_SVID),]


BV_interactive_data_HGV_filter$Group='=100%'
BV_interactive_data_HGV_filter$Group[BV_interactive_data_HGV_filter$identity<100]='>99%'
BV_interactive_data_HGV_filter$Group[BV_interactive_data_HGV_filter$identity<99]='>95%'
BV_interactive_data_HGV_filter$Group[BV_interactive_data_HGV_filter$identity<95]='>90%'
BV_interactive_data_HGV_filter$Group[BV_interactive_data_HGV_filter$identity<90]='>85%'
BV_interactive_data_HGV_filter$Group[BV_interactive_data_HGV_filter$identity<85]='>80%'

plot_bar = data.frame(table(BV_interactive_data_HGV_filter$Group))
plot_bar$Group = factor(plot_bar$Var1,levels=c("=100%",">99%",">95%",">90%",">85%",">80%"))

color <- c(brewer.pal(9,'Set1'))
p1 = ggplot(plot_bar,aes(x=Group,y=Freq)) +
    geom_bar(aes(fill=Group),width=0.6,stat='identity',color='black') +
    geom_text(aes(label=Freq,y=Freq+30)) +
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



BV_interactive_data = read.csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/BV_all_interactive_data.tsv",sep="\t")
BV_interactive_data_HGV = BV_interactive_data[BV_interactive_data$Viral_Source=='HGV',]
BV_interactive_data_HGV = BV_interactive_data_HGV[BV_interactive_data_HGV$Direction=='bacteria-to-viralSV',]
BV_interactive_data_HGV = BV_interactive_data_HGV[order(-BV_interactive_data_HGV$identity),]

BV_interactive_data_HGV_HG = BV_interactive_data_HGV[BV_interactive_data_HGV$Bacterial_Source=='HumGut',]
BV_interactive_data_HGV_filter = BV_interactive_data_HGV_HG[!duplicated(BV_interactive_data_HGV_HG$NR_SVID),]


BV_interactive_data_HGV_filter$Group='=100%'
BV_interactive_data_HGV_filter$Group[BV_interactive_data_HGV_filter$identity<100]='>99%'
BV_interactive_data_HGV_filter$Group[BV_interactive_data_HGV_filter$identity<99]='>95%'
BV_interactive_data_HGV_filter$Group[BV_interactive_data_HGV_filter$identity<95]='>90%'
BV_interactive_data_HGV_filter$Group[BV_interactive_data_HGV_filter$identity<90]='>85%'
BV_interactive_data_HGV_filter$Group[BV_interactive_data_HGV_filter$identity<85]='>80%'

plot_bar = data.frame(table(BV_interactive_data_HGV_filter$Group))
plot_bar$Group = factor(plot_bar$Var1,levels=c("=100%",">99%",">95%",">90%",">85%",">80%"))

color <- c(brewer.pal(9,'Set1'))
p2 = ggplot(plot_bar,aes(x=Group,y=Freq)) +
    geom_bar(aes(fill=Group),width=0.6,stat='identity',color='black') +
    geom_text(aes(label=Freq,y=Freq+30)) +
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

library(cowplot)

pdf("/Users/laisenying/Desktop/Fig. S4.pdf",width=12,height=6)
plot_grid(p1,p2,nrow=1)
dev.off()


###########################################################################################################
################## Genetic transfer of important functions between phages and bacteria ####################
###########################################################################################################

protein_inf = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/BacViral_protein/BacViral_all_protein_inf.tsv",sep="\t",index_col=0)
protein_inf.index = ['_'.join(x.split("~")) for x in protein_inf['geneID']]
PC_inf = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/Protein_family/BacViral_protein_family_inf_filter_new.tsv",sep="\t",index_col=0)
HGT_data = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/HGT_result/NVOG_hgt_summary.2022_09_27.tab",sep="\t")
HGT_data_PtoV = HGT_data.loc[HGT_data['direction']=='PtoV',]
HGT_data_PtoV
ViralSV_gene_inf = pd.read_csv("/share/home1/Laisenying/Data-analysis/projects/PhageSV/PacBio_SV/CAST_merge/Population/Result/SV_protein_function/SV_gene_function_all_category_noKO_inf.tsv",sep="\t",index_col=0)

PC='PC99'
HGT_data_PtoV.loc[HGT_data_PtoV['COG']==PC,]

PC_geneIDs = list(PC_inf.loc[PC_inf['PC']==PC,'protein'])
protein_inf.index = protein_inf['Rep']
selected_genes = protein_inf.loc[PC_geneIDs,'geneID']

from Bio import SeqIO

Seq = SeqIO.parse("/share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/BacViral_protein/allBac_protein.pep","fasta")
Seq = (record for record in Seq if record.id in selected_genes)
SeqIO.write(Seq,"/home1/Laisenying/Projects/PhageSV/HGT_analysis/PC99/Bac_PC99.pep","fasta")

Seq = SeqIO.parse("/share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/BacViral_protein/allViral_protein.pep","fasta")
Seq = (record for record in Seq if record.id in selected_genes)
SeqIO.write(Seq,"/home1/Laisenying/Projects/PhageSV/HGT_analysis/PC99/Viral_PC99.pep","fasta")



scp Laisenying@10.190.248.213:/share/home1/Laisenying/Data-analysis/projects/PhageSV/HGT_identify/Phylogeny_based/IQtree_output/PC99.newtreefile PC99.newtreefile

Tree_PC = PC_inf.loc[PC_inf['PC']==PC,]
Tree_PC['protein'] = ["_".join(x.split("~")) for x in Tree_PC['protein']]
map_colors = {"Viral":"#80B1D3","Prokaryotic":"#FB8072"}
Tree_PC_Group = Tree_PC.loc[:,['protein','Tax']]
Tree_PC_Group['Color'] = [map_colors[x] for x in Tree_PC_Group['Tax']]
Tree_PC_Group['type'] = 'branch'
Tree_PC_Group['style'] = 'normal'
Tree_PC_Group['size'] = 3
Tree_PC_Group = Tree_PC_Group.loc[:,['protein','type','Color','style','size']]
Tree_PC_Group.columns = ['geneID','type','Color','style','size']
Tree_PC_Group.to_csv("/home1/Laisenying/Projects/PhageSV/HGT_analysis/PC99/Tree_annotation_Group.txt",sep="\t",header=None,index=0)
scp Laisenying@10.190.248.213:/home1/Laisenying/Projects/PhageSV/HGT_analysis/PC99/Tree_annotation_Group.txt Tree_annotation_Group.txt
