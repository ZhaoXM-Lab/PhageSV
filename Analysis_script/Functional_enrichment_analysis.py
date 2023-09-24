import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.stats import fisher_exact,hypergeom
SV_gene_inf = pd.read_csv("Results/ViralSV_gene_functional_category.tsv",sep="\t",index_col=0)
SV_gene_inf = SV_gene_inf.drop_duplicates(['geneID','SVID','SVTYPE','Level2_viral_category'])
SV_gene_inf.loc[SV_gene_inf['Level2_viral_category']!=SV_gene_inf['Level2_viral_category'],'Level2_viral_category'] = 'Others'
SV_KEGG_data = SV_gene_inf
KOs = set(SV_KEGG_data['Level2_viral_category'])
SV_KEGG_data.index = SV_KEGG_data['Level2_viral_category']
All_gene_count = len(set(SV_KEGG_data['geneID']))
INSDEL_gene_count = len(set(SV_KEGG_data.loc[(SV_KEGG_data['SVTYPE']=='INS')+(SV_KEGG_data['SVTYPE']=='DEL'),'geneID']))
DUP_gene_count = len(set(SV_KEGG_data.loc[SV_KEGG_data['SVTYPE']=='DUP','geneID']))
INV_gene_count = len(set(SV_KEGG_data.loc[SV_KEGG_data['SVTYPE']=='INV','geneID']))
KO_enrich = defaultdict(list)
for ko in KOs:
    print(ko)
    subdata = pd.DataFrame(SV_KEGG_data.loc[ko,])
    if subdata.shape[1]==1:
        subdata = subdata.T
    All_count = len(set(subdata['geneID']))
    INSDEL_count = len(set(subdata.loc[(subdata['SVTYPE']=='INS')+(subdata['SVTYPE']=='DEL'),'geneID']))
    DUP_count = len(set(subdata.loc[subdata['SVTYPE']=='DUP','geneID']))
    INV_count = len(set(subdata.loc[subdata['SVTYPE']=='INV','geneID']))
    subunique = subdata.drop_duplicates(['geneID'])
    stat=pd.value_counts(subunique['Category'])
    if 'viral' in stat.index:
        viral_proportion = stat['viral']/All_count
    else:
        viral_proportion = 0
    if 'microbial' in stat.index:
        microbial_proportion=stat['microbial']/All_count
    else:
        microbial_proportion=0
    KO_enrich['Level2_category'].append(ko)
    KO_enrich['Level1_category'].append(pd.value_counts(subdata['Viral_category']).index[0])
    KO_enrich['microbial_proportion'].append(microbial_proportion)
    KO_enrich['viral_proportion'].append(viral_proportion)
    KO_enrich['GeneCount'].append(All_count)
    KO_enrich['INSDELCount'].append(INSDEL_count)
    KO_enrich['DUPCount'].append(DUP_count)
    KO_enrich['INVCount'].append(INV_count)
    KO_enrich['BgRatio'].append(All_count/All_gene_count)
    KO_enrich['INSDELRatio'].append(INSDEL_count/INSDEL_gene_count)
    KO_enrich['INVRatio'].append(INV_count/INV_gene_count)
    KO_enrich['DUPRatio'].append(DUP_count/DUP_gene_count)
    INSDELf = fisher_exact([[INSDEL_count,INSDEL_gene_count-INSDEL_count],[All_count,All_gene_count-All_count]],alternative='greater')
    INVf = fisher_exact([[INV_count,INV_gene_count-INV_count],[All_count,All_gene_count-All_count]],alternative='greater')
    DUPf = fisher_exact([[DUP_count,DUP_gene_count-DUP_count],[All_count,All_gene_count-All_count]],alternative='greater')
    KO_enrich['INSDEL_pvalue'].append(INSDELf[1])
    KO_enrich['INV_pvalue'].append(INVf[1])
    KO_enrich['DUP_pvalue'].append(DUPf[1])

SV_pathway_enrich_data = pd.DataFrame(KO_enrich)
SV_pathway_enrich_data['INSDEL_foldchange'] = SV_pathway_enrich_data['INSDELRatio']/SV_pathway_enrich_data['BgRatio']
SV_pathway_enrich_data['DUP_foldchange'] = SV_pathway_enrich_data['DUPRatio']/SV_pathway_enrich_data['BgRatio']
SV_pathway_enrich_data['INV_foldchange'] = SV_pathway_enrich_data['INVRatio']/SV_pathway_enrich_data['BgRatio']
SV_pathway_enrich_data = SV_pathway_enrich_data.loc[SV_pathway_enrich_data['Level1_category']!='Unknown_function',]
SV_pathway_enrich_data = SV_pathway_enrich_data.loc[SV_pathway_enrich_data['Level2_category']!='Others',]
SV_pathway_enrich_data = SV_pathway_enrich_data.loc[SV_pathway_enrich_data['Level1_category']!='Hypothetical protein',]

SV_pathway_enrich_data.to_csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/ViralSV_Manual_Category_greater_enriched.tsv",sep="\t")

###########################################################################################################
#################### Functional difference of GE-like SVs and noGE-like SVs ###############################
###########################################################################################################
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.stats import fisher_exact,hypergeom

ViralSV_gene_inf = pd.read_csv("Results/ViralSV_gene_functional_category.tsv",sep="\t",index_col=0)
ViralSV_gene_inf.loc[ViralSV_gene_inf['GO_name']=='GO: CRISPR-cas system','Level2_viral_category'] = 'CRISPR-cas system'
ViralSV_gene_inf.loc[['CRISPR' in str(x) for x in ViralSV_gene_inf['Pfam_function']],'Level2_viral_category'] = 'CRISPR-cas system'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Functional_category']=='Antibiotic resistance','Level2_viral_category']
ViralSV_gene_inf.loc[ViralSV_gene_inf['GO_name']=='GO: response to antibiotic','Level2_viral_category']  = 'Antibiotic resistance'
ViralSV_gene_inf.loc[ViralSV_gene_inf['GO_name']=='GO:response to antibiotic','Level2_viral_category'] = 'Antibiotic resistance'
ViralSV_gene_inf.loc[ViralSV_gene_inf['GO_name']=='GO: transposition','Level2_viral_category'] = 'Transposase'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Level2_viral_category']=='DNA integration','Level2_viral_category'] = 'Recombination & Recombinase'
ViralSV_gene_inf.loc[(ViralSV_gene_inf['GO_name']=='GO: DNA recombination')*(ViralSV_gene_inf['Level2_viral_category']!=ViralSV_gene_inf['Level2_viral_category']),'Level2_viral_category'] = 'Recombination & Recombinase'
ViralSV_gene_inf.loc[(ViralSV_gene_inf['GO_name']=='GO:DNA methylation')*(ViralSV_gene_inf['Level2_viral_category']!=ViralSV_gene_inf['Level2_viral_category']),'Level2_viral_category'] = 'Methylase'
ViralSV_gene_inf.loc[(ViralSV_gene_inf['GO_name']=='GO: rRNA base methylation')*(ViralSV_gene_inf['Level2_viral_category']=='Unknown_function'),'Level2_viral_category'] = 'rRNA base methylation'
ViralSV_gene_inf.loc[(ViralSV_gene_inf['GO_name']=='GO: rRNA base methylation')*(ViralSV_gene_inf['Level2_viral_category']!=ViralSV_gene_inf['Level2_viral_category']),'Level2_viral_category'] = 'rRNA base methylation'
ViralSV_gene_inf.loc[(ViralSV_gene_inf['GO_name']=='GO:transposition, DNA-mediated'),'Level2_viral_category'] = 'Transposition & Transposase'
ViralSV_gene_inf.loc[(ViralSV_gene_inf['GO_name']=='GO:transposition, DNA-mediated')*(ViralSV_gene_inf['Level2_viral_category']!=ViralSV_gene_inf['Level2_viral_category']),'Level2_viral_category'] = 'Transposition & Transposase'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Level2_viral_category']=='Transposase','Level2_viral_category'] = 'Transposition & Transposase'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Level2_viral_category']=='CRISPR-cas system','Viral_category']= 'CRISPR-cas system'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Level2_viral_category']=='Antibiotic resistance','Viral_category']= 'Immune evasion'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Level2_viral_category']=='rRNA base methylation','Viral_category']= 'Immune evasion'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Viral_category']=='Transduction','Viral_category'] = 'Regulation'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Level2_viral_category']=='CRISPR-cas system','Functional_category'] = 'CRISPR-cas system'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Level2_viral_category']=='Transposition & Transposase','Functional_category'] = 'Transposition & Transposase'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Functional_category']=='Transposase','Functional_category'] = 'Transposition & Transposase'
ViralSV_gene_inf.loc[(ViralSV_gene_inf['Level2_viral_category']=='Recombination & Recombinase')*(ViralSV_gene_inf['Functional_category']!=ViralSV_gene_inf['Functional_category']),'Functional_category'] = 'Recombination & Recombinase'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Functional_category']=='Integrase','Functional_category'] = 'Integration & Integrase'
ViralSV_gene_inf.loc[ViralSV_gene_inf['Functional_category']=='DNA mediated transformation','Functional_category'] = 'DNA transformation protein'
ViralSV_inf = pd.read_csv("/Users/laisenying/Desktop/Temp_file/Projects/PhageomeSV/Script v3.0/Results/SV_bed_inf.tsv",sep="\t",index_col=0)
Viral_SVLEN_inf = ViralSV_inf.loc[:,['SVID','SVLEN']]
ViralSV_gene_inf = pd.merge(ViralSV_gene_inf,Viral_SVLEN_inf,on='SVID',how='left')
ViralSV_gene_inf['SVLEN'] = np.abs(ViralSV_gene_inf['SVLEN'])
ViralSV_gene_inf['SVLEN'] = ViralSV_gene_inf['SVLEN'].fillna(0)
ViralGene_recombinase_inf = pd.read_csv("Results/proMGE_Viralprotein_recombinase.tsv",sep="\t",index_col=0)

HUH_genes = set(ViralGene_recombinase_inf.loc[ViralGene_recombinase_inf['Recombinase_category']=='HUH recombinase','geneID'])
Tyr_genes = set(ViralGene_recombinase_inf.loc[ViralGene_recombinase_inf['Recombinase_category']=='Tyr recombinase','geneID'])
Ser_genes = set(ViralGene_recombinase_inf.loc[ViralGene_recombinase_inf['Recombinase_category']=='Ser recombinase','geneID'])
Cas_genes = set(ViralGene_recombinase_inf.loc[ViralGene_recombinase_inf['Recombinase_category']=='Cas recombinase','geneID'])
DDE_genes = set(ViralGene_recombinase_inf.loc[ViralGene_recombinase_inf['Recombinase_category']=='DDE recombinase','geneID'])

ViralSV_gene_inf.loc[[x in HUH_genes for x in ViralSV_gene_inf['geneID']],'Hit_HUH'] = 'Y'
ViralSV_gene_inf.loc[[x in Tyr_genes for x in ViralSV_gene_inf['geneID']],'Hit_Tyr'] = 'Y'
ViralSV_gene_inf.loc[[x in Ser_genes for x in ViralSV_gene_inf['geneID']],'Hit_Ser'] = 'Y'
ViralSV_gene_inf.loc[[x in Cas_genes for x in ViralSV_gene_inf['geneID']],'Hit_Cas'] = 'Y'
ViralSV_gene_inf.loc[[x in DDE_genes for x in ViralSV_gene_inf['geneID']],'Hit_DDE'] = 'Y'

ViralSV_gene_inf.loc[ViralSV_gene_inf['Functional_category']!=ViralSV_gene_inf['Functional_category'],'Functional_category'] = ViralSV_gene_inf.loc[ViralSV_gene_inf['Functional_category']!=ViralSV_gene_inf['Functional_category'],'Level2_viral_category']

Viral_SV_gene_inf = ViralSV_gene_inf
Allcount=len(set(Viral_SV_gene_inf.loc[[len(x.split("~"))==1 for x in Viral_SV_gene_inf['Rep']],'geneID']))
INSDELcount = len(set(Viral_SV_gene_inf.loc[((Viral_SV_gene_inf['SVTYPE']=='INS')+(Viral_SV_gene_inf['SVTYPE']=='DEL'))*(Viral_SV_gene_inf['SV_Source']=='bacteria')*(Viral_SV_gene_inf['SVLEN']>=0),'geneID']))
DUPcount = len(set(Viral_SV_gene_inf.loc[((Viral_SV_gene_inf['SVTYPE']=='DUP')*(Viral_SV_gene_inf['SV_Source']=='bacteria'))*(Viral_SV_gene_inf['SVLEN']>=0),'geneID']))
INVcount = len(set(Viral_SV_gene_inf.loc[(Viral_SV_gene_inf['SVTYPE']=='INV')*(Viral_SV_gene_inf['SV_Source']=='bacteria')*(Viral_SV_gene_inf['SVLEN']>=0),'geneID']))
SVcount = len(set(Viral_SV_gene_inf.loc[(Viral_SV_gene_inf['SV_Source']=='bacteria')*(Viral_SV_gene_inf['SVLEN']>=0),'geneID']))
KO_enrich = defaultdict(list)

Category=set(Viral_SV_gene_inf['Level2_viral_category']).union(['HUH','DDE','Tyr','Ser','Cas','Integration & Integrase','Recombination & Recombinase','Transposition & Transposase','Conjugative system protein',
    'DNA transformation protein','Antibiotic resistance',
    'Exonuclease','Topoisomerase','Excisionase','Relaxase','Resolvase','CRISPR-cas system'])
    
    
for ko in Category:
    print(ko)
    if ko =='HUH':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_HUH']=='Y',]
    elif ko =='DDE':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_DDE']=='Y',]
    elif ko =='Ser':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Ser']=='Y',]
    elif ko =='Tyr':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Tyr']=='Y',]
    elif ko =='Cas':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Cas']=='Y',]
    elif ko == 'Transposition & Transposase':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_ISE']=='Y',]
    elif ko == 'Conjugative system protein':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['ConjScan']==Viral_SV_gene_inf['ConjScan'],]
    elif ko == 'Integration & Integrase':
        Viral_SV_gene_inf.loc[['Integrase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['integrase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['Integrase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['integrase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['Integrase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['integrase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['Integrase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['integrase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['Integrase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['integrase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_integrase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_integrase']=='Y',]
    elif ko == 'Recombination & Recombinase':
        Viral_SV_gene_inf.loc[['Recombinase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombinase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombinase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombinase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombinase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombinase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombinase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombinase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombinase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombinase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombination' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombination' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombination' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombination' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombination' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombination' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombination' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombination' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombination' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombination' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Recombinase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Recombinase']=='Y',]
    elif ko == 'Topoisomerase':
        Viral_SV_gene_inf.loc[['Topoisomerase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['topoisomerase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['Topoisomerase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['topoisomerase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['Topoisomerase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['topoisomerase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['Topoisomerase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['topoisomerase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['Topoisomerase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['topoisomerase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Topoisomerase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Topoisomerase']=='Y',]
    elif ko == 'Exonuclease':
        Viral_SV_gene_inf.loc[['Exonuclease' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['exonuclease' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['Exonuclease' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['exonuclease' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['Exonuclease' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['exonuclease' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['Exonuclease' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['exonuclease' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['Exonuclease' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['exonuclease' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Exonuclease'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Exonuclease']=='Y',]
    elif ko == 'Excisionase':
        Viral_SV_gene_inf.loc[['Excisionase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['excisionase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['Excisionase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['excisionase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['Excisionase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['excisionase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['Excisionase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['excisionase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['Excisionase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['excisionase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Excisionase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Exonuclease']=='Y',]
    elif ko == 'Relaxase':
        Viral_SV_gene_inf.loc[['Relaxase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['relaxase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['Relaxase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['relaxase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['Relaxase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['relaxase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['Relaxase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['relaxase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['Relaxase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['relaxase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Relaxase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Relaxase']=='Y',]
    elif ko == 'Resolvase':
        Viral_SV_gene_inf.loc[['Resolvase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['resolvase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['Resolvase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['resolvase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'HitResolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['Resolvase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['resolvase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['Resolvase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['resolvase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['Resolvase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['resolvase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Resolvase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Resolvase']=='Y',]
    else:
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Functional_category']==ko,]
    All_count1 = len(set(subdata['geneID']))
    All_count2 = len(set(subdata.loc[[len(x.split("~"))==1 for x in subdata['Rep']],'geneID'])) # Background
    INSDEL_count = len(set(subdata.loc[((subdata['SVTYPE']=='INS')+(subdata['SVTYPE']=='DEL'))*(subdata['SV_Source']=='bacteria')*(subdata['SVLEN']>=0),'geneID']))
    DUP_count = len(set(subdata.loc[(subdata['SVTYPE']=='DUP')*(subdata['SV_Source']=='bacteria')*(subdata['SVLEN']>=0),'geneID']))
    INV_count = len(set(subdata.loc[(subdata['SVTYPE']=='INV')*(subdata['SV_Source']=='bacteria')*(subdata['SVLEN']>=0),'geneID']))
    SV_count =  len(set(subdata.loc[(subdata['SV_Source']=='bacteria')*(subdata['SVLEN']>=0)*(subdata['SVTYPE']==subdata['SVTYPE']),'geneID']))
    KO_enrich['Functional_category'].append(ko)
    KO_enrich['INSDELCount'].append(INSDEL_count)
    KO_enrich['DUPCount'].append(DUP_count)
    KO_enrich['INVCount'].append(INV_count)
    KO_enrich['SVCount'].append(SV_count)
    KO_enrich['GeneCount'].append(All_count1)
    KO_enrich['BgRatio'].append(All_count2/Allcount)
    KO_enrich['INSDELRatio'].append(INSDEL_count/INSDELcount)
    KO_enrich['INVRatio'].append(INV_count/INVcount)
    KO_enrich['DUPRatio'].append(DUP_count/DUPcount)
    KO_enrich['SVRatio'].append(SV_count/SVcount)
    INSDELf = fisher_exact([[INSDEL_count,INSDELcount-INSDEL_count],[All_count2,Allcount-All_count2]],alternative='greater')
    INVf = fisher_exact([[INV_count,INVcount-INV_count],[All_count2,Allcount-All_count2]],alternative='greater')
    DUPf = fisher_exact([[DUP_count,DUPcount-DUP_count],[All_count2,Allcount-All_count2]],alternative='greater')
    SVf = fisher_exact([[SV_count,SVcount-SV_count],[All_count2,Allcount-All_count2]],alternative='greater')
    KO_enrich['INSDEL_pvalue'].append(INSDELf[1])
    KO_enrich['INV_pvalue'].append(INVf[1])
    KO_enrich['DUP_pvalue'].append(DUPf[1])
    KO_enrich['SV_pvalue'].append(SVf[1])

SV_enrich_data = pd.DataFrame(KO_enrich)
SV_enrich_data['INSDEL_foldchange'] = SV_enrich_data['INSDELRatio']/SV_enrich_data['BgRatio']
SV_enrich_data['DUP_foldchange'] = SV_enrich_data['DUPRatio']/SV_enrich_data['BgRatio']
SV_enrich_data['INV_foldchange'] = SV_enrich_data['INVRatio']/SV_enrich_data['BgRatio']
SV_enrich_data['SV_foldchange'] = SV_enrich_data['SVRatio']/SV_enrich_data['BgRatio']
SV_enrich_data['Group'] = 'GE-like SVs'
SV_enrich_data = SV_enrich_data.loc[SV_enrich_data['SV_foldchange']==SV_enrich_data['SV_foldchange'],]
SV_enrich_data_GE = SV_enrich_data


Allcount=len(set(Viral_SV_gene_inf.loc[[len(x.split("~"))==1 for x in Viral_SV_gene_inf['Rep']],'geneID']))
INSDELcount = len(set(Viral_SV_gene_inf.loc[((Viral_SV_gene_inf['SVTYPE']=='INS')+(Viral_SV_gene_inf['SVTYPE']=='DEL'))*(Viral_SV_gene_inf['SV_Source']!='bacteria')*(Viral_SV_gene_inf['SVLEN']>=0),'geneID']))
DUPcount = len(set(Viral_SV_gene_inf.loc[(Viral_SV_gene_inf['SVTYPE']=='DUP')*(Viral_SV_gene_inf['SV_Source']!='bacteria')*(Viral_SV_gene_inf['SVLEN']>=0),'geneID']))
INVcount = len(set(Viral_SV_gene_inf.loc[(Viral_SV_gene_inf['SVTYPE']=='INV')*(Viral_SV_gene_inf['SV_Source']!='bacteria')*(Viral_SV_gene_inf['SVLEN']>=0),'geneID']))
SVcount = len(set(Viral_SV_gene_inf.loc[(Viral_SV_gene_inf['SV_Source']!='bacteria')*(Viral_SV_gene_inf['SVLEN']>=0),'geneID']))
KO_enrich = defaultdict(list)


for ko in Category:
    print(ko)
    if ko =='HUH':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_HUH']=='Y',]
    elif ko =='DDE':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_DDE']=='Y',]
    elif ko =='Ser':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Ser']=='Y',]
    elif ko =='Tyr':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Tyr']=='Y',]
    elif ko =='Cas':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Cas']=='Y',]
    elif ko == 'Transposition & Transposase':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_ISE']=='Y',]
    elif ko == 'Conjugative system protein':
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['ConjScan']==Viral_SV_gene_inf['ConjScan'],]
    elif ko == 'Transpotision & Transposase':
        Viral_SV_gene_inf.loc[['Integrase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['integrase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['Integrase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['integrase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['Integrase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['integrase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['Integrase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['integrase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['Integrase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_integrase'] = 'Y'
        Viral_SV_gene_inf.loc[['integrase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_integrase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_integrase']=='Y',]
    elif ko == 'Recombination & Recombinase':
        Viral_SV_gene_inf.loc[['Recombinase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombinase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombinase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombinase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombinase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombinase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombinase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombinase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombinase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombinase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombination' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombination' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombination' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombination' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombination' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombination' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombination' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombination' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['Recombination' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Recombinase'] = 'Y'
        Viral_SV_gene_inf.loc[['recombination' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Recombinase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Recombinase']=='Y',]
    elif ko == 'Topoisomerase':
        Viral_SV_gene_inf.loc[['Topoisomerase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['topoisomerase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['Topoisomerase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['topoisomerase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['Topoisomerase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['topoisomerase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['Topoisomerase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['topoisomerase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['Topoisomerase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Topoisomerase'] = 'Y'
        Viral_SV_gene_inf.loc[['topoisomerase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Topoisomerase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Topoisomerase']=='Y',]
    elif ko == 'Exonuclease':
        Viral_SV_gene_inf.loc[['Exonuclease' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['exonuclease' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['Exonuclease' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['exonuclease' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['Exonuclease' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['exonuclease' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['Exonuclease' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['exonuclease' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['Exonuclease' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Exonuclease'] = 'Y'
        Viral_SV_gene_inf.loc[['exonuclease' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Exonuclease'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Exonuclease']=='Y',]
    elif ko == 'Excisionase':
        Viral_SV_gene_inf.loc[['Excisionase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['excisionase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['Excisionase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['excisionase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['Excisionase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['excisionase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['Excisionase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['excisionase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['Excisionase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Excisionase'] = 'Y'
        Viral_SV_gene_inf.loc[['excisionase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Excisionase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Exonuclease']=='Y',]
    elif ko == 'Relaxase':
        Viral_SV_gene_inf.loc[['Relaxase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['relaxase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['Relaxase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['relaxase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['Relaxase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['relaxase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['Relaxase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['relaxase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['Relaxase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Relaxase'] = 'Y'
        Viral_SV_gene_inf.loc[['relaxase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Relaxase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Relaxase']=='Y',]
    elif ko == 'Resolvase':
        Viral_SV_gene_inf.loc[['Resolvase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['resolvase' in str(x) for x in Viral_SV_gene_inf['VOGdb_function']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['Resolvase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['resolvase' in str(x) for x in Viral_SV_gene_inf['CheckV_function']],'HitResolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['Resolvase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['resolvase' in str(x) for x in Viral_SV_gene_inf['Description']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['Resolvase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['resolvase' in str(x) for x in Viral_SV_gene_inf['PFAMs']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['Resolvase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Resolvase'] = 'Y'
        Viral_SV_gene_inf.loc[['resolvase' in str(x) for x in Viral_SV_gene_inf['Pfam_function']],'Hit_Resolvase'] = 'Y'
        subdata =  Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Hit_Resolvase']=='Y',]
    else:
        subdata = Viral_SV_gene_inf.loc[Viral_SV_gene_inf['Functional_category']==ko,]
    All_count1 = len(set(subdata['geneID']))
    All_count2 = len(set(subdata.loc[[len(x.split("~"))==1 for x in subdata['Rep']],'geneID'])) # Background
    INSDEL_count = len(set(subdata.loc[((subdata['SVTYPE']=='INS')+(subdata['SVTYPE']=='DEL'))*(subdata['SV_Source']!='bacteria')*(subdata['SVLEN']>=0),'geneID']))
    DUP_count = len(set(subdata.loc[(subdata['SVTYPE']=='DUP')*(subdata['SV_Source']!='bacteria')*(subdata['SVLEN']>=0),'geneID']))
    INV_count = len(set(subdata.loc[(subdata['SVTYPE']=='INV')*(subdata['SV_Source']!='bacteria')*(subdata['SVLEN']>=0),'geneID']))
    SV_count =  len(set(subdata.loc[(subdata['SV_Source']!='bacteria')*(subdata['SVLEN']>=0)*(subdata['SVTYPE']==subdata['SVTYPE']),'geneID']))
    KO_enrich['Functional_category'].append(ko)
    KO_enrich['INSDELCount'].append(INSDEL_count)
    KO_enrich['DUPCount'].append(DUP_count)
    KO_enrich['INVCount'].append(INV_count)
    KO_enrich['SVCount'].append(SV_count)
    KO_enrich['GeneCount'].append(All_count1)
    KO_enrich['BgRatio'].append(All_count2/Allcount)
    KO_enrich['INSDELRatio'].append(INSDEL_count/INSDELcount)
    KO_enrich['INVRatio'].append(INV_count/INVcount)
    KO_enrich['DUPRatio'].append(DUP_count/DUPcount)
    KO_enrich['SVRatio'].append(SV_count/SVcount)
    INSDELf = fisher_exact([[INSDEL_count,INSDELcount-INSDEL_count],[All_count2,Allcount-All_count2]],alternative='greater')
    INVf = fisher_exact([[INV_count,INVcount-INV_count],[All_count2,Allcount-All_count2]],alternative='greater')
    DUPf = fisher_exact([[DUP_count,DUPcount-DUP_count],[All_count2,Allcount-All_count2]],alternative='greater')
    SVf = fisher_exact([[SV_count,SVcount-SV_count],[All_count2,Allcount-All_count2]],alternative='greater')
    KO_enrich['INSDEL_pvalue'].append(INSDELf[1])
    KO_enrich['INV_pvalue'].append(INVf[1])
    KO_enrich['DUP_pvalue'].append(DUPf[1])
    KO_enrich['SV_pvalue'].append(SVf[1])

SV_enrich_data = pd.DataFrame(KO_enrich)
SV_enrich_data['INSDEL_foldchange'] = SV_enrich_data['INSDELRatio']/SV_enrich_data['BgRatio']
SV_enrich_data['DUP_foldchange'] = SV_enrich_data['DUPRatio']/SV_enrich_data['BgRatio']
SV_enrich_data['INV_foldchange'] = SV_enrich_data['INVRatio']/SV_enrich_data['BgRatio']
SV_enrich_data['SV_foldchange'] = SV_enrich_data['SVRatio']/SV_enrich_data['BgRatio']
SV_enrich_data['Group'] = 'noGE-like SVs'
SV_enrich_data = SV_enrich_data.loc[SV_enrich_data['SV_foldchange']==SV_enrich_data['SV_foldchange'],]
SV_enrich_data_noGE = SV_enrich_data

SV_enrich_data_all = pd.concat([SV_enrich_data_noGE,SV_enrich_data_GE])
SV_enrich_data_all.to_csv("Results/Functional_enrichment_GEvsnoGE_SVs.tsv",sep="\t")


