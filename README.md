## Long-read sequencing reveals extensive gut phageome structural variation driven by genetic exchange with bacterial hosts

The human gut virome, consisting primarily of bacteriophages (phages), exhibits remarkable inter-personal diversity and significantly impacts intestinal ecology. While taxonomic variation is recognized, the fine-scale genetic variations in the gut virome, specially structural variations (SVs), remain underexplored. Exploring such variations is instrumental for unraveling phage evolution and deciphering their functional implications. Here, by employing virome-enriched long-read metagenomics sequencing across 91 individuals, we identified a total of 14,438 non-redundant viral SVs, and analysed their prevalence within the human gut virome. We noticed that a significant proportion of viral SVs arise from genetic exchange between phages and bacteria, with most viral SVs homologous to bacterial fragments and enriched for horizontal gene transfer (HGT) mechanisms. Further investigations showed that these SV sequences were transmitted between specific phage-bacteria pairs, particularly between phages and their respective bacterial hosts. Collectively, our findings provides novel insights into the genetic landscape of the human gut virome.

## Content 
- **PBSV_pipeline.sh**: Scripts for detecting phage structural variations across the CHGV cohort
- **Benchmarking_pipeline.sh**: Generating sumulated virome-enriched metagenomics datasets with known introduced SVs and benchmarking different long-read based SV calling methods with PSDP pipeline. 
- **Analysis_script/PhageSV_quality_Part1.R**: Overview of 14,438 non-redundant structural variants (SVs) in the human gut phageome.
- **Analysis_script/Functional_enrichment_analysis.py**: Functional enrichment analysis of phage SVs.
- **Analysis_script/PhageSV_function.R**: Visualization of enriched functions in phage SVs.
- **Analysis_script/PhageSV_GE.R**: Bacteria-to-phage HGT analysis
- **Analysis_script/PhageSV_BVinteraction.R**: Bactera-phage SV sharing network
- **Analysis_script/PhageSV_Phage_host_GE_analysis.R**: Phage-host genetical interaction analysis


## Data
- **Results/Sample_SV_common_0.8_suppl.vcf**: 14,334 non-redundant VCF file
- **Results/SV_bed_inf.tsv**: Viral SV information extracted from non-redundant VCF file
- **Results/all_genome_SNP_inf.tsv**: SNV for each viral contig extracted from the output file of inStrain
- **Results/all_gene_SNP_inf.tsv**: SNV for each viral gene extracted from the output file of inStrain
- **ViralSV_density_inf**: SV density per viral_contig per sample, taxonomy, Temperate lifestyle

## Contact
- Senying Lai: 19110850024@fudan.edu.cn
