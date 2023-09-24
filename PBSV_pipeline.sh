
export PATH=~/Tools/bcftools-1.13:$PATH
basepath=/home1/Laisenying/Data-analysis/projects/PhageSV

for sample in `cat sample.txt`
do

SubreadBam=${sample}.subreads.bam
outputdir=./${sample}
mkdir -p ${outputdir}

if [ ! -s "${SubreadBam%.bam}.fastq" ];
then
   bedtools bamtofastq \
   -i ${SubreadBam} \
    -fq ${SubreadBam%.bam}.fastq
fi

PacBio_read=${SubreadBam%.bam}.fastq
mash sketch -m 2 ${PacBio_read}
mkdir -p ${outputdir}/Mash_filter
mash screen CHGV_HQ_genome_path.txt.msh ${PacBio_read} > ${outputdir}/Mash_filter/screen.tab

python Scripts/Mash_select.py \
    --PacBio_file ${outputdir}/Mash_filter/screen.tab \
    --t 0.90 \
    --outputfile ${outputdir}/Mash_filter/Mash_filtered_ref.fasta

ViralGenome=${outputdir}/Mash_filter/Mash_filtered_ref.fasta
bamfile=${outputdir}/${sample}.mapping_Mash_ref.bam

pbmm2 align ${ViralGenome} ${SubreadBam} ${outputBam} --sort --preset CCS --sample ${sample}
samtools index ${outputBam}

# 1. pbsv
pbsv discover ${bamfile} ${outputdir}/${sample}.mapping_Mash_ref.svsig.gz
pbsv call ${ViralGenome} ${outputdir}/${sample}.mapping_Mash_ref.svsig.gz \
        ${outputdir}/${sample}.pbsv.vcf

# 2. sniffles
sniffles -i ${bamfile} -v ${outputdir}/${sample}.sniffle.vcf

# 3. cuteSV
cuteSV ${bamfile} ${ViralGenome} ${outputdir}/${sample}.cuteSV.vcf ${outputdir} \
    --max_cluster_bias_INS 1000 \
    --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5

# 4. SVIM
svim alignment ${outputdir} ${bamfile} ${ViralGenome}
mv ${outputdir}/variants.vcf ${outputdir}/${sample}.svim.vcf


# PSDP
sniffles_vcf=${outputdir}/${sample}.sniffle.vcf
cuteSV_vcf=${outputdir}/${sample}.cuteSV.vcf
pbsv_vcf=${outputdir}/${sample}.pbsv.vcf
svim_vcf=${outputdir}/${sample}.svim.vcf

# re-format
python3 Scripts/Sniffles_change_format.py \
    --raw ${sniffles_vcf} \
    --new ${outputdir}/${sniffles_vcf%.vcf}.format.vcf

python3 Scripts/cuteSV_change_format.py \
    --raw ${cuteSV_vcf} \
    --new ${outputdir}/${cuteSV_vcf%.vcf}.format.vcf

python3 Scripts/pbsv_change_format.py \
    --raw ${pbsv_vcf} \
    --new ${outputdir}/${pbsv_vcf%.vcf}.format.vcf

python3 Scripts/svim_change_format.py \
    --raw ${svim_vcf} \
    --new ${outputdir}/${svim_vcf%.vcf}.format.vcf

new_sniffles_vcf=${outputdir}/${sniffles_vcf%.vcf}.format.vcf
new_cuteSV_vcf=${outputdir}/${cuteSV_vcf%.vcf}.format.vcf
new_pbsv_vcf=${outputdir}/${pbsv_vcf%.vcf}.format.vcf
new_svim_vcf=${outputdir}/${svim_vcf%.vcf}.format.vcf

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${new_sniffles_vcf} \
    --outfile ${new_sniffles_vcf%.format.vcf}.filter.vcf \
    --support 2

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${new_cuteSV_vcf} \
    --outfile ${new_cuteSV_vcf%.format.vcf}.filter.vcf \
    --support 2

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${new_pbsv_vcf} \
    --outfile ${new_pbsv_vcf%.format.vcf}.filter.vcf \
    --support 2

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${new_svim_vcf} \
    --outfile ${new_svim_vcf%.format.vcf}.filter.vcf \
    --support 2

# Clique merge
python3 Scripts/clique_maxflow_SV.py \
    --vcf_list ${new_sniffles_vcf%.format.vcf}.filter.vcf,${new_cuteSV_vcf%.format.vcf}.filter.vcf,${new_pbsv_vcf%.format.vcf}.filter.vcf,${new_svim_vcf%.format.vcf}.filter.vcf \
    --outfile ${outputdir}/${sample}.PDSP_common.vcf \
    --overlap_perc 0.5

python3 Scripts/MergeSVOriginal.py \
    --vcf ${outputdir}/${sample}.PDSP_common.vcf \
    --sniffle ${new_sniffles_vcf%.format.vcf}.filter.vcf \
    --cuteSV ${new_cuteSV_vcf%.format.vcf}.filter.vcf \
    --pbsv ${new_pbsv_vcf%.format.vcf}.filter.vcf \
    --svim ${new_svim_vcf%.format.vcf}.filter.vcf \
    --outfile  ${outputdir}/${sample}.PDSP_original.vcf \
    --method dominant

# Filtering
python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${outputdir}/${sample}.PDSP_original.vcf \
    --outfile ${outputdir}/${sample}.PDSP_filtering.vcf \
    --support 3

bcftools sort ${outputdir}/${sample}.PDSP_filtering.vcf > ${outputdir}/${sample}.PDSP_sort.vcf
bcftools view ${outputdir}/${sample}.PDSP_sort.vcf  -Oz -o ${outputdir}/${sample}.PDSP_sort.vcf.gz
bcftools index -t ${outputdir}/${sample}.PDSP_sort.vcf.gz

done

# Merging individual SVs

PDSP_VCF_list=PDSP_VCF_list.txt

#export PATH=~/Tools/bcftools-1.13:$PATH
python3 Scripts/Population_clique_maxflow_SV.py \
    --vcf_list ${PDSP_VCF_list}  \
    --outfile Results/Sample_SV_common_0.8.vcf \
    --allele_freq 0.2 \
    --overlap_perc 0.8
    
python3 Scripts/population_CAST_VCF_suppl.py \
    --vcf Results/Sample_SV_common_0.8.vcf \
    --vcf_list ${PDSP_VCF_list}  \
    --outfile Results/Sample_SV_common_0.8_suppl.vcf \
    --endlen 500
    
python3 Scripts/CombineSVTag.py \
    --vcf Results/Sample_SV_common_0.8_suppl.vcf  \
    --outfile /share/home1/Laisenying/Data-analysis/projects/PhageSV/PacBio_SV/CAST_merge/Population/Merge/Sample_common_SV_Tag_0.8.tsv



