#########################################################################################
## Benchmarking SV calling approaches on simualted virome-enriched metagenomics datasets
#########################################################################################

## Step 1. Using sim-it to simulate viral SVs on randomly selected 1000 CHGV-HQ viral genomes

git clone https://github.com/ndierckx/Sim-it.git
cd Sim-it
# Setting: Deletions = Insertions = 3000, Tandem duplications = 1000, Inversions = 500, Inverted duplications = 200
# Randomly select 1000 CHGV-HQ genomes
mkdir -p outputdir
perl Sim-it1.3.2.pl -c config_Sim-it.txt -o ./outputdir  # output: ./outputdir/Genome_withSV.fasta, ./outputdir/Genome_withSV.vcf


## Step 2. Abundance simulation under two conditions: even and uneven distribution
python3 Scripts/Genome_write.py \
    --SimSVdir ./outputdir # output: genome_list.txt, path for each selected genome
    
# Uneven condition
MGSIM communities \
    ./outputdir/genome_list.txt \
    ./outputdir/Abundance_uneven \
    --n-comm 1 \
    --abund-dist-p mean:10,sigma:2 \
    --abund-dist lognormal \
    --perm-perc 0.5 \
    --richness 1

# Even condition
MGSIM communities \
    ./outputdir/genome_list.txt \
    ./outputdir/Abundance_even \
    --n-comm 1 \
    --abund-dist-p mean:10,sigma:0 \
    --abund-dist lognormal \
    --perm-perc 0.5 \
    --richness 1

## 3. PacBio simulation
# PacBio Simulating with even sequencing depth
python3 Scripts/PacBio_sim.py \
    --genomefile ./outputdir/genome_list.txt \
    --abund_file ./outputdir/Abundance_even_abund.txt \
    --outputdir ./outputdir/Even \
    --seq_depth 6e5


# PacBio Simulating with variable sequencing depth
python3 Scripts/PacBio_sim.py \
    --genomefile ./outputdir/genome_list.txt \
    --abund_file ./outputdir/Abundance_uneven_abund.txt \
    --outputdir ./outputdir/Uneven \
    --seq_depth 6e5

## Step 3. Generate mapping file
PacBio_even_read=./outputdir/Even/pacbio/1/R1.fq
PacBio_uneven_read=./outputdir/Uneven/pacbio/1/R1.fq

RefGenome=Select_1000genome.fasta
ViralGenome=HGV.hq.genome.fa

# Mapping reads to all viral genomes
pbmm2 align ${ViralGenome} ${PacBio_even_read} \
    ./outputdir/Even/PacBio.map_Nofilter.bam --sort --preset CCS --sample Nofilter
    
pbmm2 align ${ViralGenome} ${PacBio_uneven_read} \
    ./outputdir/Uneven/PacBio.map_Nofilter.bam --sort --preset CCS --sample Nofilter

# Mapping reads to selected 1000 viral genomes
pbmm2 align ${RefGenome} ${PacBio_even_read} \
    ./outputdir/Even/PacBio.map_Prefilter.bam --sort --preset CCS --sample Prefilter
    
pbmm2 align ${RefGenome} ${PacBio_uneven_read} \
    ./outputdir/Uneven/PacBio.map_Prefilter.bam --sort --preset CCS --sample Prefilter

## Step 4. SV calling for all approaches

### Under Uneven condition ###
condition=Uneven
# (1) PBSV
# output1: ./outputdir/${condition}/pbsv/PacBio.map_Nofilter.var.vcf
# output2: ./outputdir/${condition}/pbsv/PacBio.map_Prefilter.var.vcf
mkdir -p ./outputdir/${condition}/pbsv
pbsv discover \
    ./outputdir/${condition}/PacBio.map_Nofilter.bam \
    ./outputdir/${condition}/pbsv/PacBio.map_Nofilter.svsig.gz
pbsv call ${ViralGenome} ./outputdir/Uneven/pbsv/PacBio.map_Nofilter.svsig.gz \
    ./outputdir/${condition}/pbsv/PacBio.map_Nofilter.var.vcf

pbsv discover \
    ./outputdir/${condition}/PacBio.map_Prefilter.bam \
    ./outputdir/${condition}/pbsv/PacBio.map_Prefilter.svsig.gz
pbsv call ${RefGenome} ./outputdir/${condition}/pbsv/PacBio.map_Prefilter.svsig.gz \
    ./outputdir/${condition}/pbsv/PacBio.map_Prefilter.var.vcf

# (2) cuteSV
# output1: ./outputdir/${condition}/cuteSV/cuteSV.map_Nofilter.var.vcf
# output2: ./outputdir/${condition}/cuteSV/cuteSV.map_Prefilter.var.vcf
mkdir -p ./outputdir/${condition}/cuteSV
cuteSV ./outputdir/${condition}/PacBio.map_Nofilter.bam \
    ${ViralGenome} \
    ./outputdir/${condition}/cuteSV/cuteSV.map_Nofilter.var.vcf \
    ./outputdir/${condition}/cuteSV \
    --max_cluster_bias_INS 1000 \
    --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5

cuteSV ./outputdir/${condition}/PacBio.map_Prefilter.bam \
    ${RefGenome} \
    ./outputdir/${condition}/cuteSV/cuteSV.map_Prefilter.var.vcf \
    ./outputdir/${condition}/cuteSV \
    --max_cluster_bias_INS 1000 \
    --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5

# (3) SVIM
# output1: ./outputdir/${condition}/SVIM/Nofilter/variants.vcf
# output2: ./outputdir/${condition}/SVIM/Prefilter/variants.vcf

mkdir -p ./outputdir/${condition}/SVIM/Nofilter
mkdir -p ./outputdir/${condition}/SVIM/Prefilter

svim alignment ./outputdir/${condition}/SVIM/Nofilter \
    ./outputdir/${condition}/PacBio.map_Nofilter.bam \
     ${ViralGenome}

svim alignment ./outputdir/${condition}/SVIM/Prefilter \
    ./outputdir/${condition}/PacBio.map_Prefilter.bam \
     ${RefGenome}

# (4) Sniffles
mkdir -p ./outputdir/${condition}/Sniffles

sniffles -i ./outputdir/${condition}/PacBio.map_Nofilter.bam \
    -v ./outputdir/${condition}/Sniffles/Sniffles.map_Nofilter.var.vcf

sniffles -i ./outputdir/${condition}/PacBio.map_Prefilter.bam \
    -v ./outputdir/${condition}/Sniffles/Sniffles.map_Prefilter.var.vcf


# (5) PSDP
# output1: ./outputdir/${condition}/PSDP/PSDP.map_Prefilter.sort.vcf
# output2: ./outputdir/${condition}/PSDP/PSDP.map_Nofilter.sort.vcf

sniffles_Nofilter_vcf=./outputdir/${condition}/Sniffles/Sniffles.map_Nofilter.var.vcf
cuteSV_Nofilter_vcf=./outputdir/${condition}/cuteSV/cuteSV.map_Nofilter.var.vcf
pbsv_Nofilter_vcf=./outputdir/${condition}/pbsv/PacBio.map_Nofilter.var.vcf
svim_Nofilter_vcf=./outputdir/${condition}/SVIM/Nofilter/variants.vcf

sniffles_Prefilter_vcf=./outputdir/${condition}/Sniffles/Sniffles.map_Prefilter.var.vcf
cuteSV_Prefilter_vcf=./outputdir/${condition}/cuteSV/cuteSV.map_Prefilter.var.vcf
pbsv_Prefilter_vcf=./outputdir/${condition}/pbsv/PacBio.map_Prefilter.var.vcf
svim_Prefilter_vcf=./outputdir/${condition}/SVIM/Prefilter/variants.vcf

mkdir -p ./outputdir/${condition}/PSDP

# Formating VCF files
python3 Scripts/Sniffles_change_format.py \
    --raw ${sniffles_Nofilter_vcf} \
    --new ./outputdir/${condition}/PSDP/Sniffles.map_Nofilter.reformat.vcf

python3 Scripts/cuteSV_change_format.py \
    --raw ${cuteSV_Nofilter_vcf} \
    --new ./outputdir/${condition}/cuteSV/new_vcf/cuteSV.map_Nofilter.reformat.vcf
    
python3 Scripts/pbsv_change_format.py \
    --raw ${pbsv_Nofilter_vcf} \
    --new ./outputdir/${condition}/pbsv/new_vcf/pbsv.map_Nofilter.reformat.vcf

python3 Scripts/svim_change_format.py \
    --raw ${svim_Nofilter_vcf} \
    --new ./outputdir/${condition}/SVIM/svim.map_Nofilter.reformat.vcf
    
python3 Scripts/Sniffles_change_format.py \
    --raw ${sniffles_Prefilter_vcf} \
    --new ./outputdir/${condition}/PSDP/Sniffles.map_Prefilter.reformat.vcf

python3 Scripts/cuteSV_change_format.py \
    --raw ${cuteSV_Prefilter_vcf} \
    --new ./outputdir/${condition}/cuteSV/cuteSV.map_Prefilter.reformat.vcf
    
python3 Scripts/pbsv_change_format.py \
    --raw ${pbsv_Prefilter_vcf} \
    --new ./outputdir/${condition}/pbsv/pbsv.map_Prefilter.reformat.vcf

python3 Scripts/svim_change_format.py \
    --raw ${svim_Nofilter_vcf} \
    --new ./outputdir/${condition}/SVIM/svim.map_Prefilter.reformat.vcf

sniffles_Nofilter_vcf=./outputdir/${condition}/Sniffles/Sniffles.map_Nofilter.reformat.vcf
cuteSV_Nofilter_vcf=./outputdir/${condition}/cuteSV/cuteSV.map_Nofilter.reformat.vcf
pbsv_Nofilter_vcf=./outputdir/${condition}/pbsv/PacBio.map_Nofilter.reformat.vcf
svim_Nofilter_vcf=./outputdir/${condition}/SVIM/svim.map_Nofilter.reformat.vcf

sniffles_Prefilter_vcf=./outputdir/${condition}/Sniffles/Sniffles.map_Prefilter.reformat.vcf
cuteSV_Prefilter_vcf=./outputdir/${condition}/cuteSV/cuteSV.map_Prefilter.reformat.vcf
pbsv_Prefilter_vcf=./outputdir/${condition}/pbsv/PacBio.map_Prefilter.reformat.vcf
svim_Prefilter_vcf=./outputdir/${condition}/SVIM/svim.map_Prefilter.reformat.vcf

# Filtering SVs with < 2 supported reads
python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${sniffles_Nofilter_vcf} \
    --outfile ${sniffles_Nofilter_vcf%.reformat.vcf}.filter.vcf \
    --support 2

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${cuteSV_Nofilter_vcf} \
    --outfile ${cuteSV_Nofilter_vcf%.reformat.vcf}.filter.vcf \
    --support 2

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${pbsv_Nofilter_vcf} \
    --outfile ${pbsv_Nofilter_vcf%.reformat.vcf}.filter.vcf \
    --support 2

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${svim_Nofilter_vcf} \
    --outfile ${svim_Nofilter_vcf%.reformat.vcf}.filter.vcf \
    --support 2

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${sniffles_Prefilter_vcf} \
    --outfile ${sniffles_Prefilter_vcf%.reformat.vcf}.filter.vcf \
    --support 2

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${cuteSV_Prefilter_vcf} \
    --outfile ${cuteSV_Prefilter_vcf%.reformat.vcf}.filter.vcf \
    --support 2

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${pbsv_Prefilter_vcf} \
    --outfile ${pbsv_Prefilter_vcf%.reformat.vcf}.filter.vcf \
    --support 2

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ${svim_Prefilter_vcf} \
    --outfile ${svim_Prefilter_vcf%.reformat.vcf}.filter.vcf \
    --support 2

# Clique SV merging
python3 Scripts/clique_maxflow_SV.py \
    --vcf_list ${sniffles_Nofilter_vcf%.reformat.vcf}.filter.vcf,${cuteSV_Nofilter_vcf%.reformat.vcf}.filter.vcf,${pbsv_Nofilter_vcf%.reformat.vcf}.filter.vcf,${svim_Nofilter_vcf%.reformat.vcf}.filter.vcf \
    --outfile ./outputdir/${condition}/PSDP/PSDP.map_Nofilter.common.vcf

python3 Scripts/MergeSVOriginal.py \
    --vcf ./outputdir/${condition}/PSDP/PSDP.map_Nofilter.common.vcf \
    --sniffle ${sniffles_Nofilter_vcf%.reformat.vcf}.filter.vcf \
    --cuteSV ${cuteSV_Nofilter_vcf%.reformat.vcf}.filter.vcf \
    --pbsv ${pbsv_Nofilter_vcf%.reformat.vcf}.filter.vcf \
    --svim ${svim_Nofilter_vcf%.reformat.vcf}.filter.vcf \
    --outfile  ./outputdir/${condition}/PSDP/PSDP.map_Nofilter.original.vcf \
    --method dominant
    
python3 Scripts/clique_maxflow_SV.py \
    --vcf_list ${sniffles_Prefilter_vcf%.reformat.vcf}.filter.vcf,${cuteSV_Prefilter_vcf%.reformat.vcf}.filter.vcf,${pbsv_Prefilter_vcf%.reformat.vcf}.filter.vcf,${svim_Prefilter_vcf%.reformat.vcf}.filter.vcf \
    --outfile ./outputdir/${condition}/PSDP/PSDP.map_Prefilter.common.vcf

python3 Scripts/MergeSVOriginal.py \
    --vcf ./outputdir/${condition}/PSDP/PSDP.map_Prefilter.common.vcf \
    --sniffle ${sniffles_Prefilter_vcf%.reformat.vcf}.filter.vcf \
    --cuteSV ${cuteSV_Prefilter_vcf%.reformat.vcf}.filter.vcf \
    --pbsv ${pbsv_Prefilter_vcf%.reformat.vcf}.filter.vcf \
    --outfile  ./outputdir/${condition}/PSDP/PSDP.map_Prefilter.original.vcf \
    --method dominant
    
    
# Retaining SVs with at least three supported reads
python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ./outputdir/${condition}/PSDP/PSDP.map_Prefilter.original.vcf \
    --outfile ./outputdir/${condition}/PSDP/PSDP.map_Prefilter.filtering.vcf \
    --support 3

python3 Scripts/SVFiltReadsMultiple.py \
    --vcf ./outputdir/${condition}/PSDP/PSDP.map_Nofilter.original.vcf \
    --outfile ./outputdir/${condition}/PSDP/PSDP.map_Nofilter.filtering.vcf \
    --support 3

bcftools sort ./outputdir/${condition}/PSDP/PSDP.map_Nofilter.filtering.vcf > ./outputdir/${condition}/PSDP/PSDP.map_Nofilter.sort.vcf
bcftools sort ./outputdir/${condition}/PSDP/PSDP.map_Prefilter.filtering.vcf > ./outputdir/${condition}/PSDP/PSDP.map_Prefilter.sort.vcf

## Step 5 Performance evaluation
# Formating Simulated VCF file
SimSV_file=./outputdir/Genome_withSV.vcf
Sim_VCF_file=./outputdir/Genome_withSV_format.vcf
reference_genome=Select_1000genome.fasta

python3 Scripts/SimConvVCF.py \
    --inputfile ${SimSV_file} \
    --outputfile ${Sim_VCF_file} \
    --reference ${reference_genome}

True_VCF_file=./outputdir/Genome_withSV_format.vcf

# 1. PBSV
Com_VCF_file=./outputdir/${condition}/pbsv/PacBio.map_Nofilter.var.vcf
bcftools view ${Com_VCF_file} -Oz -o ${Com_VCF_file}.gz
bcftools index -t ${Com_VCF_file}.gz
truvari bench -b ${True_VCF_file}.gz -c ${Com_VCF_file}.gz -f ${reference_genome} -o ./outputdir/${condition}/pbsv/truvari_Nofilter --passonly -p 0 --pctsize 0 --pctovl 0

Com_VCF_file=./outputdir/${condition}/pbsv/PacBio.map_Prefilter.var.vcf
bcftools view ${Com_VCF_file} -Oz -o ${Com_VCF_file}.gz
bcftools index -t ${Com_VCF_file}.gz
truvari bench -b ${True_VCF_file}.gz -c ${Com_VCF_file}.gz -f ${reference_genome} -o ./outputdir/${condition}/pbsv/truvari_Prefilter --passonly -p 0 --pctsize 0 --pctovl 0

# 2. cuteSV
Com_VCF_file=./outputdir/${condition}/cuteSV/cuteSV.map_Nofilter.var.vcf
bcftools view ${Com_VCF_file} -Oz -o ${Com_VCF_file}.gz
bcftools index -t ${Com_VCF_file}.gz
truvari bench -b ${True_VCF_file}.gz -c ${Com_VCF_file}.gz -f ${reference_genome} -o ./outputdir/${condition}/cuteSV/truvari_Nofilter --passonly -p 0 --pctsize 0 --pctovl 0

Com_VCF_file=./outputdir/${condition}/cuteSV/cuteSV.map_Prefilter.var.vcf
bcftools view ${Com_VCF_file} -Oz -o ${Com_VCF_file}.gz
bcftools index -t ${Com_VCF_file}.gz
truvari bench -b ${True_VCF_file}.gz -c ${Com_VCF_file}.gz -f ${reference_genome} -o ./outputdir/${condition}/cuteSV/truvari_Prefilter --passonly -p 0 --pctsize 0 --pctovl 0

# 3. SVIM
Com_VCF_file=./outputdir/${condition}/SVIM/Nofilter/variants.vcf
bcftools view ${Com_VCF_file} -Oz -o ${Com_VCF_file}.gz
bcftools index -t ${Com_VCF_file}.gz
truvari bench -b ${True_VCF_file}.gz -c ${Com_VCF_file}.gz -f ${reference_genome} -o ./outputdir/${condition}/SVIM/truvari_Nofilter --passonly -p 0 --pctsize 0 --pctovl 0

Com_VCF_file=./outputdir/${condition}/SVIM/Prefilter/variants.vcf
bcftools view ${Com_VCF_file} -Oz -o ${Com_VCF_file}.gz
bcftools index -t ${Com_VCF_file}.gz
truvari bench -b ${True_VCF_file}.gz -c ${Com_VCF_file}.gz -f ${reference_genome} -o ./outputdir/${condition}/SVIM/truvari_Prefilter --passonly -p 0 --pctsize 0 --pctovl 0

# 4. Sniffles
Com_VCF_file=./outputdir/${condition}/Sniffles/Sniffles.map_Nofilter.var.vcf
bcftools view ${Com_VCF_file} -Oz -o ${Com_VCF_file}.gz
bcftools index -t ${Com_VCF_file}.gz
truvari bench -b ${True_VCF_file}.gz -c ${Com_VCF_file}.gz -f ${reference_genome} -o ./outputdir/${condition}/Sniffles/truvari_Nofilter --passonly -p 0 --pctsize 0 --pctovl 0

Com_VCF_file=./outputdir/${condition}/Sniffles/Sniffles.map_Prefilter.var.vcf
bcftools view ${Com_VCF_file} -Oz -o ${Com_VCF_file}.gz
bcftools index -t ${Com_VCF_file}.gz
truvari bench -b ${True_VCF_file}.gz -c ${Com_VCF_file}.gz -f ${reference_genome} -o ./outputdir/${condition}/Sniffles/truvari_Prefilter --passonly -p 0 --pctsize 0 --pctovl 0


# 5. PSDP
Com_VCF_file=./outputdir/${condition}/PSDP/PSDP.map_Nofilter.sort.vcf
bcftools view ${Com_VCF_file} -Oz -o ${Com_VCF_file}.gz
bcftools index -t ${Com_VCF_file}.gz
truvari bench -b ${True_VCF_file}.gz -c ${Com_VCF_file}.gz -f ${reference_genome} -o ./outputdir/${condition}/PSDP/truvari_Nofilter --passonly -p 0 --pctsize 0 --pctovl 0

Com_VCF_file=./outputdir/${condition}/PSDP/PSDP.map_Prefilter.sort.vcf
bcftools view ${Com_VCF_file} -Oz -o ${Com_VCF_file}.gz
bcftools index -t ${Com_VCF_file}.gz
truvari bench -b ${True_VCF_file}.gz -c ${Com_VCF_file}.gz -f ${reference_genome} -o ./outputdir/${condition}/PSDP/truvari_Prefilter --passonly -p 0 --pctsize 0 --pctovl 0

