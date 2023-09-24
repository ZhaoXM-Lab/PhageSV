import pandas as pd
import numpy as np
import argparse, operator, os, random, sys, time
from Bio import SeqIO

def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Genome select')
    parser.add_argument("--outputfile",help='./')
    parser.add_argument("--NGS_VCF_list",help='./NGS_VCF')
    parser.add_argument("--PacBio_VCF_list",help='./PacBio_VCF')
    args = parser.parse_args()
    return args

def VCF_process(vcf_file):
    SV_stats = defaultdict(list)
    for line in open(vcf_file,"r"):
        if line.startswith('#'):
            continue
        else:
            lines = line.strip().split()
            info=lines[7]
            info_dict={}
            for x in info.split(';'):
                if len(x.split('=')) == 2:
                    info_dict[x.split('=')[0]] = x.split('=')[1]
            if info_dict['SVTYPE'] == 'BND':
                continue
            SV_stats['SVTYPE'].append(info_dict['SVTYPE'])
            SV_stats['contig'].append(lines[0])
            SV_stats['Start'].append(int(lines[1]))
            SV_stats['SVID'].append(lines[2])
            SV_stats['FILTER'].append(lines[6])
            SV_stats['END'].append(info_dict['END'])
    return SV_stats
    
def main():
    args = parseargs()
    NGS_SV_stats = defaultdict(list)
    PB_SV_stats = defaultdict(list)
    for NGS_file in list(pd.read_csv(args.NGS_VCF_list,header=None)[0]):
        sample = NGS_file.split("/")[7]
        NGS_stats =  pd.DataFrame(VCF_process(NGS_file))
        NGS_stats = NGS_stats.loc[NGS_stats['FILTER']=='PASS',]
        NGS_SV_stat = pd.value_counts(NGS_stats['SVTYPE'])
        if 'DEL' in NGS_SV_stat.index:
            NGS_SV_stats['DEL'].append(NGS_SV_stat['DEL'])
        else:
            NGS_SV_stats['DEL'].append(0)
        if 'INS' in NGS_SV_stat.index:
            NGS_SV_stats['INS'].append(NGS_SV_stat['INS'])
        else:
            NGS_SV_stats['INS'].append(0)
        if 'DUP' in NGS_SV_stat.index:
            NGS_SV_stats['DUP'].append(NGS_SV_stat['DUP'])
        else:
            NGS_SV_stats['DUP'].append(0)
        if 'INV' in NGS_SV_stat.index:
            NGS_SV_stats['INV'].append(NGS_SV_stat['INV'])
        else:
            NGS_SV_stats['INV'].append(0)
        NGS_SV_stats['sample'].append(sample)
    for PB_file in list(pd.read_csv(args.PacBio_VCF_list,header=None)[0]):
        sample = PB_file.split("/")[7]
        PB_stats =  pd.DataFrame(VCF_process(PB_file))
        PB_stats = PB_stats.loc[PB_stats['FILTER']=='PASS',]
        PB_SV_stat = pd.value_counts(PB_stats['SVTYPE'])
        if 'DEL' in PB_SV_stat.index:
            PB_SV_stats['DEL'].append(PB_SV_stat['DEL'])
        else:
            PB_SV_stats['DEL'].append(0)
        if 'INS' in PB_SV_stat.index:
            PB_SV_stats['INS'].append(PB_SV_stat['INS'])
        else:
            PB_SV_stats['INS'].append(0)
        if 'DUP' in PB_SV_stat.index:
            PB_SV_stats['DUP'].append(PB_SV_stat['DUP'])
        else:
            PB_SV_stats['DUP'].append(0)
        if 'INV' in PB_SV_stat.index:
            PB_SV_stats['INV'].append(PB_SV_stat['INV'])
        else:
            PB_SV_stats['INV'].append(0)
        PB_SV_stats['sample'].append(sample)
    PB_SV_data = pd.DataFrame(PB_SV_stats)
    PB_SV_data['Group'] = 'PacBio'
    NGS_SV_data = pd.DataFrame(NGS_SV_stats)
    NGS_SV_data['Group'] = 'Illumina'
    SV_data = pd.concat([PB_SV_data,NGS_SV_data])
    SV_data['Total'] = SV_data['DEL'] + SV_data['INS'] + SV_data['DUP'] + SV_data['INV']
    SV_data.to_csv(args.outputfile)

if __name__=='__main__':
    main()

    
