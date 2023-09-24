import pandas as pd
import numpy as np
import argparse, operator, os, random, sys, time
import random
import os
from collections import defaultdict

def parseargs():    # handle user arguments
    # args.file='/home1/Laisenying/Data-analysis/projects/PhageSV/Test/pbsv/Mash_select/screen.tab'
    parser = argparse.ArgumentParser(description='cuteSV vcf reformat')
    parser.add_argument("--raw",help='./')  
    parser.add_argument("--new",help='./') 
    parser.add_argument("--score",default = 0, help='The score value threshold for SV.') 
    args = parser.parse_args()
    return args
    
def main():
    args=parseargs()
    # args.raw = '/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/pbsv/subread_mapping_subset.var.vcf'
    # args.new = '/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/pbsv.reformat.vcf'
    p_vcf = args.raw
    new_vcf = args.new   
    new_vcf_file = open(new_vcf,"w")  
    for line in open(p_vcf,"r").readlines():  
        if line.startswith('#'):
            new_vcf_file.write(line)
        else:
            string_list = line.strip('\n').split('\t')   
            scoreValue = string_list[5]     
            if scoreValue == ".":
                scoreValue = 0
            scoreValue = float(scoreValue)  
            if scoreValue >= args.score and string_list[6] == 'PASS':
                info_list = string_list[7].split(';')  
                if 'IMPRECISE' in info_list:
                    info_list.remove('IMPRECISE')
                if 'PRECISE' in info_list:
                    info_list.remove('PRECISE')  
                if 'SHADOWED' in info_list:
                    info_list.remove('SHADOWED')                                           
                info_dict = {}
                for i in info_list:
                    ni = i.split('=')
                    info_dict[ni[0]] = ni[1]    
                if info_dict['SVTYPE'] == 'BND':
                    continue
                if info_dict['SVTYPE'] == 'cnv':
                    continue
                if info_dict['SVTYPE'] == 'INV':
                    info_dict['SVLEN'] = str(int(info_dict['END']) - int(string_list[1]))                                                                              
                strands = '+-'
                chr2 = string_list[0]     
                format_genotype = 'GT:DR:DV'
                gt = {}
                for x,y in zip(string_list[-2].split(':'),string_list[-1].split(':')):
                    gt[x] = y  
                genotype = gt['GT']              
                if genotype == './.':
                    genotype = '0/0'                 
                dr = gt['AD'].split(',')[0]
                dv = gt['AD'].split(',')[1]                          
                new_line = string_list[0:7]
                new_info = ['SVTYPE=%s' % info_dict['SVTYPE'], 'END=%s' % info_dict['END'],'SUPPORT=%s' % dv,
                                    'SVLEN=%s' % info_dict['SVLEN'], 'CHR2=%s' % chr2, 'STRANDS=%s' % strands]
                new_line.append(';'.join(new_info))
                format_genotype = 'GT:DR:DV'
                new_line.append(format_genotype)
                new_line.append('%s:%s:%s' % (genotype, dr, dv))  
                new_vcf_file.write('\t'.join(new_line)+'\n')
    new_vcf_file.close()
    
if __name__ == '__main__':
    main()
