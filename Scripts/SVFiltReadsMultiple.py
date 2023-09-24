import pandas as pd
import numpy as np
import argparse, operator, os, random, sys, time
import random
import os
from collections import defaultdict

def parseargs():    # handle user arguments
    # args.file='/home1/Laisenying/Data-analysis/projects/PhageSV/Test/pbsv/Mash_select/screen.tab'
    parser = argparse.ArgumentParser(description='Filt SV file based on the number of support records.')
    parser.add_argument("--vcf",help='./')  
    parser.add_argument("--outfile",help='The output file.') 
    parser.add_argument("--support", type=int, default = 3, help='The read support number of SV.')     
    args = parser.parse_args()
    return args
    
def Infor_target_values(Infor, target):
    """
    Infor = "IMPRECISE;SVMETHOD=Snifflesv1.0.10;CHR2=1;END=181364;STD_quant_start=26.000000;
    STD_quant_stop=79.501572;Kurtosis_quant_start=-2.000000;Kurtosis_quant_stop=-1.999842;
    SVTYPE=DEL;RNAMES=a5c2a7ee-ce33-4dd3-8fac-3a18286ce747,b0fdea87-ced4-44a6-b6dc-b9071015fac0;SUPTYPE=AL;SVLEN=-182;STRANDS=+-;RE=2"
    target = "END"
    return:
    '181364'
    target = "END,SVTYPE"
    return:
    ["181364", "DEL"]
    """
    TagValues = {}
    Infors = Infor.split(";")
    for f in Infors:
        fs = f.split("=") 
        if len(fs) == 1:
            if fs[0] == "PRECISE":
                TagValues["PRECISE"] = True
            else:
                TagValues["PRECISE"] = False
        elif len(fs) == 2:
            tag, value = fs
            TagValues[tag] = value
        else:
            print("Please check and make sure the items of tag is no more than two.")
            sys.exit(1)
    ### get the target items
    Items = []
    targets = target.split(",")
    targets = [t.strip() for t in targets]
    for t in targets:
        try:
            value = TagValues[t]
            Items.append(value)
        except KeyError:
            print("Please check and make sure the given tag %s is in the Infor record %s." % (t, Infor))
            Items.append("0")
    if len(Items) == 1:
        return Items[0]
    else:
        return Items    
    
def main():
    args=parseargs()
    sv_file = args.vcf # args.vcf = '/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/Sniffles.reformat.vcf'
    out_file = args.outfile   
    out_h = open(out_file,'w') 
    for line in open(sv_file,'r').readlines():
        if line.startswith("#"):
            out_h.write(line)   
        else:
            string_list = line.strip('\n').split()  
            gt = {}
            for x,y in zip(string_list[8].split(':'),string_list[9].split(':')):
                gt[x] = y       
            DR = int(gt['DR'])
            DV = int(gt['DV'])
            RE = DV
            if RE >= args.support:
                out_h.write(line)
            else:
                #print(line)
                continue  
    out_h.close()                

if __name__ == '__main__':
    main()
