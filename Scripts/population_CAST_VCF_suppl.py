from collections import defaultdict, Counter
import multiprocessing as mp
import networkx as nx
import numpy as np
import argparse
import time
import re
import operator

def parseargs():    # handle user arguments
    # args.file='/home1/Laisenying/Data-analysis/projects/PhageSV/Test/pbsv/Mash_select/screen.tab'
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument("--vcf",type=str,default='',help='The input merged vcf file.')  
    parser.add_argument("--vcf_list",type=str,default='',help='The input merged vcf file.')  
    parser.add_argument("--outfile",help='path to the output vcf') 
    parser.add_argument("--endlen",default=100, type=int, help='path to the output vcf')                
    args = parser.parse_args()
    return args
    
def Infor_target_values(Infor, target):
    """
    Infor = "IMPRECISE;SVMETHOD=Snifflesv1.0.10;CHR2=1;END=181364;STD_quant_start=26.000000;STD_quant_stop=79.501572;Kurtosis_quant_start=-2.000000;Kurtosis_quant_stop=-1.999842;SVTYPE=DEL;RNAMES=a5c2a7ee-ce33-4dd3-8fac-3a18286ce747,b0fdea87-ced4-44a6-b6dc-b9071015fac0;SUPTYPE=AL;SVLEN=-182;STRANDS=+-;RE=2"

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
    
def tag_record(sv_file):
    TagRecord = {}
    sv_h = open(sv_file, "r")
    for line in sv_h:
        line = line.strip()
        if line.startswith("#"):
            continue
        else:
            lines = line.split("\t")
            Chr, Start = lines[:2]
            Infor = lines[7]
            Chr2, End, Length, Type = Infor_target_values(Infor, "CHR2,END,SVLEN,SVTYPE")
            Length = Length.lstrip("-")
            Tag = "%s~%s--%s~%s--%s--%s" % (Chr, Start, Chr2, End, Length, Type)
            TagRecord[Tag] = line
    sv_h.close()
    return TagRecord  

def parse_genotype_format(Format, Value, target):
    """
    Format = "GT:AD:DP:GQ:PL"
    Value = "1/1:0,2:2:6:90"
    target = "AD"
    return '0,2'

    target = "AD,PL"
    ['0,2', '90']
    """
    Formats = Format.split(":")
    Values = Value.split(":")
    if len(Formats) != len(Values):
        print("The number of ids and values is not identical for %s and %s." % (Format, Value))
        return None
    else:
        targets = target.split(",")
        targetValues = []
        for t in targets:
            try:
                valueIndex = Formats.index(t)
            except ValueError:
                print("Please check whether the target %s in %s." % (t, Format))
                sys.exit(1)
            v = Values[valueIndex]
            targetValues.append(v)
        ### return the values
        if len(targetValues) == 1:
            return targetValues[0]
        else:
            return targetValues   
                    
def nearst_region(Format, Start, End, geno):
    GenoDistance = {}
    genos = geno.split(";")
    genoLen = len(genos)
    if genoLen == 1:
        Tag1 = parse_genotype_format(Format, genos[0], "CO")
        Length = parse_genotype_format(Format, genos[0], "LN")
        Length = Length.lstrip("-")
        Type = parse_genotype_format(Format, genos[0], "TY")
    else:
        for g in genos:
            CO = parse_genotype_format(Format, g, "CO")
            regions = CO.split("--")
            Start2 =  int(regions[0].split("~")[-1])
            End2 =  int(regions[1].split("~")[-1])
            distance = abs(Start - Start2) + abs(End - End2)
            GenoDistance[g] = distance
        sortedDistance = sorted(GenoDistance.items(), key=operator.itemgetter(1))
        targetGeno = sortedDistance[0][0]
        Tag1 = parse_genotype_format(Format, targetGeno, "CO")
        Length = parse_genotype_format(Format, targetGeno, "LN")
        Length = Length.lstrip("-")
        Type = parse_genotype_format(Format, targetGeno, "TY")
    Tag = "%s--%s--%s" % (Tag1, Length, Type)
    return Tag        

def nearst_Tag(Format, Start, End, genos,samples):
    GenoDistance = {}
    allgenos = {}
    for sample,geno in zip(samples,genos):
        for x in geno.split(";"):
            allgenos[sample] = x     
    for sample,g in allgenos.items():
        CO = parse_genotype_format(Format, g, "CO")
        regions = CO.split("--")
        Start2 =  int(regions[0].split("~")[-1])
        End2 =  int(regions[1].split("~")[-1])
        distance = abs(Start - Start2) + abs(End - End2)
        GenoDistance[g] = (sample,distance)
    targetGeno = np.array(list(GenoDistance.keys()))[np.argsort([x[1][1] for x in GenoDistance.items()])][0]  
    targetSamp = np.array(list(GenoDistance.values()))[np.argsort([x[1][1] for x in GenoDistance.items()])][0][0]  
    #sortedDistance = sorted(GenoDistance.items(), key=operator.itemgetter(1))
    #targetGeno = sortedDistance[0][0]
    Tag1 = parse_genotype_format(Format, targetGeno, "CO")
    Length = parse_genotype_format(Format, targetGeno, "LN")
    Length = Length.lstrip("-")
    Type = parse_genotype_format(Format, targetGeno, "TY")
    Tag = "%s--%s--%s" % (Tag1, Length, Type)
    return Tag, targetSamp 

def vcf_process(vcf_list):
    vcf_lists = {}
    for line in open(vcf_list,"r"):
        vcf_file = line.strip()
        sample = line.strip().split('/')[-3]
        vcf_lists[sample] = vcf_file     
    vcfrecord_dict = {}                 
    for sample, vcf in vcf_lists.items():
        vcfrecord =  tag_record(vcf) 
        vcfrecord_dict[sample] = vcfrecord
    return vcfrecord_dict  
    
                                    
def main():
    args=parseargs()    
    # args.vcf='/home1/Laisenying/Data-analysis/projects/PhageSV/PacBio_SV/CAST_merge/Population/Merge/Sample_SV_common.vcf'
    # args.vcf_list = '/home1/Laisenying/Data-analysis/projects/PhageSV/PacBio_SV/CAST_merge/Population/Merge/VCF_list.txt'
    # args.outfile='/home1/Laisenying/Data-analysis/projects/PhageSV/PacBio_SV/CAST_merge/Population/Merge/Sample_SV_common_format.vcf'
    vcfrecord = vcf_process(args.vcf_list)
    outfile = open(args.outfile,"w")
    contig_length = {}
    for line in open(args.vcf,"r"):
        if line.startswith("##"):
            if line.startswith('##contig='):
                contig_length[line.split('ID=')[1].split(',')[0]] = int(line.split('length=')[1].split('>')[0])                
            outfile.write(line)
        elif line.startswith('#'):
            lines = line.strip("\n").split()   
            sampleNames = np.array(lines[9:])  
            outfile.write(line)                      
        else:
            lines = line.strip("\n").split()  
            Chr = lines[0] 
            contig_len = contig_length[Chr]
            Start = int(lines[1])  
            info = lines[7]  
            info_dict = {}
            Format = lines[8]
            for s in info.split(";"):
                info_dict[s.split("=")[0]] = s.split("=")[1]
            End = int(info_dict['END'])  
            if 'SVLEN' in info_dict:
                svlen = int(info_dict['SVLEN'])
            elif 'AVGLEN' in info_dict:
                svlen = int(info_dict['AVGLEN'])
            if abs(svlen) < 50:
                continue   
            if Start <= args.endlen:
                continue
            if End >= (contig_len - args.endlen):
                continue                                                                                                              
            if info_dict['SVTYPE'] == 'INV' or info_dict['SVTYPE'] == 'DUP':
                outfile.write(line)
            else:
                sample_genos = lines[9:]  
                GenoTags = np.array(["-" if g == "./.:NaN:0:0,0:--:NaN:NaN" else g for g in sample_genos])
                samples = sampleNames[GenoTags != '-']  
                genos =  GenoTags[GenoTags != '-']
                if len(samples) == 1:
                    sample = samples[0] 
                    geno = genos[0] 
                    Tag = nearst_region(Format, Start, End, geno)      
                    sample_vcf_file = vcfrecord[sample]   
                    record = sample_vcf_file[Tag]   
                    REF = record.split("\t")[3]
                    ALT = record.split("\t")[4]  
                    lines[3] = REF
                    lines[4] = ALT
                    record = '\t'.join(lines)                      
                else:
                    Tag = nearst_Tag(Format, Start, End, genos,samples)
                    sample_vcf_file = vcfrecord[sample] 
                    Tag, sample = nearst_Tag(Format, Start, End, genos,samples)    
                    sample_vcf_file = vcfrecord[sample]  
                    record = sample_vcf_file[Tag]  
                    REF = record.split("\t")[3]
                    ALT = record.split("\t")[4] 
                    lines[3] = REF
                    lines[4] = ALT
                    record = '\t'.join(lines)                                                                                                                                                                                                                                                                         
                outfile.write("%s\n" % record) 
    outfile.close()

if __name__ == "__main__":
    main()
