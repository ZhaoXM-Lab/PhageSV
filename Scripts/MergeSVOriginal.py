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
    parser.add_argument("--outfile",help='path to the output vcf')     
    parser.add_argument("--method",help="The method, such as 'common', 'dominant' or 'denovo'.")
    parser.add_argument("--sniffle",help='path to sniffle SV file') 
    parser.add_argument("--cuteSV",help='path to cuteSV SV file')
    parser.add_argument("--pbsv",help='path to pbsv SV file')  
    parser.add_argument("--svim",help='path to pbsv SV file')            
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
    
def merge_SV_original(merged_SV, sniffle_file, cuteSV_file, pbsv_file,svim_file,out_file, method):       
    SniffleRecord = tag_record(sniffle_file)  
    out_h = open(out_file, "w")
    merge_h = open(merged_SV, "r")  
    method = method.lower()
    if method == 'dominant':
        cuteSVRecord = tag_record(cuteSV_file)
        pbsvRecord = tag_record(pbsv_file)
        svimRecord = tag_record(svim_file)
        for line in merge_h:
            if line.startswith("##"):
                out_h.write(line)                
            elif line.startswith("#"):
                out_h.write("%s\n" % "\t".join(line.split("\t")[:-3]))   
            else:
                lines = line.strip().split('\t') 
                Chr = lines[0]
                Start = int(lines[1])
                Infor = lines[7]
                info_dict={}
                for s in Infor.split(";"):
                    info_dict[s.split("=")[0]] = s.split("=")[1]  
                svtype = info_dict['SVTYPE']                                      
                End = int(Infor_target_values(Infor, "END"))
                Format = lines[8]
                Genotypes = lines[9:]
                genos = [parse_genotype_format(Format, g.split(";")[0], "GT") for g in Genotypes]
                zeroNum = genos.count("./.")
                nonzeroNum = len(genos) - zeroNum
                if nonzeroNum >= 2:
                    if genos[0] != "./.":
                        Tag = nearst_region(Format, Start, End, Genotypes[0])                                                
                        try:
                            record = SniffleRecord[Tag]
                        except KeyError:
                            print("Please check whether the tag %s is in file %s." % (Tag, sniffle_file))
                            #sys.exit(1)   
                    elif genos[1] != "./.":
                        Tag = nearst_region(Format, Start, End, Genotypes[1])
                        try:
                            record = cuteSVRecord[Tag] 
                        except:
                            print("Please check whether the tag %s is in file %s." % (Tag, cuteSV_file))
                            #sys.exit(1)
                    elif genos[2] != "./.":
                        Tag = nearst_region(Format, Start, End, Genotypes[2])
                        try:
                            record = pbsvRecord[Tag]
                        except:
                            print("Please check whether the tag %s is in file %s." % (Tag, pbsv_file)) 
                            #sys.exit(1)
                    else:
                        print("Please check whether the genotypes of first two of merged SV are valid")
                        sys.exit(1)    
                    REF=record.split()[3]  
                    ALT=record.split()[4]    
                    if (svtype=='DEL' and ALT == '<DEL>') or (svtype=='INS' and ALT == '<INS>'):
                        if genos[1] != "./.":
                            Tag = nearst_region(Format, Start, End, Genotypes[1])
                            try:
                                cuteSV_record = cuteSVRecord[Tag] 
                            except:
                                print("Please check whether the tag %s is in file %s." % (Tag, cuteSV_file))
                            #sys.exit(1)  
                            if cuteSV_record.split()[4] == '<INS>' or cuteSV_record.split()[4] == '<DEL>':
                                cuteSV_record = 0
                            else:
                                cuteSV_REF = cuteSV_record.split()[3] 
                                cuteSV_ALT = cuteSV_record.split()[4]                                                                                                                                                            
                        else:
                            cuteSV_record = 0                             
                        if genos[2] != "./.":
                            Tag = nearst_region(Format, Start, End, Genotypes[2])
                            try:
                                pbsv_record = pbsvRecord[Tag] 
                            except:
                                print("Please check whether the tag %s is in file %s." % (Tag, cuteSV_file))
                            #sys.exit(1) 
                            if pbsv_record.split()[4] == '<INS>' or pbsv_record.split()[4] == '<DEL>':
                                pbsv_record = 0
                            else:
                                pbsv_REF = pbsv_record.split()[3] 
                                pbsv_ALT = pbsv_record.split()[4]                                                                                                                                
                        else:
                            pbsv_record = 0                                
                        if genos[3] != "./.":
                            Tag = nearst_region(Format, Start, End, Genotypes[3])
                            try:
                                svim_record = svimRecord[Tag] 
                            except:
                                print("Please check whether the tag %s is in file %s." % (Tag, cuteSV_file))
                            #sys.exit(1)  
                            if svim_record.split()[4] == '<INS>' or svim_record.split()[4] == '<DEL>':
                                svim_record = 0
                            else:
                                svim_REF = svim_record.split()[3] 
                                svim_ALT = svim_record.split()[4]                                                                                                                                                                                        
                        else:
                            svim_record = 0   
                        if cuteSV_record != 0:
                            REF = cuteSV_REF
                            ALT = cuteSV_ALT
                        elif pbsv_record != 0:
                            REF = pbsv_REF
                            ALT = pbsv_ALT
                        elif svim_record != 0:
                            REF = svim_REF
                            ALT = svim_ALT 
                        string_new = record.strip().split("\t")     
                        string_new[3] = REF
                        string_new[4] = ALT    
                        record = '\t'.join(string_new)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
                    out_h.write("%s\n" % record)
    elif method == "common":
        for line in merge_h:
            if line.startswith("##"):
                out_h.write(line)
            elif line.startswith("#"):
                out_h.write("%s\n" % "\t".join(line.split("\t")[:-3]))   
            else:
                lines = line.strip().split('\t') 
                Chr = lines[0]
                Start = int(lines[1])
                Infor = lines[7]
                End = int(Infor_target_values(Infor, "END"))
                Format = lines[8]
                Genotypes = lines[9:]
                genos = [parse_genotype_format(Format, g.split(";")[0], "GT") for g in Genotypes]
                zeroNum = genos.count("./.")
                nonzeroNum = len(genos) - zeroNum
                if nonzeroNum == len(genos):
                    Tag = nearst_region(Format, Start, End, Genotypes[0])
                    try:
                        record = SniffleRecord[Tag]
                    except KeyError:
                        print("Please check whether the tag %s is in file %s." % (Tag, sniffle_file))
                        sys.exit(1)
                    out_h.write("%s\n" % record)                                                                            
    elif method == "denovo":     
        for line in merge_h:
            if line.startswith("##"):
                out_h.write(line)
            elif line.startswith("#"):
                out_h.write("%s\n" % "\t".join(line.split("\t")[:-3]))   
            else:
                lines = line.strip().split('\t') 
                Chr = lines[0]
                Start = int(lines[1])
                Infor = lines[7]
                End = int(Infor_target_values(Infor, "END"))
                Format = lines[8]
                Genotypes = lines[9:]
                genos = [parse_genotype_format(Format, g.split(";")[0], "GT") for g in Genotypes]
                zeroNum = genos.count("./.")
                nonzeroNum = len(genos) - zeroNum
                if nonzeroNum == 1 and genos[0] != "./.":
                    Tag = nearst_region(Format, Start, End, Genotypes[0])
                    try:
                        record = SniffleRecord[Tag]
                    except KeyError:
                        print("Please check whether the tag %s is in file %s." % (Tag, sniffle_file))
                        sys.exit(1)
                    out_h.write("%s\n" % record)                                                
    else:
        print("Please make sure the method is 'common', 'dominant' or 'denovo'.")
        sys.exit(1)                                                   
    merge_h.close()
    out_h.close()                                                                    
                                                                                                                                                                                     
def main():
    args = parseargs()
    #args.vcf = '/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/PB_CAST_common_SV.vcf'
    #args.sniffle = '/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/Sniffles.reformat.filter.vcf'
    #args.cuteSV = '/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/cuteSV.reformat.filter.vcf'  
    #args.pbsv = '/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/pbsv.reformat.filter.vcf'
    #args.svim = '/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/svim.reformat.filter.vcf' 
    # args.outfile='/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/PB_CAST_common_original.vcf'
    merge_SV_original(args.vcf, args.sniffle, args.cuteSV, args.pbsv, args.svim, args.outfile,args.method)

if __name__ == '__main__':
    main()  
