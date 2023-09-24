import pandas as pd
import numpy as np
import argparse, operator, os, random, sys, time
from Bio import SeqIO
from collections import defaultdict

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

def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Extracting SV information from VCF file')
    parser.add_argument("--VCF_file",help='./.vcf')
        parser.add_argument("--VCF_list",help='./vcf_file_path')
    parser.add_argument("--outputfile",help='./')
    args = parser.parse_args()
    return args
    
def main():
    args = parseargs()
    vcf_list = pd.read_csv(args.VCF_list,header=None)
    vcf_list['sample'] = [x.split("/")[-1].split("_")[0] for x in vcf_list[0]]
    vcf_list.columns = ['vcf_path','sample']

    vcf_file_path={}
    for x,y in zip(vcf_list['vcf_path'],vcf_list['sample']):
        vcf_file_path[y] = x
    
    vcfrecord_dict = {}
    for sample,vcf in vcf_file_path.items():
        vcfrecord =  tag_record(vcf)
        vcfrecord_dict[sample] = vcfrecord

    SV_inf = defaultdict(list)
    vcf_file = open(args.VCF_file,"r")
    for line in vcf_file.readlines():
        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            sampleNames = np.array(line.split()[9:])
            continue
        else:
            lines = line.strip('\t').split()
            SV_inf['contig'].append(lines[0])
            SV_inf['Start'].append(int(lines[1]))
            SV_inf['SVID'].append(lines[2])
            Info = lines[7]
            Info_list={}
            for x in Info.split(';'):
                Info_list[x.split('=')[0]] = x.split('=')[1]
            SV_inf['SVTYPE'].append(Info_list["SVTYPE"])
            SV_inf['SVLEN'].append(int(Info_list["AVGLEN"]))
            SV_inf['End'].append(int(Info_list["END"]))
            SV_inf['SUPP_num'].append(int(Info_list["SUPP"]))
            sample_genos = lines[9:]
            GenoTags = np.array(["-" if g == "./.:NaN:0:0,0:--:NaN:NaN" else g for g in sample_genos])
            samples = sampleNames[GenoTags != '-']
            genos =  GenoTags[GenoTags != '-']
            Format = lines[8]
            SV_inf['Tag'].append("%s~%s--%s~%s--%s--%s" % (lines[0], lines[1],lines[0],Info_list["END"],Info_list["AVGLEN"].split("-")[-1], Info_list["SVTYPE"]))
            if Info_list['SVTYPE'] != 'INS':
                SV_inf['REF'].append(lines[3])
                SV_inf['ALT'].append(lines[4])
                SV_inf['SUPP_samples'].append(";".join(list(samples)))
            else:
                if len(samples) == 1:
                    SV_inf['REF'].append(lines[3])
                    SV_inf['ALT'].append(lines[4])
                    SV_inf['SUPP_samples'].append(";".join(list(samples)))
                else:
                    SV_inf['SUPP_samples'].append(";".join(list(samples)))
                    refs=[]
                    alts = []
                    for sample,geno in zip(samples,genos):
                        sample_vcf_file = vcfrecord_dict[sample]
                        Tag2 = parse_genotype_format(Format, geno, "CO")
                        Length = parse_genotype_format(Format, geno, "LN")
                        Type = parse_genotype_format(Format, geno, "TY")
                        Tag = "%s--%s--%s" % (Tag2, Length, Type)
                        record = sample_vcf_file[Tag]
                        REF = record.split("\t")[3]
                        ALT = record.split("\t")[4]
                        refs.append(REF)
                        alts.append(ALT)
                    SV_inf['REF'].append(";".join(refs))
                    SV_inf['ALT'].append(";".join(alts))
    SV_bed_inf = pd.DataFrame(SV_inf)
    SV_bed_inf = SV_bed_inf.loc[:,['contig', 'Start', 'End', 'SVID', 'SVTYPE', 'SVLEN', 'SUPP_num', 'Tag','REF', 'ALT', 'SUPP_samples']]
    SV_bed_inf.to_csv(args.outputfile,sep="\t")

if __name__=='__main__':
    main()
