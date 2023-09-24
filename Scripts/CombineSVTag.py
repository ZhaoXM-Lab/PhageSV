import argparse
import sys
import math
import collections
import operator
#from tinyfasta import FastaParser

def parseargs():    # handle user arguments
    # args.file='/home1/Laisenying/Data-analysis/projects/PhageSV/Test/pbsv/Mash_select/screen.tab'
    parser = argparse.ArgumentParser(description='Filt SV file based on the number of support records.')
    parser.add_argument("--vcf",type=str,default='',help='vcf file')  
    parser.add_argument("--outfile",help='path to the output vcf')             
    args = parser.parse_args()
    return args
    
def SV_tag(lines):
    ID = lines[0]
    Chr, Start = lines[:2]
    Infor = lines[7]
    ### change the ID based on the Infor
    Chr2 = Infor_target_values(Infor, "CHR2")
    End = Infor_target_values(Infor, "END")
    SVType = Infor_target_values(Infor, "SVTYPE")
    if "SVLEN" in Infor:
        SVLength = Infor_target_values(Infor, "SVLEN")
    elif "AVGLEN" in Infor:
        SVLength = Infor_target_values(Infor, "AVGLEN")        
    if SVLength.startswith("-"):
        SVLength = SVLength.lstrip("-")
    if SVType == "BND":
        SVType = "TRA"
    if SVType == "TRA":
        SVLength = 0
    Tag = "%s~%s--%s~%s--%s--%s" % (Chr, Start, Chr2, End, SVLength, SVType) ### just one chromosomes    
    return Tag
    
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

def Format_tag(Format, Geno):
    """
    Format: GT:PSV:LN:DR:ST:TY:CO
    geno: 1/1:NA:137:1,36:+-:INS,INS:1_136673-1_136891,1_136999-1_137209

    geno:
    1/1:NA:90:0,17:+-:INS:1_137053-1_137183;0/1:NA:188:29,19:+-:INS:1_136930-1_137112

    ['NaN-0-NaN', '1_10169-X_449438-0-TRA']
    """
    genos = Geno.split(";")
    Tags = []
    for geno in genos:
        LN = parse_genotype_format(Format, geno, "LN")
        TY = parse_genotype_format(Format, geno, "TY")
        CO = parse_genotype_format(Format, geno, "CO")
        ### set length of translocation as 0
        if TY == "TRA":
            LN = "0"
        tag = "%s--%s--%s" % (CO, LN, TY)
        Tags.append(tag)
    TagLen = len(Tags)
    if TagLen == 1:
        return Tags[0]
    else:
        return ";".join(Tags)
                
def vcf_position_sample(vcf_file,out_file):
    TagSampleTag = collections.defaultdict(dict)
    SampleTages = collections.defaultdict(dict)
    out_h = open(out_file,"w")
    vcf_h = open(vcf_file,"r")
    for line in vcf_h:
        line = line.strip()
        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            lines = line.split("\t")
            sampleNames = lines[9:]
            out_h.write("Tag\t%s\n" % "\t".join(sampleNames))
        else:
            lines = line.strip("\n").split("\t")
            Infor = lines[7]
            Format = lines[8]
            Genos = lines[9:]
            SVType = Infor_target_values(Infor, "SVTYPE")
            # GenoTags = [Format_tag(Format, geno) for geno in Genos]
            GenoTags = [Format_tag(Format, geno) for geno in Genos]
            ### change "NaN-0-NaN" to "-"
            GenoTags = ["-" if g == "NaN--0--NaN" else g for g in GenoTags]
            Tag = SV_tag(lines)                                                    
            out_h.write("%s\t%s\n" % (Tag, "\t".join(GenoTags)))                                  
    out_h.close()
        
def main():
    args = parseargs()    
    vcf_position_sample(args.vcf,args.outfile)
    
if __name__ == "__main__":
    main()  
