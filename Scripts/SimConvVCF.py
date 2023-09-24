# Convert Sim-it into True VCF file into 
import pandas as pd
import numpy as np
import argparse, operator, os, random, sys, time
from Bio import SeqIO
import random
import os

def parseargs():    # handle user arguments
    # args.file='/home1/Laisenying/Data-analysis/projects/PhageSV/Test/pbsv/Mash_select/screen.tab'
    parser = argparse.ArgumentParser(description='Genome select')
    parser.add_argument("--inputfile",help='/genome_select.fasta')  
    parser.add_argument("--outputfile",help='./') 
    parser.add_argument("--reference",help='./')
    args = parser.parse_args()
    return args

def main(): 
    args=parseargs()       
    file=open(args.inputfile,"r")
    Genomefasta=args.reference
    Input = SeqIO.parse(Genomefasta,"fasta")
    strings = []
    strings.append('##fileformat=VCFv4.2')
    strings.append('##FILTER=<ID=PASS,Description="All filters passed">')
    strings.append('##fileDate=2022-07-07T17:36:21.790Z')
    for record in Input:
        Length = len(record.seq)
        contig=record.id
        strings.append("##contig=<ID="+contig+',length='+str(Length))
    strings.append("##ALT=<ID=INS,Description=\"Insertion of novel sequence relative to the reference\">") 
    strings.append("##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">")  
    strings.append("##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the reference\">")  
    strings.append("##ALT=<ID=INV,Description=\"Inversion of reference sequence\">")    
    strings.append("##ALT=<ID=BND,Description=\"Breakend of translocation\">")
    strings.append("##ALT=<ID=IDUP,Description=\"Region of invered elevated copy number relative to the reference\">")      
    strings.append("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">")
    strings.append("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">")
    strings.append("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">")    
    strings.append('\t'.join(['#CHROM',"POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE"]))                    
    outputfile=open(args.outputfile,"w")
    for string in strings:
        outputfile.write(string+"\n")  
    i=1
    for line in file.readlines():
        if line[0] == '#':
            continue  
        else:
            contig=line.split()[0]
            Pos=int(line.split()[1])
            if Pos == 0:
                continue
            SVID='SimSV.'+str(i)      
            i+=1
            Ref_Base='N'
            ALT_Base="<"+line.split()[3]+">"   
            QUAL='.'
            FILTER='PASS'
            if line.split()[3] == 'DEL':
                SVLEN = -int(line.split()[2])
                END = Pos - SVLEN
                Info=';'.join(['SVTYPE='+line.split()[3],'END='+str(END),'SVLEN='+str(SVLEN)])    
            elif line.split()[3] == 'INS':                     
                SVLEN = int(line.split()[2])  
                END = Pos   
                Info=';'.join(['SVTYPE='+line.split()[3],'END='+str(END),'SVLEN='+str(SVLEN)])    
            elif line.split()[3] == 'DUP': 
                SVLEN = int(line.split()[2].split("x")[0])*int(line.split()[2].split("x")[1])   
                END = Pos + SVLEN
                Info=';'.join(['SVTYPE='+line.split()[3],'END='+str(END),'SVLEN='+str(SVLEN)])    
            elif line.split()[3] == 'INV': 
                SVLEN = int(line.split()[2]) 
                END = Pos + SVLEN
                Info=';'.join(['SVTYPE='+line.split()[3],'END='+str(END)])    
            else:
                SVLEN = int(line.split()[2])
                END = Pos + SVLEN
                Info=';'.join(['SVTYPE='+line.split()[3],'END='+str(END),'SVLEN='+str(SVLEN)])                                                                                                        
            new_line = '\t'.join([contig,str(Pos),SVID,Ref_Base,ALT_Base,QUAL,FILTER,Info,".","."])        
            outputfile.write(new_line+"\n")   
    outputfile.close()  
        
if __name__ == '__main__':
    main()  
