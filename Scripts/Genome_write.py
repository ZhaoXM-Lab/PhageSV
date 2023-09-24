import pandas as pd
import numpy as np
import argparse, operator, os, random, sys, time
from Bio import SeqIO
import random
import os

def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Genome write')
    parser.add_argument("--SimSVdir",help='./')  
    args = parser.parse_args()
    return args

def main():  
    args = parseargs()       
    Input = SeqIO.parse(args.SimSVdir+"/Genome_withSV.fasta",'fasta')
    os.system("mkdir -p "+ args.SimSVdir + "/Genomes")
    for record in Input:
        file=open(args.SimSVdir + "/Genomes/"+record.id+'.fa',"w")
        file.write(">"+record.id+"\n") 
        file.write(str(record.seq)+"\n")
        file.close()

    path_data={"Taxon":[],"Fasta":[]}        
    Input = SeqIO.parse(args.SimSVdir+"/Genome_withSV.fasta",'fasta')
    for record in Input:
        path_data['Taxon'].append(record.id)
        path_data["Fasta"].append(args.SimSVdir + "/Genomes/"+record.id+'.fa')     
    path_data = pd.DataFrame(path_data)
    path_data.to_csv(args.SimSVdir+"/genome_list.txt",sep="\t")    
    
if __name__=='__main__':
    main()
