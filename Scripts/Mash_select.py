import pandas as pd
import numpy as np
import argparse, operator, os, random, sys, time
from Bio import SeqIO

def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Genome select')
    parser.add_argument("--outputfile",help='./')
    parser.add_argument("--NGS_file",help='./screen.tab')
    parser.add_argument("--PacBio_file",help='./screen.tab')
    parser.add_argument("--t", type=float, default=0.95,help='concordant fragment size: [mu-t*dev,mu+t*dev]')    
    args = parser.parse_args()
    return args
    
def main():
    args = parseargs()     
    data1 = pd.read_csv(args.PacBio_file,sep="\t",header=None)
    data1.columns = ['Identity','Shared-hashes','median-multiplicity','p-value','Genome_fasta','Genome']
    data1['Genome'] = [x.split()[0] for x in data1['Genome']]
    selected_genome1s =  list(data1.loc[data1['Identity']>args.t,'Genome'])
    if args.NGS_file:
        data2 = pd.read_csv(args.NGS_file,sep="\t",header=None)
        data2.columns = ['Identity','Shared-hashes','median-multiplicity','p-value','Genome_fasta','Genome']
        data2['Genome'] = [x.split()[0] for x in data2['Genome']]
        selected_genome2s =  list(data2.loc[data2['Identity']>args.t,'Genome'])    
        selected_genomes = set(selected_genome1s).union(selected_genome2s)
    else:
        selected_genomes = set(selected_genome1s)           
    ViralGenome='/home1/Laisenying/Data-analysis/projects/PhageSV/ViralGenome/HGV.hq.genome.fa'
    Input = SeqIO.parse(ViralGenome,"fasta")
    Input = (record for record in Input if record.id in selected_genomes)
    SeqIO.write(Input,args.outputfile,"fasta")

if __name__=='__main__':
    main()  
