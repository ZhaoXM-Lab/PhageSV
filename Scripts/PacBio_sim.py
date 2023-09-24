import os
import sys
import re
import gzip
import time
import logging
import random
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
from multiprocessing import Pool
from distutils.spawn import find_executable
from functools import partial
from glob import glob
from shutil import rmtree
import argparse, operator, os, random, sys, time


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Pacbio Simulation')
    parser.add_argument("--genomefile",help='./genome_list.txt')
    parser.add_argument("--abund_file",help='./Abundance_abund.txt')
    parser.add_argument("--outputdir",help='./')
    parser.add_argument("--seq_depth",type=int, default=600000, help='./, about 15x') 
    args = parser.parse_args()
    return args


def tidy_taxon_names(x):
    """
    Remove special characters from taxon names
    """
    x = re.sub(r'[()\/:;, ]+', '_', x)
    return x
    
def _genome_size(x):
    """
    Get the total bp for the genome sequence of all genomes        
    Parameters
    ----------
    x : pd.Series
       Series that includes 'Fasta' in the index
    """
    bp = 0
    for record in SeqIO.parse(x['Fasta'], 'fasta'):
        bp += len(record.seq)
    return bp

def _write_reads(outFH, read_files, read_type):
    for in_file in read_files:
        taxon = os.path.split(os.path.split(in_file)[0])[1]
        for i,record in enumerate(SeqIO.parse(in_file, read_type)):
            # renaming fastq read
            name =  '{}__SEQ{}'.format(taxon, i)
            record.id = name
            record.description = name
            SeqIO.write(record, outFH, read_type)
        # delete temporary file
        os.remove(in_file)
     
def _combine_reads(read_files, output_dir, output_file, read_type='fastq',
                   gzip_out=False):
    """
    Combine fastq read files into 1 read file.
    Parameters
    ----------
    read_files : list
        All read files to combine
    output_dir : str
        Output directory path
    read_type : str
        Read file type
    gzip_out : bool
        gzip output?
    """
    if read_type == 'fastq':
        suffix = '.fq'
    elif read_type == 'fasta':
        suffix = '.fa'
    else:
        msg = 'Cannot determine suffix for read type: {}'
        raise ValueError(msg.format(read_type))
    output_file = os.path.join(output_dir, output_file + suffix)
    if gzip_out is True:
        output_file += '.gz'
        with gzip.open(output_file, 'wt') as outFH:
            _write_reads(outFH, read_files, read_type)
    else:
        with open(output_file, 'w') as outFH:
            _write_reads(outFH, read_files, read_type)
                                       
def combine_reads_by_sample(sample, fq_files, temp_dir, file_prefix, output_dir,
                            read_type='fastq', gzip_out=False, debug=False):
    """
    Concat all sample-taxon read files into per-sample read files
    Parameters
    ----------
    sample : str
        Sample ID
    temp_dir : str
        Temporary directory path
    file_prefix : str
        Output file prefix
    output_dir : str
        Output directory path
    debug : bool
        gzip output?
    debug : bool
        Debug mode
    """
    # sorting fastq files by read pair
    sample = str(sample)
    R1_files = [x[1] for x in fq_files if x[0] == sample]
    try:
        R2_files = [x[2] for x in fq_files if x[0] == sample]
    except IndexError:
        R2_files = None
    # writing files
    output_dir = os.path.join(output_dir, str(sample))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    ## read1
    R1_files = _combine_reads(R1_files, output_dir, 'R1',
                              read_type=read_type, gzip_out=gzip_out)
    ## read2
    if R2_files is not None:
        R2_files = _combine_reads(R2_files, output_dir, 'R2',
                                  read_type=read_type, gzip_out=gzip_out)    
    return [R1_files, R2_files]

def load_genome_table(in_file, nproc=1):
    """
    Load genome table
    Parameters
    ----------
    in_file : str
        input file path
    """
    nproc = int(nproc)
    df = pd.read_csv(in_file, sep='\t')    
    ## check headers
    diff = set(['Taxon','Fasta']) - set(df.columns.values)
    if len(diff) > 0:
        diff = ','.join(diff)
        raise ValueError('Cannot find table columns: {}'.format(diff))
    # getting genome sizes
    logging.info('  Getting genome sizes (no. of threads: {})'.format(nproc))
    if nproc > 1:
        p = Pool(nproc)
        df['Genome_size'] = p.map(_genome_size, [x for i,x in df.iterrows()])
        p.close()
    else:
        df['Genome_size'] = [_genome_size(x) for i,x in df.iterrows()]
    # tidy taxon names
    df['Taxon'] = df['Taxon'].astype(str).apply(tidy_taxon_names)    
    return df

def load_abund_table(in_file):
    """
    Load abundance table
    Parameters
    ----------
    in_file : str
        input file path
    """
    df = pd.read_csv(in_file, sep='\t')    
    ## check headers
    diff = set(['Community','Taxon','Perc_rel_abund']) - set(df.columns.values)
    if len(diff) > 0:
        diff = ','.join(diff)
        raise ValueError('Cannot find table columns: {}'.format(diff))
    # tidy taxon names
    df['Taxon'] = df['Taxon'].astype(str).apply(tidy_taxon_names)
    return df   
    
def sample_taxon_list(genome_table, abund_table):
    """
    Create [sample,taxon] list of lists
    Parameters
    ----------
    genome_table : pd.DataFrame
    abund_table : pd.DataFrame
    """
    # joining tables
    df = abund_table.merge(genome_table, on=['Taxon'])
    df['Perc_rel_abund'] = df['Perc_rel_abund']*df['Genome_size']/np.mean(df['Genome_size'])
    df['Perc_rel_abund'] = df['Perc_rel_abund']*100/np.sum(df['Perc_rel_abund'])
    # convert to a list of lists
    sample_taxon = []
    cols = ['Community', 'Taxon', 'Genome_size', 'Fasta', 'Perc_rel_abund']
    for i,x in df[cols].iterrows():
        sample_taxon.append(x.tolist())
    return sample_taxon
    
                 
def sim_pacbio(sample_taxon,output_dir, seq_depth, sl_params, temp_dir,
               nproc=1, rndSeed=None, gzip_out=False, debug=False):
    """
    Simulate pacbio reads
    Parameters
    ----------
    sample_taxon : list
        [Community,Taxon,Genome_size,Fasta,Perc_rel_abund,Fold]
    output_dir : str
        Output director for all read files
    seq_depth : int
        Sequencing depth per sample
    sl_params : dict
        Parameters provided to SimLord
    temp_dir : str
        Temporary file directory
    nproc : int
        Number of parallel processes
    debug : bool
        Debug mode
    """
    # setting seed (also passed to read simulator)
    if rndSeed is not None:
        random.seed(rndSeed)        
    # check that simulator exists
    exe = 'simlord'
    if find_executable(exe) == '':
        raise IOError('Cannot find executable: {}'.format(exe))    
    # calculating fraction of total reads to simulate
    for x in sample_taxon:
        perc_rel_abund = float(x[4])
        genome_size = int(x[2])
        x.append(int(perc_rel_abund/100.0 * seq_depth))
    # directories
    ## temp dir
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
    ## output dir
    output_dir = os.path.join(output_dir, 'pacbio')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)        
    # simulate per sample
    logging.info('Simulating PacBio reads...')
    func = partial(sim_simlord,
                   sl_params=sl_params,
                   temp_dir=temp_dir,
                   debug=debug)
    if debug is True:
        fq_files = map(func, sample_taxon)
    else:
        p = Pool(nproc)
        fq_files = p.map(func, sample_taxon)
    fq_files = list(fq_files)
    # combining all reads by sample
    logging.info('Combining simulated reads by sample...')    
    comms = list(set([x[0] for x in sample_taxon]))
    func = partial(combine_reads_by_sample,
                   fq_files=fq_files,
                   temp_dir=temp_dir,
                   file_prefix='pacbio',
                   output_dir=output_dir,
                   gzip_out=gzip_out,
                   debug=debug)
    if debug is True:
        res = map(func, comms)
    else:
        p = Pool(nproc)
        res = p.map(func, comms)
    res = list(res)
    # removing temp dir
    logging.info('Removing temp directory...')
    rmtree(temp_dir)    
    # status
    for sample_list in res:
        for file_name in sample_list:
            if file_name is not None:
                logging.info('File written: {}'.format(file_name))                                           

def sim_simlord(x, sl_params, temp_dir, debug=False):
    """
    Simulate pacbio reads with simlord
    Parameters
    ----------
    x : list
        [Community,Taxon,Genome_size,Fasta,Perc_rel_abund,Fold]
    """
    community = str(x[0])
    taxon = str(x[1])
    fasta = str(x[3])
    perc_rel_abund = float(x[4])
    n_reads = float(x[5])    
    # output
    ## temporary directories
    out_dir = os.path.join(temp_dir, str(community), taxon)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    ## temporary files
    output_prefix = os.path.join(out_dir, 'pacbio')    
    # art command
    sl_params = ' '.join(['{} {}'.format(k,v) for k,v in sl_params.items()])
    cmd = 'simlord {sl_params} --num-reads {n_reads}'
    cmd += ' --read-reference {input} {output_prefix}'
    cmd = cmd.format(sl_params=sl_params,
                     n_reads=int(n_reads),
                     input=fasta,
                     output_prefix=output_prefix)    
    ## system call
    if debug is True:
        sys.stderr.write('CMD: ' + cmd + '\n')
    try:
        res = subprocess.run(cmd, check=True, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise e
    if debug is True:
        sys.stderr.write(res.stderr.decode() + '\n')
        sys.stderr.write(res.stdout.decode() + '\n')
        
    # check that files have been created
    R0_file = output_prefix + '.fastq'
    if os.path.isfile(R0_file):
        return [community, R0_file]
    else:
        msg = 'Cannot find simlord output fastq file!'
        raise ValueError(msg)        

def main():
    args=parseargs()
    genome_file=args.genomefile
    abund_file=args.abund_file
    output_dir=args.outputdir
    seq_depth=args.seq_depth
    sl_params = {'--min-readlength' : 500,
        '--sample-readlength-from-fastq':'A01.subreads.fastq'}
    temp_dir=output_dir+'/PacBio_tmp'
    nproc = 5        
    genome_table = load_genome_table(genome_file,nproc=2)
    abund_table = load_abund_table(abund_file)
    sample_taxon = sample_taxon_list(genome_table,abund_table)
    sim_pacbio(sample_taxon,
                            output_dir=output_dir,
                            seq_depth=seq_depth,
                            sl_params=sl_params,
                            temp_dir=temp_dir,
                            nproc=nproc)

if __name__ == '__main__':
    main()
