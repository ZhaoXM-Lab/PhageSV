from collections import defaultdict, Counter
import multiprocessing as mp
import networkx as nx
import numpy as np
import argparse
import time
import re

def parseargs():    # handle user arguments
    # args.file='/home1/Laisenying/Data-analysis/projects/PhageSV/Test/pbsv/Mash_select/screen.tab'
    parser = argparse.ArgumentParser(description='Filt SV file based on the number of support records.')
    parser.add_argument("--vcf_list",type=str,default='',help='list of vcf (separated by commas)')  
    parser.add_argument("--outfile",help='path to the output vcf')     
    parser.add_argument("--allele_freq",default=0.2, type=float,help='use allele frequency to determine genotype (default: 0.2)')
    parser.add_argument("--interval",default=0, type=int,help='segregated all SVs into non-overlapping blocks, interval of block (default: 0)')   
    parser.add_argument("--extended_dis",default=100, type=int,help='extended distance between the <TRA> start and end upstream and downstream (default: 100)')  
    parser.add_argument("--overlap_perc",default=0.5, type=float,help='except for <TRA>, overlap percentage of the two SVs (default: 0.5)')  
    parser.add_argument("--workers",default=5, type=int,help='number of CPUs to use (default: 5)')           
    args = parser.parse_args()
    return args
    
class SvInfo(object):
    def __init__(self):
        self.patient = ''
        self.sv_id = ''
        self.chr1 = ''
        self.chr2 = ''
        self.start = 0
        self.end = 0
        self.svlen = 0
        self.strands = ''
        self.sv_type = ''
        self.genotype = ''
        self.dr = 0
        self.dv = 0    
    
def read_vcf(vcf_list,af):
    contigs = set()
    vcf_file_list = []
    for line in open(vcf_list,"r").readlines():
        vcf_file_list.append(line.strip("\n"))        
    all_sv_form = {'<DEL>': defaultdict(list),
                   '<INS>': defaultdict(list),
                   '<INV>': defaultdict(list),
                   '<DUP>': defaultdict(list),
                   '<INVDUP>': defaultdict(list)}                 
    for per_file in vcf_file_list:        
        with open(per_file,'r') as f:
            for line_raw in f.readlines():
                if not line_raw.startswith('#'):
                    line = line_raw.strip('\n').split('\t')
                    info_list = line[7].split(';')                                          
                    info_dict = {}
                    for x in info_list:
                        x = x.split('=')
                        if len(x) > 1:
                            info_dict[x[0]] = x[1]
                    temp_sv = SvInfo()
                    temp_sv.patient = per_file
                    temp_sv.sv_id = line[2]
                    temp_sv.chr1 = line[0]
                    temp_sv.chr2 = info_dict['CHR2']
                    temp_sv.start = int(line[1])
                    temp_sv.end = int(info_dict['END'])
                    temp_sv.svlen = int(info_dict['SVLEN'])
                    temp_sv.strands = info_dict['STRANDS']
                    temp_sv.sv_type = info_dict['SVTYPE']
                    dr = int(line[-1].split(':')[1])
                    dv = int(line[-1].split(':')[2])   
                    # custom genotype
                    if dv/(dr + dv) <= af:
                        temp_sv.genotype = '0/0'
                    elif af < dv/(dr + dv) <= 1 - af:
                        temp_sv.genotype = '0/1'
                    else:
                        temp_sv.genotype = '1/1'  
                    temp_sv.dr = dr
                    temp_sv.dv = dv
                    all_sv_form['<'+info_dict['SVTYPE']+'>'][line[0]].append(temp_sv)
                re_pattern_len = re.compile(r'length=(\d+)>')  
                if line_raw.startswith('##contig'):
                    contig_len = int(re.findall(pattern=re_pattern_len, string=line_raw)[0])
                    contigs.add((contig_len, line_raw)) 
    for type_sv in all_sv_form:
        for chrom in all_sv_form[type_sv]: 
            start = np.array([i.start for i in all_sv_form[type_sv][chrom]]) 
            # end = start + svlen
            end = np.array([i.start + abs(i.svlen) for i in all_sv_form[type_sv][chrom]]) 
            order = np.lexsort((end, start))   
            per_chrom = np.array(all_sv_form[type_sv][chrom])[order]
            all_sv_form[type_sv][chrom] = list(per_chrom) 
    return vcf_file_list, contigs, all_sv_form     
    
# find non-overlapping groups
def break2block(all_sv_form, interval):
    all_sv_block = {'<DEL>': defaultdict(list),
                    '<INS>': defaultdict(list),
                    '<INV>': defaultdict(list),
                    '<DUP>': defaultdict(list),
                    '<INVDUP>': defaultdict(list)}
    for type_sv in all_sv_form:
        if type_sv == '<TRA>':  # TRA don't have to divide
            for chrom in all_sv_form[type_sv]:
                all_sv_block[type_sv][chrom].append(all_sv_form[type_sv][chrom])
        else:
            for chrom in all_sv_form[type_sv]:
                start = np.array([i.start for i in all_sv_form[type_sv][chrom]])
                # end = start + svlen
                end = np.array([i.start + abs(i.svlen) for i in all_sv_form[type_sv][chrom]])
                break_idx = []
                if len(start) != 1:
                    for s, e in zip(range(1, len(start)), range(len(end) - 1)):
                        if start[s] - end[e] > interval:
                            break_idx.append(s)
                    if len(break_idx) == 0:
                        print('The %s on the %s chromosome belongs to a block!' % (type_sv, chrom))
                else:
                    break_idx = []
                    print('%s chromosome have only one %s!' % (chrom, type_sv))
                block = []
                if len(break_idx) != 0:
                    block.append(all_sv_form[type_sv][chrom][0:break_idx[0]])
                    if len(break_idx) > 1:
                        for b in range(len(break_idx) - 1):
                            block.append(all_sv_form[type_sv][chrom][break_idx[b]:break_idx[b + 1]])
                    block.append(all_sv_form[type_sv][chrom][break_idx[-1]:len(all_sv_form[type_sv][chrom])])
                else:
                    block = [all_sv_form[type_sv][chrom]]
                all_sv_block[type_sv][chrom] = block
    return all_sv_block   
    
# sub-process
def pro_block2clique2merge(all_sv_block, type_sv, chrom, overlap_perc):
    per_chrom_merge = []
    for per_block in all_sv_block[type_sv][chrom]:
        g_non_tra = nx.Graph()
        for per_sv in per_block:
            g_non_tra.add_node(per_sv)
            for other_sv in per_block:
                start1 = per_sv.start
                end1 = per_sv.start + abs(per_sv.svlen)  # end = start + svlen
                start2 = other_sv.start
                end2 = other_sv.start + abs(other_sv.svlen)  # end = start + svlen
                overlap = overlap_fun(start1, end1, start2, end2)
                if overlap != 0:
                    rate1 = overlap / (end1 - start1)
                    rate2 = overlap / (end2 - start2)
                    if rate1 > overlap_perc and rate2 > overlap_perc:
                        g_non_tra.add_edge(per_sv, other_sv, weight=rate1 + rate2)
        # to find the largest clique
        max_clique_list = []
        temp_nodes = set(g_non_tra.nodes)
        while temp_nodes:
            temp_graph = nx.Graph.subgraph(g_non_tra, list(temp_nodes))
            clique_list = [(i, len(i)) for i in nx.algorithms.clique.find_cliques(temp_graph)]
            clique_list = sorted(clique_list, key=lambda x: x[1], reverse=True)
            # take largest clique
            clique_max_num = clique_list[0][1]
            clique_max_list = [i for i in clique_list if i[1] == clique_max_num]
            # sort the elements in the group
            clique_max_list_2 = []
            for per in clique_max_list:
                per = list(per)
                per[0] = sorted(per[0], key=lambda x: x.patient + x.sv_id)
                per = tuple(per)
                clique_max_list_2.append(per)
            # sort group
            clique_max_list = sorted(clique_max_list_2, key=for_list_sort)
            # take the first largest
            clique_max = clique_max_list[0][0]
            max_clique_list.append(tuple(clique_max))
            temp_nodes -= set(clique_max)
        # merge clique
        final_merge = []
        if len(max_clique_list) > 1:
            last_clique_list, stab_clique_list = [], []
            stab_clique_list.extend(max_clique_list)
            while last_clique_list != stab_clique_list:
                last_clique_list = []
                last_clique_list.extend(stab_clique_list)
                stab_clique_list = []
                stab_clique_list.extend(find_stab_clique(last_clique_list, g_non_tra))
            # a group containing only one element
            clique_list_one = [i[0] for i in stab_clique_list if len(i) == 1]
            clique_list_one = sorted(clique_list_one, key=lambda x: x.patient + x.sv_id)
            # the elements in the group are greater than two
            clique_list_two = [i for i in stab_clique_list if len(i) >= 2]
            # sort the elements in the group
            clique_list_two_2 = []
            for per in clique_list_two:
                per = list(per)
                per = sorted(per, key=lambda x: x.patient + x.sv_id)
                per = tuple(per)
                clique_list_two_2.append(per)
            # sort group
            clique_list_two_2 = sorted(clique_list_two_2, key=lambda x: x[0].patient + x[0].sv_id)
            # add the SV to the nearest group
            for one_sv in clique_list_one:
                max_weight, max_index = 0, 'none'
                for idx_cl, per_cli in enumerate(clique_list_two_2):
                    num_weight = 0
                    for sv in per_cli:
                        if g_non_tra.get_edge_data(sv, one_sv):
                            num_weight += g_non_tra.get_edge_data(sv, one_sv)['weight']
                    if num_weight > max_weight:
                        max_index, max_weight = idx_cl, num_weight
                if max_index == 'none':
                    clique_list_two_2.append(tuple([one_sv]))
                else:
                    new_per_cli = list(clique_list_two_2[max_index])
                    new_per_cli.append(one_sv)
                    clique_list_two_2[max_index] = tuple(new_per_cli)
            final_merge.extend(clique_list_two_2)
            per_chrom_merge.append(final_merge)
        else:
            final_merge.append(max_clique_list[0])
            per_chrom_merge.append(final_merge)
    return per_chrom_merge, chrom    
    
# host process
def block2clique2merge(all_sv_block, extended_dis, overlap_perc, worker):
    all_sv_merge = {'<DEL>': defaultdict(list),
                    '<INS>': defaultdict(list),
                    '<INV>': defaultdict(list),
                    '<DUP>': defaultdict(list),
                    '<INVDUP>': defaultdict(list)}        
    for type_sv in all_sv_block:
        if type_sv != '<TRA>':
            print('Start %s: block2clique2merge' % type_sv)
            ppp = mp.Pool(worker)
            multi_res = [ppp.apply_async(func=pro_block2clique2merge,
                                         args=(all_sv_block, type_sv, chrom, overlap_perc))
                         for chrom in all_sv_block[type_sv]]
            ppp.close()
            ppp.join()
            for res in multi_res:
                all_sv_merge[type_sv][res.get()[1]] = res.get()[0]
    return all_sv_merge   
    
def overlap_fun(start1, end1, start2, end2):
    return max(max((end2 - start1), 0) - max((end2 - end1), 0) - max((start2 - start1), 0), 0)

def find_stab_clique(clique_list_raw, sv_graph):
    clique_list = []
    clique_list.extend(clique_list_raw)
    gc = nx.Graph()
    for per1 in clique_list:
        for per2 in clique_list:
            if per1 != per2:
                num = 0
                for node1 in per1:
                    for node2 in per2:
                        if sv_graph.get_edge_data(node1, node2):
                            num += 1
                if num >= 2:  # if two cliques have two SV edges in the original SV graph
                    gc.add_edge(per1, per2)
    clusters = list(nx.connected_components(gc))
    for cis in clusters:
        combine = set()
        for ci in cis:
            clique_list.remove(ci)
            combine |= set(ci)
        clique_list.append(tuple(combine))
    clique_list = sorted(clique_list, key=lambda x: x[0].patient + x[0].sv_id)
    return clique_list

# sort to make results reproducible
def for_list_sort(per):
    per_list = []
    for i in per[0]:
        per_list.append(i.patient + i.sv_id)
    return tuple(per_list)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                

def output_vcf(all_sv_merge, vcf_file_list, contigs, output):
    f_output = open(output, 'w')
    # annotation information for vcf
    f_output.write('##fileformat=VCFv4.1\n')
    f_output.write('##source=cliqueSV_maxflow\n')
    f_output.write('##fileDate=%s\n' % time.strftime("%Y%m%d", time.localtime()))
    f_output.write('##ALT=<ID=INV,Description="Inversion">\n')
    f_output.write('##ALT=<ID=DUP,Description="Duplication">\n')
    f_output.write('##ALT=<ID=INV,Description="Inversion of reference sequence">\n')
    f_output.write('##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">\n')
    f_output.write('##ALT=<ID=DEL,Description="Deletion relative to the reference">\n')
    f_output.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    f_output.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
    f_output.write('##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Number of detected tools">\n')
    f_output.write('##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of detected tools">\n')
    f_output.write('##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Vector of detected tools">\n')
    f_output.write('##INFO=<ID=AVGLEN,Number=1,Type=Integer,Description="Average length of Merged SVs">\n')
    f_output.write('##INFO=<ID=CIEND,Number=1,Type=String,Description="# CIEND">\n')     
    f_output.write('##INFO=<ID=SVMETHOD,Number=1,Type=Integer,Description="Methods for merging SVs">\n')
    f_output.write('##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">\n')
    f_output.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')  
    f_output.write('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">\n')  
    f_output.write('##INFO=<ID=CILEN,Number=2,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends">\n')      
    f_output.write('##INFO=<ID=RE,Number=1,Type=Integer,Description="Number of read support this record">\n')
    f_output.write('##INFO=<ID=STRANDS,Number=A,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">\n') 
    f_output.write('##INFO=<ID=STRAND,Number=A,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">\n') 
    f_output.write('##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names of SVs (comma separated)">\n') 
    f_output.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency.">\n')
    f_output.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n') 
    f_output.write('##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# High-quality reference reads">\n')
    f_output.write('##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# High-quality variant reads">\n') 
    f_output.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="# Phred-scaled genotype likelihoods rounded to the closest integer">\n')
    f_output.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="# Genotype quality">\n')  
    f_output.write('##FORMAT=<ID=PSV,Number=1,Type=String,Description="# PSV">\n')  
    f_output.write('##FORMAT=<ID=LN,Number=1,Type=Integer,Description="# SV length">\n')
    f_output.write('##FORMAT=<ID=ST,Number=1,Type=String,Description="# Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">\n')     
    f_output.write('##FORMAT=<ID=TY,Number=1,Type=String,Description="# SV type">\n') 
    f_output.write('##FORMAT=<ID=CO,Number=1,Type=String,Description="# Contig name">\n')     
    contigs = sorted(contigs, key=lambda x: x[0], reverse=True)
    ###sample names
    vcf_file_names = [n.split("/")[-1].split("_")[0].split(".")[0] for n in vcf_file_list]
    print(vcf_file_names)
    for con in contigs:
        f_output.write(con[1])
    f_output.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +'\t'.join(vcf_file_names) + '\n') #vcf_file_list
    total_wait_dict = defaultdict(list)
    for type_sv in all_sv_merge:
        type_sv_name = type_sv[1:-1]
        for chrom in all_sv_merge[type_sv]:
            for per_block in all_sv_merge[type_sv][chrom]:
                for per_clique in per_block:
                    per_clique_list = list(per_clique)
                    per_clique_list = sorted(per_clique_list, key=lambda x: x.patient + x.sv_id)
                    per_clique_start = [i.start for i in per_clique_list]
                    per_clique_end = [i.end for i in per_clique_list]
                    # take the most frequent
                    max_common_start_idx = per_clique_start.index(Counter(per_clique_start).most_common(1)[0][0])
                    max_read = per_clique_list[max_common_start_idx]
                    start = max_read.start
                    end = max_read.end
                    strands_max = max_read.strands
                    chr1 = max_read.chr1
                    chr2 = max_read.chr2
                    cipos = [0, 0]
                    max_start = int(np.max(np.array(per_clique_start)))
                    min_start = int(np.min(np.array(per_clique_start)))
                    cipos[0] = min_start - start
                    cipos[1] = max_start - start
                    cipos = ','.join([str(i) for i in cipos])
                    ciend = [0, 0]
                    max_end = int(np.max(np.array(per_clique_end)))
                    min_end = int(np.min(np.array(per_clique_end)))
                    ciend[0] = min_end - end
                    ciend[1] = max_end - end
                    ciend = ','.join([str(i) for i in ciend])
                    if type_sv_name != 'TRA':
                        avglen = int(np.array([abs(i.svlen) for i in per_clique_list]).mean())
                        if type_sv_name == 'DEL':
                            avglen = -avglen
                    else:
                        avglen = 0
                    supp_vec = ['0'] * len(vcf_file_list)
                    supp_info = defaultdict(list)
                    for per in per_clique_list:
                        if per.patient in vcf_file_list:
                            idx = vcf_file_list.index(per.patient)
                            supp_vec[idx] = '1'
                            genotype = per.genotype
                            length = abs(per.svlen)
                            dr = per.dr
                            dv = per.dv
                            strands = per.strands
                            start_end = '%s~%s--%s~%s' % (per.chr1, per.start, per.chr2, per.end)
                            supp_line = '%s:NA:%s:%s,%s:%s:%s:%s' % \
                                        (genotype, length, dr, dv, strands, type_sv_name, start_end)
                            supp_info[per.patient].append(supp_line)
                    supp_line_list = []
                    for i in vcf_file_list:
                        if i not in supp_info:
                            supp_line_list.append('./.:NaN:0:0,0:--:NaN:NaN')
                        else:
                            supp_line_list.append(';'.join(supp_info[i]))
                    supp_vec = ''.join(supp_vec)
                    supp = supp_vec.count('1')
                    supp_info1 = '\t'.join(supp_line_list)
                    info_all = ['SUPP=%s' % supp, 'SUPP_VEC=%s' % supp_vec, 'AVGLEN=%s' % avglen,
                                'SVTYPE=%s' % type_sv_name, 'SVMETHOD=cliqueSV_maxflow', 'CHR2=%s' % chr2,
                                'END=%s' % end, 'CIPOS=%s' % cipos, 'CIEND=%s' % ciend, 'STRANDS=%s' % strands_max]
                    info_all = ';'.join(info_all)
                    newline = [chr1, start, 'sv_id', 'N', type_sv, '.', 'PASS', info_all,
                               'GT:PSV:LN:DR:ST:TY:CO', supp_info1]
                    total_wait_dict[chrom].append(newline)
    # sort by start
    for total_wait_idx in total_wait_dict:
        total_wait_chrom = total_wait_dict[total_wait_idx]
        total_wait_chrom = sorted(total_wait_chrom, key=lambda x: x[1])
        total_wait_dict[total_wait_idx] = total_wait_chrom
    n = 0
    for i in total_wait_dict:
        for i_per in total_wait_dict[i]:
            sv_id = '%s00%sSUR' % (i_per[4][1:-1], n)
            n += 1
            i_per[2] = sv_id
            i_per[1] = str(i_per[1])
            f_output.write('\t'.join(i_per) + '\n')
    f_output.close()

def log_output(all_sv_form, all_sv_block, all_sv_merge, output):
    # type    chr1    chr2    start    end    patient    sv_id    genotype
    if output.endswith('vcf'):
        output = output[:-4]
    shuchu = open('%s_clique_maxflow.log' % output, 'w')
    num = 0
    for type_sv in all_sv_merge:
        for chrom in all_sv_merge[type_sv]:
            for per_block in all_sv_merge[type_sv][chrom]:
                shuchu.write('block_%s_%s\n' % (type_sv, num))
                for idx, per_clique in enumerate(per_block):
                    shuchu.write('block_%s_%s_clique_%s\n' % (type_sv, num, idx))
                    for per_sv in per_clique:
                        newline = [per_sv.sv_type, per_sv.chr1, per_sv.chr2, str(per_sv.start),
                                   str(per_sv.end), per_sv.patient, per_sv.sv_id, per_sv.genotype]
                        newline = '\t'.join(newline) + '\n'
                        shuchu.write(newline)
                num += 1
    shuchu.close()
    shuchu = open('%s_block_maxflow.log' % output, 'w')
    num = 0
    for type_sv in all_sv_block:
        for chrom in all_sv_block[type_sv]:
            for per_block in all_sv_block[type_sv][chrom]:
                shuchu.write('block_%s_%s\n' % (type_sv, num))
                num += 1
                for per_sv in per_block:
                    newline = [per_sv.sv_type, per_sv.chr1, per_sv.chr2, str(per_sv.start),
                               str(per_sv.end), per_sv.patient, per_sv.sv_id, per_sv.genotype]
                    newline = '\t'.join(newline) + '\n'
                    shuchu.write(newline)
    shuchu.close()
    shuchu = open('%s_total_maxflow.log' % output, 'w')
    for type_sv in all_sv_form:
        for chrom in all_sv_form[type_sv]:
            for per_sv in all_sv_form[type_sv][chrom]:
                newline = [per_sv.sv_type, per_sv.chr1, per_sv.chr2, str(per_sv.start),
                           str(per_sv.end), per_sv.patient, per_sv.sv_id, per_sv.genotype]
                newline = '\t'.join(newline) + '\n'
                shuchu.write(newline)
    shuchu.close()   
     
def main():
    args=parseargs()
    start_time1 = time.time()
    #args.vcf_list = '/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/Sniffles.reformat.filter.vcf,/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/cuteSV.reformat.filter.vcf,/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/pbsv.reformat.filter.vcf,/home1/Laisenying/Data-analysis/projects/PhageSV/Test/Simulator/SimSV/Test2_1000/PBSim_reads/replicate2_6e5/CAST/new_vcf/svim.reformat.filter.vcf'
    vcf_file_list, contigs, all_sv_form = read_vcf(args.vcf_list, args.allele_freq)    
    print('-----------------[Read end!]', time.time() - start_time1, 's-----------------')
    all_sv_block = break2block(all_sv_form, args.interval)
    print('-----------------[Block end!]', time.time() - start_time1, 's-----------------')
    all_sv_merge = block2clique2merge(all_sv_block, args.extended_dis, args.overlap_perc, args.workers)
    print('-----------------[Final end!]', time.time() - start_time1, 's-----------------')
    output_vcf(all_sv_merge, vcf_file_list, contigs, args.outfile)
    print('-----------------[Output end!]', time.time() - start_time1, 's-----------------')
    log_output(all_sv_form, all_sv_block, all_sv_merge, args.outfile)    
    
if __name__ == '__main__':
    main()
