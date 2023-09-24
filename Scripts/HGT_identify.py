import sys
from collections import Counter
from datetime import datetime
from glob import glob
from ete3 import PhyloTree, TreeStyle, NCBITaxa, NodeStyle, faces, AttrFace, CircleFace
from numpy import mean,std
import numpy as np
from itertools import chain
import numpy as np
import pandas as pd
import os


# 4c. Find hgts with unknown directionality
def unkn_detect(t,sub_PC_inf):
    for clade in chain(t.get_monophyletic(values=["p"], target_attr="group"),t.get_monophyletic(values=["p","PtoV"], target_attr="group")):
        # look for eukaryotic clades that are sister to viral clades but did not get called as hgts
        try:
            s = []
            for sister in clade.get_sisters():
                s.append(clade_identifier([leaf.group for leaf in sister.get_leaves()]))
            sister_group = group_identifier(s)
            p = []
            for sister in clade.up.get_sisters():
                p.append(clade_identifier([leaf.group for leaf in sister.get_leaves()]))
            parent_group = group_identifier(p)
        except:
            continue
        if (sister_group == 'v' or sister_group == 'prob_v') and (parent_group == 'unkn'):
            note = []
            # if found, change the variable found from 'no' to 'yes'
            global found
            found = 'yes'
            # record the sister group so that we can see the donor
            sisters = []
            for sister in clade.get_sisters():
                for leaf in sister.get_leaves():
                    if sub_PC_inf.loc[leaf.name,'Tax']=='Viral':
                        sisters.append(leaf.name)
            #sister_tree=ncbi.get_topology(sisters)
            recipient_taxa = ''
            for leaf in clade.get_leaves():
                recipient_taxa = recipient_taxa + ',' + leaf.name
            donor_taxa = ''
            sister_branch_length = []
            for sister in clade.get_sisters():
                sister_branch_length.append(sister.dist)
                for leaf in sister.get_leaves():
                    donor_taxa = donor_taxa + ',' + leaf.name
            # make notes:
            if (clade.dist > long_branch):
                note.append('Long branching recipient')
            elif (mean(sister_branch_length) > long_branch):
                note.append('Long branching donor')
            elif sister_group == 'prob_p':
                note.append('Sister unclear')
            elif parent_group == 'prob_p':
                note.append('Cousin unclear')
            if len(note) == 1:
                notes = note[0]
            elif len(note) > 1:
                notes = note[0]
                for n in note[1:]:
                    notes = notes + ',' + n
            else:
                notes = 'NA'
            # replace group information from E to unkn
            for leaf in clade.get_leaves():
                leaf.add_feature("group","unkn")
            # output results: COG\tdirection\ttaxid_recip\tname_recip\trecipients\ttaxid_don\tname_don\tdonors\tnotes\n
            if (clade.dist > long_branch) or (mean(sister_branch_length) > long_branch):
                out.write(fname.split('/')[-1].split("_")[0]+'\tunkn\t'+recipient_taxa.strip(',')+'\tUnresolved\t'+notes+'\t'+str(clade.dist)+'\t'+str(clade.up.support)+'\t'+str(clade.up.up.support)+'\n')
            else:
                out.write(fname.split('/')[-1].split("_")[0]+'\tunkn\t'+recipient_taxa.strip(',')+'\t'+donor_taxa.strip(',')+'\t'+notes+'\t'+str(clade.dist)+'\t'+str(clade.up.support)+'\t'+str(clade.up.up.support)+'\n')


# 4a. Define group identifier functions incase there are polytomies
# Use this to assign an identification to a sister clade
def clade_identifier(g):
    if (g.count('v')+g.count('VtoE'))/float(len(g)) == 1.0:
        id = 'v'
    elif (g.count('v')+g.count('VtoE'))/float(len(g)) >= 0.8:
        id = 'prob_v'
    elif (g.count('p')+g.count('PtoV'))/float(len(g)) == 1.0:
        id = 'p'
    elif (g.count('p')+g.count('PtoV'))/float(len(g)) >= 0.8:
        id = 'prob_p'
    else:
        id = 'unkn'
    return id
    
# Use this to assign an identification to a sister group that is comprised of a polytomy
def group_identifier(g):
    if (g.count('v'))/float(len(g)) == 1.0:
        id = 'v'
    elif (g.count('v')+g.count('prob_v'))/float(len(g)) > 0.66:
        id = 'prob_v'
    elif (g.count('p'))/float(len(g)) == 1.0:
        id = 'p'
    elif (g.count('p')+g.count('prob_p'))/float(len(g)) > 0.66:
        id = 'prob_p'
    else:
        id = 'unkn'
    return id
    
# 4b. EtoV detection
def PtoV_detect(t,sub_PC_inf):
    for clade in chain(t.get_monophyletic(values=["v"], target_attr="group"),t.get_monophyletic(values=["v","VtoP"], target_attr="group")):
        try:
            s = []
            for sister in clade.get_sisters():
                s.append(clade_identifier([leaf.group for leaf in sister.get_leaves()]))
            sister_group = group_identifier(s)
            p = []
            for sister in clade.up.get_sisters():
                p.append(clade_identifier([leaf.group for leaf in sister.get_leaves()]))
            parent_group = group_identifier(p)
        except:
            continue
        if (sister_group == 'p' or sister_group == 'prob_p') and (parent_group == 'p' or parent_group == 'prob_p'):
            note = []
            # if found, change the variable found from 'no' to 'yes'
            global found
            found = 'yes'
            # record the sister group so that we can see the donor
            sisters = []
            for sister in clade.get_sisters():
                for leaf in sister.get_leaves():
                    if sub_PC_inf.loc[leaf.name,'Tax']=='Prokaryotic':
                        sisters.append(leaf.name)
            #sister_tree=ncbi.get_topology(sisters)
            # write results to summary file
            if 'PtoE' in [leaf.group for leaf in clade.get_leaves()]:
                l = []
                for leaf in clade.get_leaves():
                    if leaf.group == 'v':
                        l.append(leaf.name)
                clade.taxid = 10239
                note.append('Prokaryotic transduction')
            else:
                clade.taxid = 10239
            # write results to summary file
            recipient_taxa = ''
            for leaf in clade.get_leaves():
                recipient_taxa = recipient_taxa + ',' + leaf.name
            donor_taxa = ''
            sister_branch_length = []
            for sister in clade.get_sisters():
                sister_branch_length.append(sister.dist)
                for leaf in sister.get_leaves():
                    donor_taxa = donor_taxa + ',' + leaf.name
            # make notes:
            if (clade.dist > long_branch):
                note.append('Long branching recipient')
            elif (mean(sister_branch_length) > long_branch):
                note.append('Long branching donor')
            elif sister_group == 'prob_p':
                note.append('Sister unclear')
            elif parent_group == 'prob_p':
                note.append('Cousin unclear')
            if len(note) == 1:
                notes = note[0]
            elif len(note) > 1:
                notes = note[0]
                for n in note[1:]:
                    notes = notes + ',' + n
            else:
                notes = 'NA'
            # replace group information from V to EtoV
            for leaf in clade.get_leaves():
                leaf.add_feature("group","PtoV")
            #output results: COG\tdirection\ttaxid_recip\tname_recip\trecipients\ttaxid_don\tname_don\tdonors\tnotes\n
            if (clade.dist > long_branch) or (mean(sister_branch_length) > long_branch):
                out.write(fname.split('/')[-1].split("_")[0]+'\tPtoV\t' + recipient_taxa.strip(',')+'\tUnresolved\t'+notes+'\t'+str(clade.dist)+'\t'+str(clade.up.support)+'\t'+str(clade.up.up.support)+'\n')
            elif sister_group == 'prob_p':
                out.write(fname.split('/')[-1].split("_")[0]+'\tPtoV\t'+recipient_taxa.strip(',')+'\tUnresolved\t'+notes+'\t'+str(clade.dist)+'\t'+str(clade.up.support)+'\t'+str(clade.up.up.support)+'\n')
            else:
                out.write(fname.split('/')[-1].split("_")[0]+'\tPtoV\t'+recipient_taxa.strip(',')+'\t'+ donor_taxa.strip(',')+'\t'+notes+'\t'+str(clade.dist)+'\t'+str(clade.up.support)+'\t'+str(clade.up.up.support)+'\n')


# 4a. VtoE detection
def VtoP_detect(t,sub_PC_inf):
    for clade in chain(t.get_monophyletic(values=["p"], target_attr="group"),t.get_monophyletic(values=["p","PtoV"], target_attr="group")):
        try:
            s = []
            for sister in clade.get_sisters():
                s.append(clade_identifier([leaf.group for leaf in sister.get_leaves()]))
            sister_group = group_identifier(s)
            p = []
            for sister in clade.up.get_sisters():
                p.append(clade_identifier([leaf.group for leaf in sister.get_leaves()]))
            parent_group = group_identifier(p)
        except:
            continue
        if (sister_group == 'v' or sister_group == 'prob_v') and (parent_group == 'v' or parent_group == 'prob_v'):
            note = []
            # if found, change the variable found from 'no' to 'yes'
            global found
            found = 'yes'
            # record the sister group so that we can see the donor
            sisters = []
            for sister in clade.get_sisters():
                for leaf in sister.get_leaves():
                    if sub_PC_inf.loc[leaf.name,'Tax']=='Viral':
                        sisters.append(leaf.name)
            #sister_tree=ncbi.get_topology(sisters)
            # recipient taxa
            if 'PtoV' in [leaf.group for leaf in clade.get_leaves()]:
                l = []
                for leaf in clade.get_leaves():
                    if leaf.group == 'p':
                        l.append(leaf.name)
                clade.taxid = 2
            else:
                clade.taxid = 2
            # get results for summary file
            recipient_taxa = ''
            for leaf in clade.get_leaves():
                recipient_taxa = recipient_taxa + ',' + leaf.name
            donor_taxa = ''
            sister_branch_length = []
            for sister in clade.get_sisters():
                sister_branch_length.append(sister.dist)
                for leaf in sister.get_leaves():
                    donor_taxa = donor_taxa + ',' + leaf.name
            # check the origin of the gene - eukaryotic transduction, bacterial, or viral?
            og_clade = clade
            x = 'v'
            while x == 'v':
                try:
                    for sister in clade.get_sisters():
                        s = []
                        s.append(clade_identifier([leaf.group for leaf in sister.get_leaves()]))
                        sister_group = group_identifier(s)
                        if sister_group == 'prob_v':
                            x = 'v'
                        else:
                            x = sister_group
                    if x == 'v':
                        if clade.up:
                            clade = clade.up
                        else:
                            break
                except:
                    x = 'unkn'
            origin = 'origin = ' + x
            clade = og_clade
            # make notes:
            note.append(origin)
            if (clade.dist > long_branch):
                note.append('Long branching recipient')
            elif (mean(sister_branch_length) > long_branch):
                note.append('Long branching donor')
            elif sister_group == 'prob_v':
                note.append('Sister unclear')
            elif parent_group == 'prob_v':
                note.append('Cousin unclear')
            if len(note) == 1:
                notes = note[0]
            elif len(note) > 1:
                notes = note[0]
                for n in note[1:]:
                    notes = notes + ',' + n
            else:
                notes = 'NA'
            # replace group information from E to VtoE
            for leaf in clade.get_leaves():
                leaf.add_feature("group","VtoP")
            # output results: COG\tdirection\ttaxid_recip\tname_recip\trecipients\ttaxid_don\tname_don\tdonors\tnotes\n
            if (clade.dist > long_branch) or (mean(sister_branch_length) > long_branch):
                out.write(fname.split('/')[-1].split("_")[0]+'\tVtoP\t'+recipient_taxa.strip(',')+'\tUnresolved\t'+notes+'\t'+str(clade.dist)+'\t'+str(clade.up.support)+'\t'+str(clade.up.up.support)+'\n')
            elif sister_group == 'prob_v':
                out.write(fname.split('/')[-1].split("_")[0]+'\tVtoP\t'+recipient_taxa.strip(',')+'\tUnresolved\t'+notes+'\t'+str(clade.dist)+'\t'+str(clade.up.support)+'\t'+str(clade.up.up.support)+'\n')
            else:
                out.write(fname.split('/')[-1].split("_")[0]+'\tVtoP\t'+recipient_taxa.strip(',')+'\t'+donor_taxa.strip(',')+'\t'+notes+'\t'+str(clade.dist)+'\t'+str(clade.up.support)+'\t'+str(clade.up.up.support)+'\n')
                                 
# 3. Define group annotation function
def annotate_leaf_group(t,sub_PC_inf):
    for leaf in t.get_leaves():
         #check if it's a cellular organism, if not - v
         sub_PC_inf.loc[leaf.name,]
         if sub_PC_inf.loc[leaf.name,'Tax'] == 'Prokaryotic':
             leaf.add_feature("group", "p")
         elif sub_PC_inf.loc[leaf.name,'Tax'] == 'Viral':
             leaf.add_feature("group", "v")
        
def custom_check_monophyly(t, values, attr):
    '''Returns True if all attributes (attr) in a tree (t) are monophyletic for the provided set of values. The most monophyletic partition is also returned.'''
    content = t.get_cached_content()
    alltips = content[t]
    targets = set([n for n in alltips if getattr(n, attr) in values])
    smallest = None
    for n, leaves in content.items():
        if targets.issubset(leaves) and (not smallest or len(leaves) < len(smallest)):
            smallest = leaves
        else:
            # if the node itself is not defining the monophyly, break on through
            # to the other side... and check if the monophyly is there
            other_side = alltips - leaves
            if targets.issubset(other_side) and (not smallest or len(other_side) < len(smallest)):
                smallest = other_side
    return len(smallest) == len(targets), smallest

# 2. Define layout parameters for tree visualization
def my_layout(node):
    # Layout function
    node.img_style["size"] = 0
    if node.is_leaf():
        name = faces.TextFace("%s {%s}" % (node.name, node.sci_name))
        # colour viruses yellow
        if node.group == "v":
            name.background.color = "yellow"
        # colour eukaryotes green
        elif node.group == "p":
            name.background.color = "green"
        elif node.group == 'VtoP':
        # colour virus to euk Hgt red
            name.background.color = 'red'
        elif node.group == 'PtoV':
        # colour euk to virus Hgt blue
            name.background.color = 'blue'
        # colour hgts with unknown direction orange
        elif node.group == 'unkn':
            name.background.color = 'orange'
        # colour bacteria cyan
        else:
            name.background.color = "cyan"
        faces.add_face_to_node(name, node, 0)

date = str(datetime.now()).split(' ')[0].replace('-','_')
out = open('Phylogeny_based/HGT_result/NVOG_hgt_summary.' + date + '.tab','w')
failed_out=open('Phylogeny_based/HGT_result/failed_hgt_PC.' + date + '.tab','w')
Analyzied_out=open('Phylogeny_based/HGT_result/Parsed_hgt_PC.' + date + '.tab','w')
out.write('COG\tdirection\trecipients\tdonors\tnotes\tbranch_length\tsupport1\tsupport2\n')
                                                  
PC_inf = pd.read_csv("BacViral_protein_family_inf_filter_new.tsv",sep="\t",index_col=0)
n=1
PC_list = pd.read_csv("all_PC_new.list",header=None)
PC_list = list(PC_list[0])
for PC in PC_list:
    fname="Phylogeny_based/IQtree_output/"+PC+"_align_filter.fa.treefile"
    if not os.path.exists(fname):
        continue
    print(fname)
    sub_PC_inf = PC_inf.loc[PC_inf['PC']==PC,]
    sub_PC_inf['protein_name'] = ["_".join(x.split("~")) for x in sub_PC_inf['protein']]
    sub_PC_inf.index = sub_PC_inf['protein_name']
    i = int(((n)/float(len(glob("Phylogeny_based/IQtree_output/*treefile")))*100))
    sys.stdout.write("[%-100s] %d%%" % ('='*i, i) + '\t' + str(n)+'/'+str(len(glob("*treefile"))))
    sys.stdout.flush()
    # load in tree
    try:
        t = PhyloTree(fname, sp_naming_function=lambda name: int(name.split(".")[0]))
    except:
        failed_out.write(PC+"\n")
        continue
    # collapse nodes with support less than 80 (SH-alrt IQ-Tree)
    collapse_value = 80
    for node in t.get_descendants():
        if not node.is_leaf() and node.support <= collapse_value:
            node.delete()
    # if mad rooting available, use it, otherwise use midpoint rooting
    if fname+'.rooted' in glob('Phylogeny_based/IQtree_output/*.rooted'):
        try:
            MAD = PhyloTree(fname+'.rooted', sp_naming_function=lambda name: int(name.split(".")[0]))
        except:
            failed_out.write(PC+"\n")
            continue
        try:
            madroot = [otu.name for otu in MAD.get_children()[0].get_leaves()]
            if len(madroot) > 1:
                t.set_outgroup(t.get_common_ancestor(madroot))
            else:
                t.set_outgroup(madroot[0])
        except:
            t.set_outgroup([otu.name for otu in MAD.get_children()[1].get_leaves()][0])
            t.set_outgroup(t.get_common_ancestor(madroot))
    else:
        midpoint = t.get_midpoint_outgroup()
        if midpoint: t.set_outgroup(midpoint)
    # annotate taxa and group (v, e, p)
    annotate_leaf_group(t,sub_PC_inf)
    # check for monophyly - if only E and V present and they are monophyletic, call as an unknown because rooting is ambiguous - proportions may suggest direction though
    # If #Euk sps / # Viral sps > 80% of all alignments (>49) or < 20% (0.15) of all alignments continue. Otherwise call as unknown.
    tree_otus = [leaf.name for leaf in t.get_leaves()]
    taxa = list(sub_PC_inf.loc[tree_otus,'Tax'])
    sps = ['P' if x=='Prokaryotic' else 'V' for x in taxa]
    tally = Counter(sps)
    if tally['V'] == 0:
        continue
    PoV = float(tally['P'])/tally['V']
    # if the tree only contains euks and viruses AND they can be split into monophyletic groups AND the species proportions are not extreme (between 20-80 percentile #E/#V) - call unknowns
    if (set(sps) == {'P','V'}) and (custom_check_monophyly(t,values=set(['v']), attr='group')[0] == True) and (PoV > 0.15384615384615385 and PoV < 49.0):
        for clade in t.get_monophyletic(values=["p"], target_attr="group"):
            pro_clade = clade
        t.set_outgroup(pro_clade)
    # determine long branch as extreme outlier value (Outlier = Q3 + 3*IQR)
    branches = np.asarray(list(set([node.dist for node in t.iter_descendants()])))
    long_branch = (np.percentile(branches,75))+(3*(np.percentile(branches,75))-(np.percentile(branches,25)))
    # second, look for virus to/from eukaryote hgts
    i = 0
    while i <= 3:
        PtoV_detect(t,sub_PC_inf)
        VtoP_detect(t,sub_PC_inf)
        i += 1
    # lastly, call any unknown hgts
    unkn_detect(t,sub_PC_inf)
    #if found == 'yes' and print_trees =='yes':
    #    t.render(fname+'.pdf', w=183, tree_style = S, layout = my_layout)
    n += 1
    Analyzied_out.write(PC+"\n")
    new_file = open("Phylogeny_based/IQtree_output/"+PC+".newtreefile","w")
    new_file.write(t.write()+"\n")
    new_file.close()
out.close()
