#!/usr/bin/env python
# coding: utf-8

# In[71]:


import glob
from Bio import SeqIO, Seq
from itertools import combinations


# In[2]:


def read_filter_write(file_list, final_directory = "filtered_yeast"):
    for file in file_list:
        with open(file, 'rt') as file_to_check:
            proteins = list(SeqIO.parse(file_to_check, "fasta"))
            if len(proteins) == 9:
                genomes = []
                for protein in proteins:
                    genome_name = protein.id[0:4]
                    genomes.append(genome_name)
                    protein.id = protein.id[0:4]
                    
                if len(genomes) == len(set(genomes)):
                    with open("{}/{}".format(final_directory, file[7:]), 'w') as file_to_write:
                        SeqIO.write(proteins, file_to_write, "fasta")
        file_to_check.close()

files = glob.glob("famdir/*.fasta")
read_filter_write(files)


# In[2]:


from Bio.Align.Applications import ClustalwCommandline


# In[ ]:


handler_clustal = '/usr/bin/clustalw'

filtered_files = glob.glob("filtered_yeast/*.fasta")

for file in filtered_files:
    clustal_alignment = ClustalwCommandline(handler_clustal, 
                                            infile= file)
    out_log, err_log = clustal_alignment()


# In[20]:


from Bio.Phylo.TreeConstruction import *
from Bio import AlignIO, Phylo
from matplotlib import pyplot as plt


# In[21]:


def create_tree(infile, outfile = None, algorithm = 'upgma'):
    aln = AlignIO.read(infile, 'clustal')
    constructor = DistanceTreeConstructor()
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    
    if algorithm == 'upgma':
        tree = constructor.upgma(dm)
    elif algorithm == 'nj':
        tree = constructor.nj(dm)
        
    return(tree)


# In[22]:


def tree_to_dictionary(tree):
    inner_clade_dic = {}

    for inner_clade in tree.get_nonterminals():
        inner_clade_dic[inner_clade.name] = set([])
        
    leaves_lst = tree.get_terminals()
    
    for l1 in range(len(leaves_lst)):
        leaves_lst_without_lth = leaves_lst[0:l1] + leaves_lst[(l1+1):]
        for l2 in range(len(leaves_lst_without_lth)):
            leaf1 = leaves_lst[l1]
            leaf2 = leaves_lst_without_lth[l2]
            clade = tree.common_ancestor(leaf1, leaf2)
            inner_clade_dic[clade.name].update([leaf1.name])
            inner_clade_dic[clade.name].update([leaf2.name])
            
    #add leaves
    clade_dic = {d:list(inner_clade_dic[d]) for d in inner_clade_dic}
    
    for leaf in tree.get_terminals():
        clade_dic[leaf.name] = [leaf.name]
        
    for clade in clade_dic:
        clade_dic[clade].sort()
            
    return(clade_dic)


# In[23]:


def clade_count(aligned_files):
    
    counts = {}
    tree_number = 0
    
    for f in aligned_files:
        
        tree_number += 1
        
        tr = create_tree(f)
        t_dictionary = tree_to_dictionary(tr)
        
        for k in t_dictionary:
            key = ""
            for el in t_dictionary[k]:
                key += el + '_'
            
            if key in counts:
                counts[key] += 1
            else:
                counts[key] = 1
                
    return(counts, tree_number)


# In[34]:


def dictionary_converter(dictionary):
    lst = []
    
    for key in dictionary:
        new_key = key.split('_')[:-1]
        lst.append((new_key, len(new_key)))
        
    return(sorted(lst, key = lambda x: x[1]))


# In[117]:


def generate_consensus(counts, no, typeof):
    
    frequency_dictionary = {}
    dic = {}
    
    for c in counts:
        frequency_dictionary[c] = counts[c]/no
        
    if typeof == 'majority':
        dic = {k:frequency_dictionary[k] for k in frequency_dictionary if frequency_dictionary[k] > 0.5}
    elif typeof == 'strict':
        dic = {k:frequency_dictionary[k] for k in frequency_dictionary if frequency_dictionary[k] == 1}
    
    sorted_dic = dictionary_converter(dic)
    
    tree = [s[0] for s in sorted_dic]
    
    return(tree)


# In[25]:


aln_files = glob.glob('filtered_yeast/*.aln')


# In[26]:


cnts, number = clade_count(aln_files)


# In[115]:


majority_consensus = generate_consensus(cnts, number, 'majority')
for maj in majority_consensus:
    print(maj)


# In[121]:


majority_tree = "(((DEHA)), ((CAGL), (ERGO), (KLLA), (SACE), (KLTH, SAKL), (ZYRO)), ((YALI)))"


# In[118]:


strict_consensus = generate_consensus(cnts, number, 'strict')
for stri in strict_consensus:
    print(stri)


# In[123]:


strict_tree = "CAGL, DEHA, ERGO, KLLA, KLTH, SACE, SAKL, YALI, ZYRO"


# In[129]:


from io import StringIO
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [15, 8]


# In[139]:


tree = Phylo.read(StringIO(majority_tree), "newick")
Phylo.draw(tree)


# In[140]:


tree = Phylo.read(StringIO(strict_tree), "newick")
Phylo.draw(tree)

