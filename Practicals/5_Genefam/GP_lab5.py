#!/usr/bin/env python
# coding: utf-8

# # Genomika porównawcza
# ## Laboratorium 5
# ### Igor Filipiuk

# In[1]:


from Bio import SeqIO, Seq, Phylo, AlignIO
from Bio.Align.Applications import ClustalwCommandline


# ### Process protein sequences of 18 species

# In[2]:


hemoglobins = list(SeqIO.parse(open("Hemoglobina.fasta"), "fasta"))

for h in hemoglobins:
    h.id = h.description.split()[-2][1::] + "_" + h.description.split()[-1][0:-1]
    h.name = h.description.split()[-2][1::] + "_" + h.description.split()[-1][0:-1]

with open("Hemoglobina_speciesnames.fasta", "w") as text_file:
    SeqIO.write(hemoglobins, text_file, "fasta")


# ### Align sequences usin clustal algorithm

# In[3]:


handler_clustal = '/usr/bin/clustalw'
clustal_alignment = ClustalwCommandline(handler_clustal, 
                                        infile='Hemoglobina_speciesnames.fasta')


# In[4]:


out_log, err_log = clustal_alignment()


# In[5]:


aln = AlignIO.read('Hemoglobina_speciesnames.aln', 'clustal')
tree = Phylo.read('Hemoglobina_speciesnames.dnd', 'newick')


# ### Build tree using Maximum Likelihoog (ML)

# In[6]:


from Bio.Phylo.TreeConstruction import *
from Bio import AlignIO

from Bio.Phylo.Applications import PhymlCommandline

AlignIO.convert("Hemoglobina_speciesnames.aln", "clustal", 
                "Hemoglobina_speciesnames.phy", "phylip-relaxed")
phyml_cl = PhymlCommandline(input='Hemoglobina_speciesnames.phy', 
                            datatype = 'aa', alpha = 'e', bootstrap = '120')
out_log, err_log = phyml_cl()


# In[7]:


mltree = Phylo.read('Hemoglobina_speciesnames.phy_phyml_tree.txt', 'newick')

Phylo.draw_ascii(mltree)


# ### Build tree using Neighbour Joining (NJ)

# In[8]:


aln = AlignIO.read('Hemoglobina_speciesnames.aln', 'clustal')
constructor = DistanceTreeConstructor()
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)

njtree = constructor.nj(dm)
Phylo.draw_ascii(njtree)


# ### Build tree using Maximum Parsimny (NJ)

# In[9]:


aln = AlignIO.read('Hemoglobina_speciesnames.aln', 'clustal')
starting_tree = njtree
scorer = ParsimonyScorer()
searcher = NNITreeSearcher(scorer)
constructor = ParsimonyTreeConstructor(searcher, starting_tree)

pars_tree = constructor.build_tree(aln)
Phylo.draw_ascii(pars_tree)


# In[ ]:


"""All trees generated be each of the methods (ML, NJ, MP) well resemble general taxonomy of species. 
Horse species (Equus Burchelli, Equus przewalskii, Equus hemionus, Equus Zebra, Equus asinus) are 
clustered together in every of three trees. Fish like cod (Gadus morhua), tuna (Thunnus thynuns) 
are clustered together, are the furtherst from any other species in the tree and at the same time the closest
species to them are frog (Xenopus tropicalis) and viper (Crotalus_adamanteus).

One wierd anomaly that one can find is that beaver (Castor fiber) are closer related to fish, reptiles and 
amphibians than to other mammals (especially do other rodents like mice)."""

