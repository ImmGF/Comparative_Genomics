#!/usr/bin/env python
# coding: utf-8


import os
import pprint
import shutil
from os import mkdir
import re
from subprocess import Popen, PIPE
from io import StringIO
from tqdm import tqdm
from Bio import Entrez
from Bio import SeqIO


class Protein(object):
    def __init__(self, locus, accession):
        self.locus = locus
        self.accession = accession

    def __str__(self):
        return self.accession

    def __repr__(self):
        return self.accession


families = dict()
with open("gene_families.txt") as gene_families_file:
    lines = gene_families_file.readlines()
    length = len(lines)
    header = ""
    for i in range(1, length):
        if re.match("^=+$", lines[i]):
            header = "_".join(lines[i - 1].strip().split())
            families[header] = []
        else:
            match = re.match("(\S+_\S+)\s+\((\S{6})\)\s+(.*)", lines[i])
            if match:
                families[header].append(
                    Protein(match.group(1), match.group((2))))


print(families)



familiesTrimmedTo20Proteins = {k: v[:20]
                               for k, v in families.items() if len(v) >= 20}

selectedFamilies = {tup[0]: tup[1]
                    for tup in list(familiesTrimmedTo20Proteins.items())[:5]}



Entrez.email = "if290707@student.mimuw.edu.pl"


def getFasta(selectedFamilies):
    for family, proteins in tqdm(selectedFamilies.items()):
        for protein in tqdm(proteins):
            handle = Entrez.efetch(db="protein", id=protein.accession,
                                   rettype="fasta", retmode="text")
            record = SeqIO.read(StringIO(handle.read()), "fasta")
            record.id = f"{protein.accession}|{protein.locus}|{family}\n"
            record.description = ""
            yield record


with open("proteins.fasta", "w") as fasta:
    SeqIO.write(getFasta(selectedFamilies), fasta, "fasta")



print("exprected clusters:")
pp = pprint.PrettyPrinter(indent=4)
pp.pprint(selectedFamilies)



shutil.rmtree("uclust", True)
mkdir("uclust")
Popen("usearch -cluster_fast proteins.fasta -id 0.25 -clusters uclust/cluster".split(),).wait()

print("UCLUST:")
for file in sorted(os.listdir(os.fsencode("uclust"))):
    filename = os.fsdecode(file)
    print(filename + ":")
    pp.pprint([seqRecord.id for seqRecord in SeqIO.parse(
        "uclust/" + filename, "fasta")])


# %%
shutil.rmtree("cd-hit", True)
mkdir("cd-hit")
Popen("cd-hit -i proteins.fasta -o cd-hit/output -c 0.4 -n 2 -d 40 -T 4".split()).wait()

print("CD-HIT:")
with open("cd-hit/output.clstr") as output:
    pp.pprint(output.readlines())


# %%
shutil.rmtree("mcl", True)
mkdir("mcl")
Popen("makeblastdb -in proteins.fasta -dbtype prot".split()).wait()
Popen("blastp -db proteins.fasta -query proteins.fasta -out mcl/blast_results.out -outfmt 6".split()).wait()

with open("mcl/blast_results.abc", "w") as blast_results_abc:
    Popen("cut -f 1,2,11 mcl/blast_results.out".split(),
          stdout=blast_results_abc).wait()
Popen("mcxload -abc mcl/blast_results.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o mcl/blast_results.mci -write-tab mcl/blast_results.tab".split()).wait()
Popen("mcl mcl/blast_results.mci -I 5 -use-tab mcl/blast_results.tab".split()).wait()

# %%
print("mcl:")
with open("out.blast_results.mci.I50") as output:
    pp.pprint(output.readlines())


# %%
shutil.rmtree("blastclust", True)
mkdir("blastclust")
Popen("blastclust -i proteins.fasta -d proteins.fasta -o blastclust/out -S 50 -a 4".split()).wait()
print("blastclust:")
with open("blastclust/out") as output:
    pp.pprint(output.readlines())
os.chdir("..")

