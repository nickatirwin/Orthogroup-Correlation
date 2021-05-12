####################################
### PROTEIN CORRELATION ANALYSIS ###
####################################

"""

Correlate a given protein/phenotype 
with protein families.

Arguments:

-q query protein (.fasta) (temp: sys.argv[1])
-p curated protein family (.fasta)
-c manually generated clade file (.tab)
-d custom mcl dumpfile/orthogroup (.tab)
-t number of threads (temp: sys.argv[2])
-n number of permutations (temp: sys.argv[3])

Usage: 

Requirements:
BLAST
faSomeRecords
mafft
trimAl
hmmer
iqtree
ete3

"""

import sys
import os
from subprocess import call
from ete3 import PhyloTree, NCBITaxa
from collections import Counter
from difflib import SequenceMatcher

ncbi = NCBITaxa()

### 1. Map the query protein to the provided protein families

call('mkdir temp', shell = True)

# BLAST the query protein to the provided dataset

print('\nMapping query to dataset [1/6]')
call('diamond blastp --quiet --sensitive --query ' + sys.argv[1] + ' --db files/eukaryotes.fasta.dmnd --max-target-seqs 20 --outfmt 6 --evalue 1e-5 --threads ' + sys.argv[2] + ' --out temp/query.blastp', shell = True, stdout=open(os.devnull, 'wb'))
try:
    mapped_seq = open('temp/query.blastp','r').readlines()[0].split('\t')[1]
except:
    sys.exit('\nError: The query did not have a match in the provided dataset. Try a custom dataset.\n')

# identify the COG with the mapped protein

dump = open('files/dump.out.eukaryotes.allVSall.dmnd.blastp.mci.I14','r').readlines()

m = 0
while m < 1:
    for line in dump:
        if mapped_seq in line.split('\t'):
            COG = line.split('\t')
            m += 1

out = open('temp/queryCOG','w')
sep = '\n'
out.write(sep.join(COG))
out.close()

### 2. Find orthologues using HMMs

# extract the COG, build an HMM, search the provided dataset, and make a phylogeny

call('faSomeRecords files/eukaryotes.fasta temp/queryCOG temp/queryCOG.fasta', shell = True, stdout=open(os.devnull, 'wb'))
print('\nAligning intial cluster and conducting HMM search [2/6]')
call('mafft --quiet --auto --thread ' + sys.argv[2] + ' temp/queryCOG.fasta > temp/queryCOG.fasta.aln', shell = True, stdout=open(os.devnull, 'wb'))
call('hmmbuild temp/queryCOG.hmm temp/queryCOG.fasta.aln', shell = True, stdout=open(os.devnull, 'wb'))
call('hmmsearch -E 1e-5 --incE 1e-5 --domE 1e-5 --cpu ' + sys.argv[2] + ' -o /dev/null --tblout temp/queryCOG.hmmhits temp/queryCOG.hmm files/eukaryotes.fasta', shell = True, stdout=open(os.devnull, 'wb'))
call('faSomeRecords files/eukaryotes.fasta temp/queryCOG.hmmhits temp/queryCOG.hmmhits.fasta', shell = True, stdout=open(os.devnull, 'wb'))
print('\nGenerating phylogeny of the HMM hits [3/6]')
call('mafft --quiet --auto --thread ' + sys.argv[2] + ' temp/queryCOG.hmmhits.fasta > temp/queryCOG.hmmhits.fasta.aln', shell = True, stdout=open(os.devnull, 'wb'))
call('trimal -in temp/queryCOG.hmmhits.fasta.aln -out temp/queryCOG.hmmhits.fasta.aln.trimal -gt 0.3', shell = True, stdout=open(os.devnull, 'wb'))
call('iqtree -m LG+F+R5 -fast -s temp/queryCOG.hmmhits.fasta.aln.trimal -alrt 1000 -nt ' + sys.argv[2], shell = True, stdout=open(os.devnull, 'wb'))

# annotate phylogeny using BLAST

print('\nAnnotating sequences [4/6]')
call('diamond blastp --quiet --sensitive --query temp/queryCOG.hmmhits.fasta --db files/uniprot_sprot.90.fasta.dmnd --outfmt 6 --evalue 1e-20 --out temp/queryCOG.sprot.blastp --max-target-seqs 100 --max-hsps 1 --threads ' + sys.argv[2], shell = True, stdout=open(os.devnull, 'wb'))
blastfile = open('temp/queryCOG.sprot.blastp','r').readlines()
blast_annotation = {}
seen = []
for l in blastfile:
    if l.split('\t')[0] not in seen:
        blast_annotation[l.split('\t')[0]] = l.split('\t')[1].split('|')[2].split('_')[0]
        seen.append(l.split('\t')[0])

allseqs = []
hmmhits = open('temp/queryCOG.hmmhits.fasta','r').readlines()
for line in hmmhits:
    if line.startswith('>'):
        allseqs.append(line.strip('>').strip())

seq_annotation = {}
for s in allseqs:
    s = s.strip()
    try:
        seq_annotation[s] = blast_annotation[s]
    except:
        seq_annotation[s] = 'NoHit'
query_annot = seq_annotation[mapped_seq]

# load phylogeny into ETE3 and identify the orthologue clade

# if there are not hits to swissprot, just use the POG for presence/absence
print('\nIdentifying orthologues [5/6]')

if query_annot == 'NoHit':
    print('\nProtein identity unknown. Resorting to using original COG.\n')
    call('cp temp/queryCOG.fasta temp/orthologues.fasta', shell = True, stdout=open(os.devnull, 'wb'))
else:
    # for comparing uniprot gene names - use a similarity measure to account for subtle variation in gene identifiers   
    def similar(a, b):
        return SequenceMatcher(None, a, b).ratio()
    # for annotating leafs with their annotations
    def annotate_leaf_proteins(t):         
         for leaf in t.get_leaves():
             leaf.add_feature("annot", protein_annotations[leaf.name])
    annotation_similarity_threshold = 0.8
    a = 0
    x = 0
    while a == 0:
        protein_annotations = {}
        for s in seq_annotation:
            if seq_annotation[s] == query_annot:
                protein_annotations[s] = 'a'
                x += 1
            elif similar(query_annot, seq_annotation[s]) >= annotation_similarity_threshold:
                protein_annotations[s] = 'a'
                x += 1
            elif seq_annotation[s] == 'NoHit':
                protein_annotations[s] = 'b'
            else:
                protein_annotations[s] = 'c'
        # load in the tree         
        t = PhyloTree('temp/queryCOG.hmmhits.fasta.aln.trimal.treefile')
        # try to root with largest paralog clade
        try:
            m = 0
            for clade in t.get_monophyletic(values=["c","b"], target_attr="annot"):
                if len([l for l in clade.get_leaves()]) > m:
                    rooting_clade = clade
                    m = len([l for l in clade.get_leaves()])
            t.set_outgroup(clade)
        except:  
            midpoint = t.get_midpoint_outgroup()
            if midpoint: t.set_outgroup(midpoint)
        annotate_leaf_proteins(t)
        # find the clade with the mapped sequence
        for l in t.iter_descendants():
            if l.name == mapped_seq:
                clade = l
        n = 0
        while n == 0:
            if clade.up:
                annotations = [leaf.annot for sister in clade.get_sisters() for leaf in sister.get_leaves()]
                counts = dict(Counter(annotations))
                if 'c' not in counts.keys():
                    counts['c'] = 0
                if ('a' in annotations) and (counts['a']/(counts['a']+counts['c']) > 0.15):
                    clade = clade.up
                elif 'c' not in annotations:
                    clade = clade.up
                elif clade.up.up:
                    up_annotations = [leaf.annot for sister in clade.up.get_sisters() for leaf in sister.get_leaves()]
                    if ('a' in annotations) and (counts['a']/(counts['a']+counts['c']) > 0.15):
                        clade = clade.up
                    else:
                        n += 1
                else:
                    n += 1
            else:
                n += 1
        # if the identified proteins are too few, try decreasing the gene name similarity threshold
        seqs = [leaf.name for leaf in clade.get_leaves()]
        # remove large clades of paralogs if present
        for paralog in clade.get_monophyletic(values=["c"], target_attr="annot"):
            if len([leaf.name for leaf in paralog.get_leaves()]) >= 5:
                for leaf in paralog.get_leaves():
                    seqs.remove(leaf.name) 
        # if we only get 50% of sequences with the query annotation, try reducing name similarity threshold
        if (len(seqs) <= 0.5*x) and (annotation_similarity_threshold == 0.8):
            annotation_similarity_threshold = 0.66
        elif (len(seqs) <= 0.5*x) and (annotation_similarity_threshold == 0.66):
            annotation_similarity_threshold = 0.5
        elif(len(seqs) <= 0.5*x) and (annotation_similarity_threshold == 0.5):
            a += 1
        else:
            a += 1
    # write out a file with the identified proteins
    out = open('temp/treeseqs','w')
    sep = '\n'
    out.write(sep.join(seqs))
    out.close()

# extract the final identified proteins
call('mkdir results', shell = True, stdout=open(os.devnull, 'wb'))
call('faSomeRecords files/eukaryotes.fasta temp/treeseqs results/' + query_annot +'.fasta', shell = True, stdout=open(os.devnull, 'wb'))

# generate a clade file using NCBI taxonomy and the orthologues identified previously

print('\nMaking clade file [6/6]')

# create NCBI taxonomy tree

# make NCBI taxonomy tree with all species  
taxa = []  
db = open('files/eukaryotes.fasta','r').readlines()  
for line in db: 
    if line.startswith('>'): 
        taxa.append(line.split('.')[0].strip('>')) 
taxa = list(set(taxa)) 
tree = ncbi.get_topology(taxa)   

# load in NCBI tree and annotate supergroups and protein presence

def annotate_supergroup(t):         
    for leaf in t.get_leaves():
        leaf.add_feature("group", ncbi.get_lineage(int(leaf.name))[3])

def protein_presence(t,f):
    seqs = open(f,'r').readlines()
    taxa = [t.split('.')[0] for t in seqs]
    for leaf in t.get_leaves():
        if leaf.name in taxa:
            leaf.add_feature("presence", "1")
        else:
            leaf.add_feature("presence", "0")


annotate_supergroup(tree)
protein_presence(tree,'temp/treeseqs')
supergroups = list(set([leaf.group for leaf in tree.get_leaves()]))

out = open('temp/cladefile.tab','w')
sep = '\t'

absents = 0
for sg in supergroups:
    for clade in tree.get_monophyletic(values=[sg], target_attr="group"):
        taxa = [leaf.name for leaf in clade.get_leaves()]
        absentees = []
        for subclade in clade.get_monophyletic(values=["0"], target_attr="presence"):
            out.write('absent\t'+sep.join([leaf.name for leaf in subclade.get_leaves()])+'\n')
            absents += 1
            absentees = absentees + [leaf.name for leaf in subclade.get_leaves()]
        for t in absentees:
            taxa.remove(t)
        if len(taxa) > 0:
            out.write('present\t'+sep.join(taxa)+'\n')
out.close()

# Run the correlation analysis

if absents > 0:
    call('python files/orthogroup_enrichment_V2.py files/dump.out.eukaryotes.allVSall.dmnd.blastp.mci.I14.filt temp/cladefile.tab ' + sys.argv[3], shell = True)     
else:
   print('\n Error: Your gene is conserved across all species. Try another gene or dataset.')
   
# Assign metadata to correlated clusters

# record metadata
uniprot = open('files/UniProt_ids.metadata.tab','r').readlines()
accession_d = {}

sep = ','
for line in uniprot:
    accession = line.split('\t')[0]
    gene = line.split('\t')[2].split('(')[0].strip().replace(' ','_')
    try:
        loc = [l.split('[')[0].strip() for l in line.split('\t')[3].split(';')]
    except:
        loc = ['NA']
    accession_d[accession] = [gene, sep.join(loc)]
    
# load in the clusters
clusters = open('correlated_clusters.tab','r').readlines()
out = open('correlated_clusters_annotated.tab','w')
out.write(clusters[0].strip()+'\tGene_name\tLocalization\tHuman_presence\tHuman_sequences\tArabidopsis_presence\tArabidopsis_sequences\n')
sep = ', '
for c in clusters[1:]:
    try:
        cluster = open(c.split('\t')[0],'r').readlines()
    except:
        continue
    # record key species
    hs = []
    at = []
    for seq in cluster:
        if seq.split('.')[0] == '3702':
            at.append(seq.strip())
        elif seq.split('.')[0] == '9606':
            hs.append(seq.strip())
    hs_presence = 'no'
    at_presence = 'no'
    if len(hs) > 0:
        hs_presence = 'yes'
    else:
        hs = ['NA']
    if len(at) > 0:
        at_presence = 'yes'
    else:
        at = ['NA']
    # record gene names using arabidopsis and human
    gene_names = []
    for seq in cluster:
        if (seq.split('.')[0] == '3702') or (seq.split('.')[0] == '9606'):
            try:
                gene_names.append(accession_d[seq.split('.')[1].strip()][0])
            except:
                pass
    gene_names = [i[0] for i in Counter(gene_names).most_common(5)]
    if len(gene_names) == 0:
        gene_names = ['NA']
    # record localizations using arabidopsis and human
    localizations = []
    for seq in cluster:
        if (seq.split('.')[0] == '3702') or (seq.split('.')[0] == '9606'):
            try:
                if accession_d[seq.split('.')[1].strip()][1] != 'NA':
                    localizations = localizations + accession_d[seq.split('.')[1].strip()][1].split(',')
            except:
                pass
    localizations = [i[0] for i in Counter(localizations).most_common(5)]
    if len(localizations) == 0:
        localizations = ['NA']
    # output results
    out.write(c.strip()+'\t'+sep.join(gene_names).strip(',')+'\t'+sep.join(localizations).strip(',')+'\t'+hs_presence+'\t'+sep.join(hs[0:11])+'\t'+at_presence+'\t'+sep.join(at[0:11])+'\n')
out.close()