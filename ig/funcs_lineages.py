import os, pandas as pd, numpy as np, atriegc

#read example igBlast-aligned repertoire
filename = 'healthy/BZR10_117_t0_IgG_1_HC_unique.df'
usecols = ['seq_ID',
           'igBlast_V_best_align_family',
           'igBlast_J_best_align_family',
           'igBlast_CDR3_identity',
           'noHyperIndels_seq_nt']
df = pd.read_table(filename, usecols=usecols)
df['v_gene'] = df['igBlast_V_best_align_family']
df['j_gene' = df['igBlast_J_best_align_family']
df['cdr3'] = df['igBlast_CDR3_identity']
df['cdr3_length'] = df['cdr3'].str.len()

#infer lineage structure
threshold_l = 29
threshold_x = 0.1
by = ['v_gene','j_gene','cdr3_length']
gdfs = []
for (v,j,l),group in df.loc[df['cdr3_length']>threshold_l].groupby(by):
    gdf = group.copy()
    trie = atriegc()
    for cdr3 in gdf['cdr3']:
        trie.insert(cdr3)
    t = int(threshold_x*l)
    dct = trie.clusters(t)
    gdf['cluster'] = gdf['cdr3'].map(dct)
    gdfs.append(gdf)

new_df = pd.concat(gdfs)
new_df['clone'] = new_df.groupby(by + ['cluster']).ngroup()+1

#distribution of lineage sizes
#use large lineages to reconstruct phylogenies
sizes = new_df.groupby(['clone']).size()
sizes_bins = np.arange(1,int(1e4)+1,1)
sizes_histogram,_ = np.histogram(sizes,bins=bins,density=True)
threshold_size = 50
#indices of lineages to be analyzed with RAxML:
indices = sizes[sizes>=threshold_size].index
#phylogeny and ancestral states reconstruction under GTRGAMMA
#raxml commands $clone is the lineage index
#raxml -s "lineage_$clone.fasta" -w healthy/ -g "start_tree_$clone.newick" -o "$clone*0" -m GTRGAMMA -n "$clone" -p 12345 
#raxml -f A -s "lineage_$clone.fasta" -w healthy/ -t RAxML_bestTree.$clone" -m GTRGAMMA -n "$1a" -p 12345 
   
#analyze phylogenies in directories
#newick files in directory correspond to a given cohort
#estimate skewedness within cohorts
from ete3 import Tree
def get_weights(dirname):
    weights = []
    for filename in os.listdir(dirname):
        if 'RAxML_bestTree' in filename:
            tree = open(os.path.join(dirname,filename), 'r').read()
            t = Tree(tree,format=1)
            t.set_outgroup('0')
            for c in t.children:
                if not c.name=='0':
                    for c2 in c.children:
                        weights.append(len(c2.get_leaves())/(len(t)-1))
    return weights
   
def get_fractions(dirname):
    fractions = []
    for filename in os.listdir(dirname):
        if 'RAxML_bestTree' in filename:
            tree = open(os.path.join(dirname,filename), 'r').read()
            t = Tree(tree,format=1)
            t.set_outgroup('0')
            for node in list(t.search_nodes(dist=0)):
                if node.name!='0':
                    parent = node.up
                    if parent:
                        for c in list(parent.get_children()):
                            c.detach()
                            if (c == node):
                                for d in c.get_children():
                                    parent.add_child(d)
                            else: parent.add_child(c)  

            branch_lengths,branch_lengths_terminal = [],[]
            for node in t.traverse('preorder'):
                if not node.name=='0':
                    if node.is_leaf():
                        branch_lengths_terminal.append(node.dist)
                    branch_lengths.append(node.dist)
            fractions.append(np.mean(branch_lengths_terminal)/np.mean(branch_lengths))
    return fractions

weights_healthy = get_weights('healthy/')
weight_unit = 0.05
weights_bins = np.arange(0,1+weight_unit,weight_unit)
weights_histogram,_ = np.histogram(weights_healthy,bins=weights_bins,density=True)
fractions_healthy = get_fractions('healthy/')
fractions_unit = 0.03
fractions_bins = np.arange(0,2+fractions_unit,fractions_unit)
fractions_histogram,_ = np.histogram(fractions_healthy,bins=fractions_bins,density=True)
