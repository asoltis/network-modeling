import sys, os, csv, _pickle
import networkx as nx

Gname = 'ReactomeFI_2020_remove_UBC_signed_exclude_transcriptional_edges_undirected_graph'
Gin = _pickl.load(open('', 'rb'))
nodes0 = list(Gin.nodes)
Goname = Gname + '_CCLE_metab_PMI_RPPA_phosphosites.pkl'

# CCLE metabolomic info to add
metabInfo = 'data_files/metabolite_HMDB_conversion_with_protein_interactions.txt'

with open(metabInfo) as MI:
    reader = csv.reader(MI, delimiter = '\t')
    for ri, row in enumerate(reader):
        if ri == 0: continue
        hid = row[1]
        if hid == 'NA': continue
        prots = row[-1].split(';')

        for p in prots:
            if p == '': continue
            Gin.add_edge(hid, p, weight = 1.0, interaction = '-', annotation = 'protein-metabolite')
            if Gin.is_directed():
                Gin.add_edge(p, hid, weight = 1.0, interaction = '-', annotation = 'protein-metabolite')

# CCLE RPPA antibodies
phos_sites = [x.strip() for x in open('CCLE_RPPA_phosphosites.txt')]
for ps in phos_sites:
    g, site = ps.split('_')
    Gin.add_edge(ps, g, weight = 1.0, interaction = '->', annotation = 'phosphosite-protein')
    if Gin.is_directed():
        Gin.add_edge(g, ps, weight = 1.0, interaction = '->', annotation = 'phosphosite-protein')


print('Total nodes: %d' % (len(Gin.nodes)))
print('Total undirected edges: %d' % (len(Gin.edges)))

_pickle.dump(Gin, open(Goname, 'wb'))


