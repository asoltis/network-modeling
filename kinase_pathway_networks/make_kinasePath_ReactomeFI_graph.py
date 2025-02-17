import sys, os, cPickle
import networkx as nx

## Data sources
# kinase-substrate links (gmt format)
ks_fn = 'PhosphoSitePlus/parsed_files/kinase_substrate_sets.human_to_human.20190304.gmt'

# protein-pathway links (gmt format)
path_fn = 'MSigDB/downloads/Hallmark_sets/h.all.v7.0.symbols.gmt'

# ReactomeFI PPI
ppi_fn = 'ReactomeFI/version_2017/parsed_versions/symbol_update/ReactomeFI_071718.all_data.no_UBC.HGNC_symbol_update_20190524.txt'

# TF list
tf_fn = 'TFs/human_TFs/Lambert_etal_human_TFs/parsed_files/Human_TFs_HGNC_symbols.txt'
TFs = [x.strip() for x in open(tf_fn)]

# whether or not to filter interactome edges labelled only as 'predicted'
no_predict = True

## Create graph
G = nx.DiGraph()

# add kinase-substrate links and substrate-protein links
kinases = []
nodeAttributes = {}
edgeAttributes = []
for line in open(ks_fn).readlines():
    l = line.strip().split('\t')
    k = l[0]
    kinases.append(k)
    for subs in l[2:]:
        if subs.startswith('_'): continue # Issue with no gene symbol for substrate
        edge = (k, subs)
        G.add_edge(*edge, etype = 'PSP_KS', score = 0.99, directed = 'D', direction = '->')
        ea = [k, subs, 'D', '->', 0.99, 'PSP_KS']
        edgeAttributes.append(ea)

        prot = subs.split('_')[0]
        edge2 = (subs, prot)
        G.add_edge(*edge2, etype = 'SubsProt', score = 0.99, directed = 'D', direction = '->')
        ea2 = [subs, prot, 'D', '->', 0.99, 'SubsProt']
        edgeAttributes.append(ea2)

        nodeAttributes[k] = 'Kinase'
        nodeAttributes[subs] = 'Phospho'
        if prot not in nodeAttributes:
            if prot in TFs:
                nodeAttributes[prot] = 'TransFac'
            else:
                nodeAttributes[prot] = 'Protein'

# add protein to pathway links
pathways = []
for line in open(path_fn).readlines():
    l = line.strip().split('\t')
    path = l[0]
    pathways.append(path)
    for prot in l[2:]:
        edge = (prot, path)
        G.add_edge(*edge, etype = 'ProtPath', score = 0.99, directed = 'D', direction = '->')
        ea = [prot, path, 'D', '->', 0.99, 'ProtPath']
        edgeAttributes.append(ea)

        if prot not in nodeAttributes:
            if prot in TFs:
                nodeAttributes[prot] = 'TransFac'
            else:
                nodeAttributes[prot] = 'Protein'
        nodeAttributes[path] = 'Pathway'

# PPI links
DSET = set()
for line in open(ppi_fn):
    if line.startswith('Gene1'): continue
    ga, gb, annotation, direction, score = line.strip().split('\t')
    score = float(score)
    if ga == 'UBC' or gb == 'UBC': continue # skip ubiquitin

    ga = ga.replace(' ','_')
    gb = gb.replace(' ','_')

    if no_predict:
        if annotation == 'predicted': continue

    #if 'expression' in annotation:
    #    if direction == '-':
    #        print ga, gb, annotation, direction
    #    DSET.add(direction)

    if direction == '-': # non-directed edge
        edge1 = (ga, gb)
        edge2 = (gb, ga)
        ea1 = [ga, gb, 'U', direction, score, annotation]
        ea2 = [gb, ga, 'U', direction, score, annotation]
        G.add_edge(*edge1, etype = annotation, score = score, directed = 'U', direction = direction)
        G.add_edge(*edge2, etype = annotation, score = score, directed = 'U', direction = direction)
        edgeAttributes.append(ea1)
        edgeAttributes.append(ea2)

    elif direction == '<-' or direction == '|-':
        if direction == '<-': direction = '->'
        elif direction == '|-': direction = '-|'
        edge = (gb, ga)
        ea = [gb, ga, 'D', direction, score, annotation]
        G.add_edge(*edge, etype = annotation, score = score, directed = 'D', direction = direction)
        edgeAttributes.append(ea)
         
    elif direction == '->' or direction == '-|':
        edge = (ga, gb)
        ea = [ga, gb, 'D', direction, score, annotation]
        G.add_edge(*edge, etype = annotation, score = score, directed = 'D', direction = direction)
        edgeAttributes.append(ea)
   
    elif direction in ['|-|', '<-|', '|->', '<->']:
        dir1 = direction[1:]
        edge1 = (ga, gb)
        ea1 = [ga, gb, 'D', dir1, score, annotation]
        G.add_edge(*edge1, etype = annotation, score = score, directed = 'D', direction = dir1)
        edgeAttributes.append(ea1)

        dir2 = direction[0:2][::-1]
        if dir2 == '-<': dir2 = '->'
        edge2 = (gb, ga)
        ea2 = [gb, ga, 'D', dir2, score, annotation]
        G.add_edge(*edge2, etype = annotation, score = score, directed = 'D', direction = dir2)
        edgeAttributes.append(ea2)

    if ga not in nodeAttributes:
        if ga in TFs:
            nodeAttributes[ga] = 'TransFac'
        else:
            nodeAttributes[ga] = 'Protein'
    if gb not in nodeAttributes:
        if gb in TFs:
            nodeAttributes[gb] = 'TransFac'
        else:
            nodeAttributes[gb] = 'Protein'

# Write out graph
if no_predict:
    odir = '../graphs/ReactomeFI_noPredicted_PSP-KS_MSigDB-v7-Hallmark/'
else:
    odir = '../graphs/ReactomeFI_PSP-KS_MSigDB-v7-Hallmark/'
os.system('mkdir -p %s'%(odir))
ofbase = odir + 'ReactomeFI_PSP-KS_MSigDB-v7-Hallmark'
ofpkl = open(ofbase + '.pkl','w')
cPickle.dump(G, ofpkl)

# write out PCSF format interactome
ofpcsf = open(ofbase + '_for_PCSF.txt','w')
Gcheck = nx.Graph()
fmtstr = '%s\t%s\t%0.2f\t%s\n'
for e in G.edges():
    score = G[e[0]][e[1]]['score']
    directed = G[e[0]][e[1]]['directed']
    if directed == 'D':
        ol = fmtstr % (e[0], e[1], score, directed)
    elif directed == 'U':
        edge = (e[0], e[1])
        if edge in Gcheck.edges(): 
            #print edge
            continue
        ol = fmtstr % (e[0], e[1], score, directed)
        Gcheck.add_edge(*edge)
    ofpcsf.writelines(ol)
ofpcsf.close()

# write out node attributes
ofattr = open(ofbase + '_nodeAttributes.txt','w')
ofattr.writelines('Node\tType\n')
for n in nodeAttributes:
    ofattr.writelines('%s\t%s\n'%(n, nodeAttributes[n]))
ofattr.close()

# write edge attributes
ofea = open(ofbase + '_edgeAttributes.txt', 'w')
ofea.writelines('Edge\tScore\tDirection\tAnnotation\n')
for edata in edgeAttributes:
    n1, n2, dflag, direction, score, annotation = edata
    if dflag == 'D':
        oea = '%s (pd) %s\t%0.2f\t%s\t%s\n'%(n1, n2, score, direction, annotation)
    elif dflag == 'U':
        oea = '%s (pp) %s\t%0.2f\t%s\t%s\n'%(n1, n2, score, direction, annotation)

    ofea.writelines(oea)
ofea.close()

