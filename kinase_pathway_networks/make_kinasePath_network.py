import sys, os, cPickle
import networkx as nx

# Inputs
args = sys.argv[1:]
netfn = args[0]
termfn = args[1]
ofbase = args[2]
deprot_fn = None
if len(args) > 3:
    deprot_fn = args[3]

kin_list_fn = 'kinases/human_protein_kinases_pkinfam_20190605.txt'
kin_list = [x.strip() for x in open(kin_list_fn)]
tf_list_fn = 'TFs/human_TFs/Lambert_etal_human_TFs/parsed_files/Human_TFs_HGNC_symbols.txt'
tf_list = [x.strip() for x in open(tf_list_fn)]

# Load graph
G = cPickle.load(open(netfn))
gattr_fn = netfn.replace('.pkl','_nodeAttributes.txt')
Gattr = {}
for line in open(gattr_fn).readlines()[1:]:
    n, attr = line.strip().split('\t')
    Gattr[n] = attr

# Load terminals
terminals = set()
kinases, phosphos, TFs, mutations, pathways = set(), set(), set(), set(), set()
terminalInfo = {}
for line in open(termfn).readlines()[1:]:
    l = line.strip().split('\t')
    term, typ = l[0:2]
    vals = [float(x) for x in l[3].split(';')]
    prize = float(l[4])

    val = None
    for v in vals:
        if abs(v) == prize: val = v

    if 'Kinase' in typ: 
        kinases.add(term)
        terminalInfo[term] = ['Kinase',val]
    elif 'Pathway' in typ: 
        pathways.add(term)
        terminalInfo[term] = ['Pathway',val]
    elif 'phos' in typ:
        phosphos.add(term)
        terminalInfo[term] = ['Phospho',val]
    elif 'TransFac' in typ:
        TFs.add(term)
        terminalInfo[term] = ['TransFac',val]
    elif 'Mutation' in typ:
        mutations.add(term)
        terminalInfo[term] = ['Mutation',val]
    else:
        if abs(val) < 0.25: continue ## effect size filter on proteins (DEFAULT: 0.5)
        attr = 'Protein'
        if term in Gattr:
            gattr = Gattr[term]
            if gattr != 'Protein':
                attr = 'Protein-%s'%(gattr)
        terminalInfo[term] = [attr, val]

    terminals.add(term)

deprots = {}
if deprot_fn != None:
    for line in open(deprot_fn).readlines()[1:]:
        l = line.strip().split('\t')
        prot, beta, pfdr = l[1], float(l[2]), float(l[4])
        deprots[prot] = [beta, pfdr]

# annotation terms in PPI to skip
skip_annots = ['expression regulated by', 'expression regulates']

##
GOUT = nx.DiGraph()
for k in kinases:
    #if terminalInfo[k][1] < 0: continue # take only positive
    for phos in phosphos:
        if (k, phos) in G.edges:
            prot = phos.split('_')[0]
            #if prot in terminals: # add all kinase->phosphosite->protein edges if all in terminals
            # add all enriched kinase -> enriched phosphosite -> protein edges
            kpedge = (k, phos)
            ppedge = (phos, prot)
            GOUT.add_edge(*kpedge)
            GOUT.add_edge(*ppedge)

            for tf in TFs:
                #if terminalInfo[tf][1] < 0: continue # take only positive
                if tf not in G.nodes: continue
                
                #for pw in pathways:
                #    if (tf, pw) not in G.edges: 
                #        continue

                paths = nx.all_simple_paths(G, source = prot, target = tf, cutoff = 2)
                pscores = []

                try:
                    
                    for path in paths:
                        inTerm = []
                        for p in path[1:-1]:
                            if p in terminals: inTerm.append(1)
                            else: inTerm.append(0)
                        if sum(inTerm) == len(inTerm): # - 1:
                            skip = False
                            for i in range(0, len(path)-1):
                                annot = G[path[i]][path[i+1]]['etype']
                                if annot in skip_annots: skip = True
                            for i in range(0, len(path)):
                                n = path[i]
                                if '_' in n and n not in terminals:
                                    skip = True ## Skip un-measured phosphosite paths

                            if not skip:
                                if 'expression' in annot:
                                    if not ('activate' in annot or 'catalyze' in annot or 'complex' in annot):
                                        continue
                                    #else:
                                    #    print path[i], path[i+1], annot
                                
                                path2 = [k, phos] + path
                                #for i in range(0, len(path2)-1):
                                #    edge = (path2[i], path2[i+1])
                                #    GOUT.add_edge(*edge)
                                pscore = 0
                                for n in path2:
                                    if n in terminalInfo:
                                        pscore += abs(terminalInfo[n][1])
                                pscore = pscore / len(path2)
                                pscores.append([path2, pscore])
                                
                except nx.NetworkXNoPath: continue

                #
                pscores.sort(key = lambda x: x[1], reverse = True)
                if len(pscores) == 0: continue
                for ps in pscores: #[0:5]:
                    path = ps[0]
                    print path
                    for i in range(0, len(path)-1):
                        edge = (path[i], path[i+1])
                        GOUT.add_edge(*edge)

#nodesC = list(GOUT.nodes)
for m in mutations:
    for node in list(kinases)+list(TFs):
        #if terminalInfo[node][1] < 0: continue # taking positive
        if node not in G.nodes: continue
        paths = nx.all_simple_paths(G, source = m, target = node, cutoff = 3)
        try:
            for path in paths:
                inTerm = []
                for p in path[1:-1]:
                    if p in terminals: inTerm.append(1)
                    else: inTerm.append(0)
                if sum(inTerm) == len(inTerm): # - 1:
                    skip = False
                    for i in range(0, len(path) - 1):
                        annot = G[path[i]][path[i+1]]['etype']
                        if annot in skip_annots: skip = True

                    if not skip:
                        if 'expression' in annot:
                            if not ('activate' in annot or 'catalyze' in annot or 'complex' in annot):
                                continue
                        pscore = 0
                        for i in range(0, len(path)-1):
                            pscore += G[path[i]][path[i+1]]['score']
                            edge = (path[i], path[i+1])
                            GOUT.add_edge(*edge)
                        print path

        except nx.NetworkXNoPath: continue

# Remove cycle edges resulting from kinase to substrate links from PPI network already accounted for with kinase->psite->substrate
EDGES = GOUT.edges()
#kinIn = [x for x in list(GOUT.nodes) if Gattr[x] == 'Kinase'] ## In future, should run with these kinases as extras can be added
for k in kinases:
    for phos in phosphos:
        prot = phos.split('_')[0]
       
        if (k, phos) in EDGES and ( (k, prot) in EDGES or (prot, k) in EDGES ):
            if (k, prot) in EDGES:
                edge = (k, prot)
                GOUT.remove_edge(*edge)
                print 'Removed: ',edge,G[k][prot]
            elif (prot, k) in EDGES:
                directed = G[prot][k]['directed']
                if directed == 'U': 
                    edge = (prot, k)
                    GOUT.remove_edge(*edge)
                    print 'Removed: ',edge, G[prot][k]

# Add pathway links
for term in terminals:
    for pw in pathways:
        edge = (term, pw)
        if edge in G.edges():
            GOUT.add_edge(*edge)
            print edge

# Add additional significant phosphosite links for nodes in network
NODES = list(GOUT.nodes)
for node in NODES:
    for ps in phosphos:
        prot = ps.split('_')[0]
        if prot == node:
            if (ps, prot) not in GOUT.edges():
                edge = (ps, prot)
                GOUT.add_edge(*edge)
                G.add_edge(*edge)
                G[ps][prot]['directed'] = 'D'
                G[ps][prot]['etype'] = 'SubsProt'
                G[ps][prot]['score'] = 0.99
                G[ps][prot]['direction'] = '->'
                print edge

# TESTING adding back edges
if False:
    NODES = list(GOUT.nodes)
    for i1 in range(0, len(NODES)-1):
        for i2 in range(i1 + 1, len(NODES)):
            n1 = NODES[i1]
            n2 = NODES[i2]
            if (n1, n2) not in GOUT.edges() and (n2, n1) not in GOUT.edges():
                if (n1, n2) in G.edges():
                    annot = G[n1][n2]['etype']
                    valid = True
                    if 'expression' in annot:
                        if not ('activate' in annot or 'catalyze' in annot or 'complex' in annot):
                            valid = False
                    if valid:
                        print n1, n2, G[n1][n2]
                elif (n2, n1) in G.edges():
                    annot = G[n2][n1]['etype']
                    valid = True
                    if 'expression' in annot:
                        if not ('activate' in annot or 'catalyze' in annot or 'complex' in annot):
                            valid = False
                    if valid:
                        print n2, n1, G[n2][n1]

#for node in list(GOUT.nodes):
#    ti = None
#    try: ti = terminalInfo[node]
#    except: pass
#    if node in kin_list:
#        print node, 'kinase', ti
#    elif node in tf_list:
#        print node, 'TF', ti

# write out
of = open(ofbase + '.sif', 'w')
ofattr = open(ofbase + '_node_attributes.txt','w')
ofattr.writelines('Node\tTerminal\tType\tValue\n')
ofeattr = open(ofbase + '_edge_attributes.txt', 'w')
ofeattr.writelines('Edge\tScore\tDirection\tAnnotation\n')
terminalOut = set()
for edge in GOUT.edges:
    # write out edge info to sif
    e1, e2 = edge
    directed = G[e1][e2]['directed']
    if directed == 'D': dirstr = 'pd'
    elif directed == 'U': dirstr = 'pp'
    of.writelines('%s\t%s\t%s\n'%(e1.replace('HALLMARK_',''), dirstr, e2.replace('HALLMARK_','')))
    
    # get edge attribute info
    edgestr = '%s (%s) %s' % (e1.replace('HALLMARK_',''), dirstr, e2.replace('HALLMARK_',''))
    escore = G[e1][e2]['score']
    eannot = G[e1][e2]['etype']
    edirection = G[e1][e2]['direction']
    ofeattr.writelines('%s\t%0.2f\t%s\t%s\n' % (edgestr, escore, edirection, eannot))

    # get node attribute info
    if e1 in terminalInfo:
        tout = '%s\t%s\t%s\t%0.5f\n'%(e1.replace('HALLMARK_',''), 'Terminal', terminalInfo[e1][0], terminalInfo[e1][1])
        terminalOut.add(tout)
    elif e1 in deprots:
        #tout = '%s\t%s\t%s\t%0.5f\n'%(e1.replace('HALLMARK_',''), 'Interactome', Gattr[e1], deprots[e1]) # original
        beta, pfdr = deprots[e1]
        ptype = 'Interactome'
        attr = Gattr[e1]
        if pfdr < 0.1: 
            ptype = 'Terminal'
            if attr != 'Protein': attr = 'Protein-%s'%(Gattr[e1])
        tout = '%s\t%s\t%s\t%0.5f\n'%(e1.replace('HALLMARK_',''), ptype, attr, beta)
        terminalOut.add(tout) 
    else:
        tout = '%s\t%s\t%s\tNA\n'%(e1.replace('HALLMARK_',''), 'Interactome', Gattr[e1])
        terminalOut.add(tout)
    if e2 in terminalInfo:
        tout = '%s\t%s\t%s\t%0.5f\n'%(e2.replace('HALLMARK_',''), 'Terminal', terminalInfo[e2][0], terminalInfo[e2][1])
        terminalOut.add(tout)
    elif e2 in deprots:
        #tout = '%s\t%s\t%s\t%0.5f\n'%(e2.replace('HALLMARK_',''), 'Interactome', Gattr[e2], deprots[e2]) # original
        beta, pfdr = deprots[e2]
        ptype = 'Interactome'
        attr = Gattr[e2]
        if pfdr < 0.1: 
            ptype = 'Terminal'
            if attr != 'Protein': attr = 'Protein-%s'%(Gattr[e2])
        tout = '%s\t%s\t%s\t%0.5f\n'%(e2.replace('HALLMARK_',''), ptype, attr, beta)
        terminalOut.add(tout)
    else:
        tout = '%s\t%s\t%s\tNA\n'%(e2.replace('HALLMARK_',''), 'Interactome', Gattr[e2])
        terminalOut.add(tout)
for tout in terminalOut:
    ofattr.writelines(tout)
of.close()
ofattr.close()
ofeattr.close()

