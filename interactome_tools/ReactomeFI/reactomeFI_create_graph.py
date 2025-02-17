import sys, os, csv, _pickle
import numpy as np
import networkx as nx

# Input file
infile = 'FIsInGene_070323_with_annotations.txt'
baseo = 'ReactomeFI_2022'

# Transcription factors list
TFlist = [x.strip() for x in open('data_files/TF_names_v_1.01.txt')]
regulators = [x.strip() for x in open('data_files/GO_transcriptional_regulators.txt')]

# Filters
remove_nodes = ['UBC']
pred_scoreT = 0.95 # minimum score for predicted interactions
filt_predict = False
signed = True

# Tail to append to gene transcripts
gtail = '-gene'

# output file base
if filt_predict:
    baseo += '_remove_%s_predScore_ge_%s' % ('_'.join(remove_nodes), str(pred_scoreT))
else:
    baseo += '_remove_%s' % ('_'.join(remove_nodes))
if signed:
    baseo += '_signed'

# Read file and create directed graph
dtypes = set()
DG = nx.DiGraph() # main directed graph with all interactions
DGtf = nx.DiGraph() # directed graph with transcriptional interactions modified
DGnt = nx.DiGraph() # directed graph with transcriptional interactions excluded
with open(infile) as IF:
    reader = csv.reader(IF, delimiter = "\t")
    for ri, row in enumerate(reader):
        if ri == 0: continue # skip header
        n1, n2, anno, direction, score = row
        n1 = n1.replace(' ', '_')
        n2 = n2.replace(' ', '_')
        score = float(score)
        dtypes.add(direction)

        if n1 in remove_nodes or n2 in remove_nodes: continue
        if filt_predict:
            if anno == 'predicted' and score < pred_scoreT: continue

        expInteraction = False
        # expression regulation reactions
        isTF1, isTF2 = n1 in TFlist, n2 in TFlist
        isReg1, isReg2 = n1 in regulators, n2 in regulators
        if 'expression' in anno:
            expInteraction = True
            if ';' not in anno: # only expression regulation is reported
                if len(direction)==1 or len(direction)==3: # should not be any of these in this case
                    print(row)
                if isTF1 and (direction == '->' or direction == '-|'):
                    if signed and direction == '-|': sc = -score
                    else: sc = score
                    DGtf.add_edge(n1, n2+gtail, weight = sc, interaction = direction,
                                  annotation = anno, etype = 'transcription')
                elif isTF2 and (direction == '<-' or direction == '|-'):
                    if direction == '<-': dflip = '->'
                    elif direction == '|-': dflip = '-|'
                    if signed and dflip == '-|': sc = -score
                    else: sc = score
                    DGtf.add_edge(n2, n1+gtail, weight = sc, interaction = dflip,
                                  annotation = anno, etype = 'transcription')
                else: # regulation by chromatin modifiers, others generally
                    if (direction == '->' or direction == '-|'):
                        if signed and direction == '-|': sc = -score
                        else: sc = score

                        if isReg1:
                            DGtf.add_edge(n1, n2+gtail, weight = sc, interaction = direction,
                                          annotation = anno, etype = 'transcription')
                        else:
                            DGtf.add_edge(n1, n2, weight = sc, interaction = direction,
                                           annotation = anno, etype = 'expression')
                            DGnt.add_edge(n1, n2, weight = sc, interaction = direction,
                                           annotation = anno, etype = 'expression')
                    elif (direction == '<-' or direction == '|-'):
                        if direction == '<-': dflip = '->'
                        elif direction == '|-': dflip = '-|'
                        if signed and dflip == '-|': sc = -score
                        else: sc = score

                        if isReg2:
                            DGtf.add_edge(n2, n1+gtail, weight = sc, interaction = dflip,
                                          annotation = anno, etype = 'transcription')
                        else:
                            DGtf.add_edge(n2, n1, weight = sc, interaction = dflip,
                                           annotation = anno, etype = 'expression')
                            DGnt.add_edge(n2, n1, weight = sc, interaction = dflip,
                                           annotation = anno, etype = 'expression')

            # multiple annotations
            else:
                non_exp = []
                for an in anno.split('; '):
                    if 'expression' not in an and 'expressed' not in an: non_exp.append(an)

                if len(direction)==2:
                    if len(non_exp)==0:
                        if (direction == '->' or direction == '-|'):
                            if signed and direction == '-|': sc = -score
                            else: sc = score

                            if isTF1 or isReg1:
                                DGtf.add_edge(n1, n2+gtail, weight = sc, interaction = direction,
                                           annotation = anno, etype = 'transcription')
                            else:
                                DGtf.add_edge(n1, n2, weight = sc, interaction = direction,
                                              annotation = anno, etype = 'expression')
                                DGnt.add_edge(n1, n2, weight = sc, interaction = direction,
                                              annotation = anno, etype = 'expression')

                        elif (direction == '<-' or direction == '|-'):
                            if direction == '<-': dflip = '->'
                            elif direction == '|-': dflip = '-|'
                            if signed and dflip == '-|': sc = -score
                            else: sc = score

                            if isTF2 or isReg2:
                                DGtf.add_edge(n2, n1+gtail, weight = sc, interaction = dflip,
                                              annotation = anno, etype = 'transcription')
                            else:
                                DGtf.add_edge(n2, n1, weight = sc, interaction = dflip,
                                              annotation = anno, etype = 'expression')
                                DGnt.add_edge(n2, n1, weight = sc, interaction = dflip,
                                              annotation = anno, etype = 'expression')

                    else: # generally these are complex/activation interactions, often between two TFs
                        if direction == '->' or direction == '-|':
                            if signed and direction == '-|': sc = -score
                            else: sc = score
                            DGtf.add_edge(n1, n2, weight = sc, interaction = direction,
                                          annotation = anno, etype = 'interaction')
                            DGnt.add_edge(n1, n2, weight = sc, interaction = direction,
                                          annotation = anno, etype = 'interaction')
                        elif direction == '<-' or direction == '|-':
                            if direction == '<-': dflip = '->'
                            elif direction == '|-': dflip = '-|'
                            if signed and dflip == '-|': sc = -score
                            else: sc = score
                            DGtf.add_edge(n2, n1, weight = sc, interaction = direction,
                                          annotation = anno, etype = 'interaction')
                            DGnt.add_edge(n2, n1, weight = sc, interaction = direction,
                                          annotation = anno, etype = 'interaction')

                # direction 1 or 3 edges
                else:
                    if len(direction) == 1:
                        d1 = '-'
                        d2 = '-'
                    elif len(direction) == 3:
                        d1 = direction[1:]
                        d2 = direction[0:2]
                    s1, s2 = score, score
                    if d2 == '<-': d2 = '->'
                    elif d2 == '|-': d2 = '-|'
                    if signed and d1 == '-|': s1 = -s1
                    if signed and d2 == '-|': s2 = -s2

                    if len(non_exp) == 0:
                        if isTF1 or isReg1:
                            DGtf.add_edge(n1, n2+gtail, weight = s1, interaction = d1,
                                          annotation = anno, etype = 'transcription')
                        else:
                            DGtf.add_edge(n1, n2, weight = s1, interaction = d1,
                                          annotation = anno, etype = 'expression')
                            DGnt.add_edge(n1, n2, weight = s1, interaction = d1,
                                          annotation = anno, etype = 'expression')

                        if isTF2 or isReg2:
                            DGtf.add_edge(n2, n1+gtail, weight = s2, interaction = d2,
                                          annotation = anno, etype = 'transcription')
                        else:
                             DGtf.add_edge(n2, n1, weight = s2, interaction = d2,
                                           annotation = anno, etype = 'expression')
                             DGnt.add_edge(n2, n1, weight = s2, interaction = d2,
                                           annotation = anno, etype = 'expression')

                    else:
                        DGtf.add_edge(n1, n2, weight = s1, interaction = d1,
                                      annotation = anno, etype = 'interaction')
                        DGtf.add_edge(n2, n1, weight = s2, interaction = d2,
                                      annotation = anno, etype = 'interaction')

                        DGnt.add_edge(n1, n2, weight = s1, interaction = d1,
                                      annotation = anno, etype = 'interaction')
                        DGnt.add_edge(n2, n1, weight = s2, interaction = d2,
                                      annotation = anno, etype = 'interaction')

        ## Normal graph creation + non-expression edges
        if direction == '-':
            DG.add_edge(n1, n2, weight = score, interaction = direction, annotation = anno)
            DG.add_edge(n2, n1, weight = score, interaction = direction, annotation = anno)
            if not expInteraction:
                DGtf.add_edge(n1, n2, weight = score, interaction = direction, annotation = anno, etype = 'interaction')
                DGtf.add_edge(n2, n1, weight = score, interaction = direction, annotation = anno, etype = 'interaction')

                DGnt.add_edge(n1, n2, weight = score, interaction = direction, annotation = anno, etype = 'interaction')
                DGnt.add_edge(n2, n1, weight = score, interaction = direction, annotation = anno, etype = 'interaction')

        elif len(direction) == 2:
            if direction == '->' or direction == '-|':
                if signed and direction == '-|': sc = -score
                else: sc = score
                DG.add_edge(n1, n2, weight = sc, interaction = direction, annotation = anno)
                if not expInteraction:
                    DGtf.add_edge(n1, n2, weight = sc, interaction = direction, annotation = anno, etype = 'interaction')
                    DGnt.add_edge(n1, n2, weight = sc, interaction = direction, annotation = anno, etype = 'interaction')

            elif direction == '<-' or direction == '|-':
                if direction == '<-': direction = '->'
                elif direction == '|-': direction = '-|'
                if signed and direction == '-|': sc = -score
                else: sc = score
                DG.add_edge(n2, n1, weight = sc, interaction = direction, annotation = anno)
                if not expInteraction:
                    DGtf.add_edge(n2, n1, weight = sc, interaction = direction, annotation = anno, etype = 'interaction')
                    DGnt.add_edge(n2, n1, weight = sc, interaction = direction, annotation = anno, etype = 'interaction')

        elif len(direction) == 3:
            d1 = direction[1:]
            d2 = direction[0:2]
            s1, s2 = score, score
            if d2 == '<-': d2 = '->'
            elif d2 == '|-': d2 = '-|'
            if signed and d1 == '-|': s1 = -s1
            if signed and d2 == '-|': s2 = -s2

            DG.add_edge(n1, n2, weight = s1, interaction = d1, annotation = anno)
            DG.add_edge(n2, n1, weight = s2, interaction = d2, annotation = anno)
            if not expInteraction:
                DGtf.add_edge(n1, n2, weight = s1, interaction = d1, annotation = anno, etype = 'interaction')
                DGtf.add_edge(n2, n1, weight = s2, interaction = d2, annotation = anno, etype = 'interaction')
                DGnt.add_edge(n1, n2, weight = s1, interaction = d1, annotation = anno, etype = 'interaction')
                DGnt.add_edge(n2, n1, weight = s2, interaction = d2, annotation = anno, etype = 'interaction')

# Create undirected graph
if not signed:
    G = DG.to_undirected()
    Gtf = DGtf.to_undirected()
    Gnt = DGnt.to_undirected()
else:
    G, Gtf, Gnt = nx.DiGraph(), nx.DiGraph(), nx.DiGraph()
    for edge in DG.edges:
        e1,e2 = edge
        w = DG[e1][e2]['weight']
        inter = DG[e1][e2]['interaction']
        anno = DG[e1][e2]['annotation']
        #et = DG[e1][e2]['etype']
        G.add_edge(e1, e2, weight = w, interaction = inter, annotation = anno)
        if (e2, e1) not in DG.edges:
            G.add_edge(e2, e1, weight = w, interaction = inter, annotation = anno)
    for edge in DGtf.edges:
        e1,e2 = edge
        w = DGtf[e1][e2]['weight']
        inter = DGtf[e1][e2]['interaction']
        anno = DGtf[e1][e2]['annotation']
        et = DGtf[e1][e2]['etype']
        Gtf.add_edge(e1, e2, weight = w, interaction = inter, annotation = anno, etype = et)
        if (e2, e1) not in DGtf.edges:
            Gtf.add_edge(e2, e1, weight = w, interaction = inter, annotation = anno, etype = et)
    for edge in DGnt.edges:
        e1,e2 = edge
        w = DGnt[e1][e2]['weight']
        inter = DGnt[e1][e2]['interaction']
        anno = DGnt[e1][e2]['annotation']
        et = DGnt[e1][e2]['etype']
        Gnt.add_edge(e1, e2, weight = w, interaction = inter, annotation = anno, etype = et)
        if (e2, e1) not in DGnt.edges:
            Gnt.add_edge(e2, e1, weight = w, interaction = inter, annotation = anno, etype = et)

# Write graphs to file
_pickle.dump(DG, open(baseo + '_directed_graph.pkl', 'wb'))
_pickle.dump(G, open(baseo + '_undirected_graph.pkl', 'wb'))

_pickle.dump(DGtf, open(baseo + '_with_transcriptional_edges_directed_graph.pkl', 'wb'))
_pickle.dump(Gtf, open(baseo + '_with_transcriptional_edges_undirected_graph.pkl', 'wb'))

_pickle.dump(DGnt, open(baseo + '_exclude_transcriptional_edges_directed_graph.pkl', 'wb'))
_pickle.dump(Gnt, open(baseo + '_exclude_transcriptional_edges_undirected_graph.pkl', 'wb'))

print('Total nodes: %d' % (len(G.nodes)))
print('Total directed edges: %d' % (len(DG.edges)))
print('Total undirected edges: %d' % (len(G.edges)))

print()
print('Total nodes: %d' % (len(Gtf.nodes)))
print('Total directed edges: %d' % (len(DGtf.edges)))
print('Total undirected edges: %d' % (len(Gtf.edges)))

print()
print('Total nodes: %d' % (len(Gnt.nodes)))
print('Total directed edges: %d' % (len(DGnt.edges)))
print('Total undirected edges: %d' % (len(Gnt.edges)))

