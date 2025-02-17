# Code repository for various network modeling methods/functions

- Testing signed and directed propagation

# 2022-04-28
- Major bug fix in network propagation function:
    - indexing of nodes from subgraphs was not being performed correctly, so assignments random to some extent
    - Old line:
        - sub_inds = [i for i,x in enumerate(nodes) if x in subgraph_nodes] # 
    - New fixed lines:
        - sub_indsD = {x: i for i,x in enumerate(nodes) if x in subgraph_nodes}
          sub_inds = [sub_indsD[x] for x in subgraph_nodes]

# 2024-06-06
- Include interactome_tools sub-folder with code for handling interactomes

