## Repo for creating networkx graph objects from ReactomFI data (source: https://reactome.org/download-data)

- reactomeFI_create_graph.py
    - script to parse ReactomeFI text files and generate several networkx graph objects
    - handles sign and direction information on edges when present

- adjust_graph_*directed.py scripts
    - scripts to add additional edges to graphs
    - set up for CCLE RPPA phosphosites and CCLE metabolomics data metabolite PMIs

data_files/
    - data files used to help graph construction as well as additions
    - Files:
        - GO_transcriptional_regulators.txt - Gene ontology annotated transcriptional regulators
        - TF_names_v_1.01.txt - transcription factor names from Human TFs DB
            - Source file: Google Drive/Shared drives/Hunter Bio Team/Anthony/databases/HumanTFsDB/TF_names_v_1.01.txt
        - CCLE_RPPA_phosphosites.txt - RPPA phosphosites included in CCLE data
        - metabolite_HMDB_conversion_with_protein_interactions.txt - CCLE metabolomics data metabolite HMDB mappings with PMIs

