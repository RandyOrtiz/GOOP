# GOOP
GOOP - Gene Ontology Overlap Profiler

Gene Ontology (GO) analysis serves as a foundational tool in bioinformatics for interpreting functional differences among organisms by categorizing genes into standardized terms associated with primary metabolic processes, molecular functions, and cellular components. Despite its utility, comparative GO analyses between groups of microorganisms remain computationally intensive due to the lack of streamlined, user-friendly tools. Here, we introduce Gene Ontology Overlap Profiler (GOOP), a novel bioinformatics pipeline designed to facilitate the comparison of GO term distributions between two groups of microorganisms. GOOP automates the extraction, normalization, comparison, and visualization of GO terms directly from Bacterial and Viral Bioinformatics Resource Center (BV-BRC) data tables and local genome datasets, providing a cohesive workflow for detecting functional divergences. Key features of GOOP include customizable GO term filtering, integration of comprehensive statistical analyses for identifying significantly enriched or depleted GO categories, and the generation of detailed reports accompanied by an informative heatmap-like graphic. 

Sample GOOP Run with remote subject, local query, and remote UniProt databases.
Also can be used to test GOOP installation.

python GOOP.py \
    -i /path_to_subject_BVBRC_genome_tsv_file/BVBRC_test.csv \
    -p /path_to_prodigal_command/prodigal.linux \
    -d 194 \ # taxid of subject
    -g 0005975 \ # outer GO id of interest in subject
    --filtering-go 0044238 \ # inner GO id of interest in analysis
    -m /path_to_ncbi_blast_makeblastdb_command/makeblastdb \
    -b /path_to_blastx_command/blastx \
    -t 24 \ # number of threads for GOOP
    --process-folder /path_to_query_genomes_folder_in_fasta_format/query \
    --process-folder-taxid 194 \ # taxid of query
    --process-folder-go 0005975 \ # outer GO id of interest in query

Sample GOOP run with remote subject, local query, and local UniProt databases

python GOOP.py \
    -i /path_to_subject_BVBRC_genome_tsv_file/BVBRC_test.csv \
    -p /path_to_prodigal_command/prodigal.linux \
    -d 194 \ # taxid of subject
    -g 0005975 \ # outer GO id of interest in subject
    --filtering-go 0044238 \ # inner GO id of interest in analysis
    -m /path_to_ncbi_blast_makeblastdb_command/makeblastdb \
    -b /path_to_blastx_command/blastx \
    -t 24 \ # number of threads for GOOP
    --process-folder /path_to_query_genomes_folder_in_fasta_format/query \
    --process-folder-taxid 194 \ # taxid of query
    --process-folder-go 0005975 \ # outer GO id of interest in query
    --uniprot-fasta /home/velox/Documents/GOOP_PeerJ_Paper/Downloads/Final_Test/uniprot_194_proteins.fasta \
    --uniprot-tsv /home/velox/Documents/GOOP_PeerJ_Paper/Downloads/Final_Test/uniprot_194_proteins_filtered_by_go_0005975.tsv \
    --obo-file /home/velox/Documents/GOOP_PeerJ_Paper/Downloads/Final_Test/go-basic.obo \
    --process-folder-uniprot-fasta /home/velox/Documents/GOOP_PeerJ_Paper/Downloads/Final_Test/uniprot_194_proteins.fasta \
    --process-folder-uniprot-tsv /home/velox/Documents/GOOP_PeerJ_Paper/Downloads/Final_Test/uniprot_194_proteins_filtered_by_go_0005975.tsv
