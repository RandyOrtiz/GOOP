# GOOP - Gene Ontology Overlap Profiler

**GOOP** is a comprehensive pipeline for downloading, processing, and analyzing genomic data. It automates the steps of genome retrieval, gene prediction, BLAST searches, Gene Ontology (GO) term processing, and visualization through heatmaps.

Gene Ontology (GO) analysis serves as a foundational tool in bioinformatics for interpreting functional differences among organisms by categorizing genes into standardized terms associated with primary metabolic processes, molecular functions, and cellular components. Despite its utility, comparative GO analyses between groups of microorganisms remain computationally intensive due to the lack of streamlined, user-friendly tools. Here, we introduce Gene Ontology Overlap Profiler (GOOP), a novel bioinformatics pipeline designed to facilitate the comparison of GO term distributions between two groups of microorganisms. GOOP automates the extraction, normalization, comparison, and visualization of GO terms directly from Bacterial and Viral Bioinformatics Resource Center (BV-BRC) data tables and local genome datasets, providing a cohesive workflow for detecting functional divergences. Key features of GOOP include customizable GO term filtering, integration of comprehensive statistical analyses for identifying significantly enriched or depleted GO categories, and the generation of detailed reports accompanied by an informative heatmap-like graphic. 

## Table of Contents

- [Overview](#overview)
- [Dependencies](#dependencies)
  - [Python Packages](#python-packages)
  - [External Tools](#external-tools)
- [Installation](#installation)
  - [Python Environment Setup](#python-environment-setup)
  - [Installing External Tools](#installing-external-tools)
- [Sample Results](#sample-results)
- [Usage](#usage)
- [Example Command](#examples)
- [License](#license)

## Overview

The pipeline performs the following steps:

1. **Download and Rename Genomes**: Retrieves genome assemblies based on provided accession numbers and renames them for consistency.
2. **Process with Prodigal**: Unzips, linearizes, and predicts genes using Prodigal.
3. **Download UniProt Data**: Retrieves protein sequences and annotations from UniProt based on Taxonomy ID and optional GO terms.
4. **Run BLAST**: Creates a BLAST database and runs BLASTX to find homologous proteins.
5. **Process BLAST Results**: Matches BLAST results with GO terms and generates summary CSV files.
6. **Comparison and Visualization**: Compares GO terms between datasets and generates heatmaps for visualization.

## Dependencies

### Python Packages

- Python 3.6 or higher
- `argparse`
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `goatools`
- `requests`
- `Biopython` (optional, for advanced FASTA handling)

### External Tools

- **Prodigal**: Gene prediction tool.

  - Website: [Prodigal](https://github.com/hyattpd/Prodigal)
  - Version: 2.6.3 or higher

- **NCBI BLAST+**: Sequence alignment tool suite.

  - Website: [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  - Tools Used: `makeblastdb`, `blastx`
  - Version: 2.10.0 or higher

- **Entrez Direct (EDirect)**: Command-line tool for accessing NCBI's databases.

  - Installation Instructions: [EDirect Documentation](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
  - Tools Used: `esearch`, `esummary`, `xtract`

- **wget**: Network downloader for retrieving files from the web.

  - Installation: Available in most package managers (`apt`, `yum`, `brew`, etc.)

## Installation

### Python Environment Setup

It's recommended to use a virtual environment to manage the Python packages.

1. **Create a Virtual Environment**

   ```bash
   python3 -m venv goop_env
   ```

2. **Activate the Virtual Environment**

   - On Linux/macOS:

     ```bash
     source goop_env/bin/activate
     ```

   - On Windows:

     ```bash
     goop_env\Scripts\activate
     ```

3. **Upgrade `pip`**

   ```bash
   pip install --upgrade pip
   ```

4. **Install Python Packages**

   ```bash
   pip install pandas numpy matplotlib seaborn goatools requests biopython
   ```

### Installing External Tools

#### Prodigal

1. **Download Prodigal**

   ```bash
   wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux
   ```

2. **Make the Binary Executable**

   ```bash
   chmod +x prodigal.linux
   ```

3. **Move to a Directory in PATH**

   ```bash
   sudo mv prodigal.linux /usr/local/bin/prodigal
   ```

#### NCBI BLAST+

1. **Download BLAST+**

   - Visit the [NCBI BLAST+ Download Page](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and download the appropriate package for your operating system.

2. **Extract and Install**

   ```bash
   tar -xzf ncbi-blast-*.tar.gz
   sudo cp ncbi-blast-*/bin/* /usr/local/bin/
   ```

#### Entrez Direct (EDirect)

1. **Install EDirect**

   ```bash
   sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
   ```

2. **Update PATH Environment Variable**

   ```bash
   export PATH=${PATH}:${HOME}/edirect
   ```

   - Add the above line to your `~/.bashrc` or `~/.bash_profile` to make it permanent.

#### wget

- **Install via Package Manager**

  - On Ubuntu/Debian:

    ```bash
    sudo apt-get install wget
    ```

  - On CentOS/RHEL:

    ```bash
    sudo yum install wget
    ```

  - On macOS (with Homebrew):

    ```bash
    brew install wget
    ```

## Sample Results

- Final_Results_001_query_counts_report.tsv (This should be used to validate final results png)
  - This file describes all GO Terms that are unique to the Query (and thus completely absent from all Subject)
    - First Column: GO Terms
    - Second Column: Count (Number of files in which the GO Term is present)
    - Third Column: Percentage (Percentage of files in which the GO Term is present)
    - Fourth Column: Similar GOs (Similar GOs utilizing GO node proximity)
    - Fifth Column: Similar Strings (Similar GOs utilizing similar strings)
      
- Final_Results_001_subject_counts_report.tsv (This should be used to validate final results png)
  - This file describes all GO Terms that are unique to the Subject (and thus completely absent from all Query)
    - First Column: GO Terms
    - Second Column: Count (Number of files in which the GO Term is present)
    - Third Column: Percentage (Percentage of files in which the GO Term is present)
    - Fourth Column: Similar GOs (Similar GOs utilizing GO node proximity)
    - Fifth Column: Similar Strings (Similar GOs utilizing similar strings)
      
- GOs_present_in_both_query.csv
  - GO Terms that are present in both  
  - First row: GO Terms
  - Second row: Genome comparison - Query sided
  - Last row:
    
- GOs_present_in_both_subject.csv
  - GO Terms that are present in both
  - First row: GO Terms
  - Second row: Genome comparison - Subject sided
  - Last row:
    
- combined_heatmaps.png
  - Final results comparison between query and subject (utilize Final_Results csv files to validate these results)

![image](https://github.com/user-attachments/assets/4389ab9e-dfb1-459c-b0df-f17a15330c2f)

## Usage

### Arguments

- Subject Genomes (choose -c or -s)
  - `-c`, `--bvbrc-csv`: Path to the input BV-BRC CSV file containing genome accession numbers (this will download genomes from the NCBI database)
  - `-s`, `--subject-genomes-folder`: Path to folder containing subject genomes
- Subject Information (choose -st and -sg; or else choose -sc)
  - `-st`, `--subject-genomes-taxid`: Taxonomy ID for the subject organism (e.g., `194` for *Campylobacter*)
  - `-sg`, `--subject-genomes-go`: GO term to filter the subject organism (without "GO:" prefix)
  - `-sc`, `--subject-folder-complete`: Subject folder already has been run with GOOP, skip subject processing steps (must choose 1b as well)
- Query Genomes
  - `-q`, `--query-genomes-folder`: Path to folder containing query genomes
- Query Information (choose -qt and -qg; or else choose -qc)
  - `-qt`, `--query-genomes-taxid`: Taxonomy ID for the query organism (e.g., `194` for *Campylobacter*)
  - `-qg`, `--query-genomes-go`: GO term to filter the query organism (without "GO:" prefix)
  - `-qc`, `--query-folder-complete`: Query folder already has been run with GOOP, skip subject processing steps
- Other Mandatory Arguments
  - `-f`, `--comparison-go`: GO term to be used for filtering both query and subject (without the "GO:" prefix)
  - `-p`, `--prodigal`: Path to the `prodigal` command
  - `-m`, `--makeblastdb`: Path to the `makeblastdb` command
  - `-b`, `--blastx`: Path to the `blastx` command
  - `-t`, `--threads`: Number of threads to use for BLASTX (default: `24`).
- Annotation Optional Arguments
  - `--evalue`: Annotation e-value (default: 1e-10)
  - `--qcov`: Anotation minimum query coverage (default: 80) %
  - `--min-seq-id`: Anotation minimum query coverage for MMseqs2 (default: 0.8)
  - `--sensitivity`: MMseqs2 sensitivity level (default: 8.0) very sensitive
- Other Optional Arguments (to avoid downloads)
  - `-sf`, `--subject-uniprot-fasta`: Path to a local UniProt FASTA file to use instead of downloading
  - `-sv`, `--subject-uniprot-tsv`: Path to a local UniProt TSV file to use instead of downloading
  - `-qf`, `--query-uniprot-fasta`: Path to a local UniProt FASTA file to use instead of downloading
  - `-qv`, `--query-uniprot-tsv`: Path to a local UniProt TSV file to use instead of downloading
  - `-o`, `--obo-file`: Path to a local OBO file to use instead of downloading
- Other Optional Modes
  - `--metagenome-subject`: Use this if your subject genomes are metagenome (prodigal -p meta)
  - `--metagenome-query`: Use this if your query genomes are metagenome (prodigal -p meta)
  - `-d`, `--diamond`: Path to the `diamond` command
  - `--use-diamond`:Use DIAMOND instead of BLASTX
  - `-mm`, `--mmseqs`: Path to the `mmseqs2` command
  - `--use-mmseqs`: Use MMseqs2 instead of BLASTX

## Examples

## Sample GOOP Runs with:
1. remote subject (will download subject genomes from NCBI)
2. local query
3. remote UniProt databases (will download databases from UniProt)

```bash
python GOOP.py \
    --bvbrc-csv /path_to_subject_BVBRC_genome_tsv_file/BVBRC_test.csv \
    --subject-genomes-taxid 194 \
    --subject-genomes-go 0008150 \
    --query-genomes-folder /path_to_query_genomes_folder_in_fasta_format/query \
    --query-genomes-taxid 194 \
    --query-genomes-go 0008150 \
    --comparison-go 0044238 \
    --prodigal /path_to_prodigal_command/prodigal \
    --makeblastdb /path_to_ncbi_blast_makeblastdb_command/makeblastdb \
    --blastx /path_to_blastx_command/blastx \
    --threads 24
```

1. local subject
2. local query
3. local UniProt databases

```bash
python GOOP.py \
    --subject-genomes-folder /path_to_subject_genomes_folder_in_fasta_format/genomes \
    --subject-genomes-taxid 194 \
    --subject-genomes-go 0008150 \
    --query-genomes-folder /path_to_query_genomes_folder_in_fasta_format/query \
    --query-genomes-taxid 194 \
    --query-genomes-go 0008150 \
    --comparison-go 0044238 \
    --prodigal /path_to_prodigal_command/prodigal \
    --makeblastdb /path_to_ncbi_blast_makeblastdb_command/makeblastdb \
    --blastx /path_to_blastx_command/blastx \
    --threads 24 \
    --subject-uniprot-fasta /path_to_uniprot_database_in_fasta_format/uniprot_194_proteins.fasta \
    --subject-uniprot-tsv /path_to_uniprot_tsv/uniprot_194_proteins_filtered_by_go_0005975.tsv \
    --query-uniprot-fasta /path_to_uniprot_database_in_fasta_format/uniprot_194_proteins.fasta \
    --query-uniprot-tsv /path_to_uniprot_tsv/uniprot_194_proteins_filtered_by_go_0005975.tsv \
    --obo-file /path_to_obo_file/go-basic.obo
```

1. Subject folder has been run with GOOP
2. Query folder has been run with GOOP

```bash
python GOOP.py \
    --subject-genomes-folder /path_to_subject_genomes_folder_in_fasta_format/genomes \
    --query-genomes-folder /path_to_query_genomes_folder_in_fasta_format/query \
    --comparison-go 0044238 \
    --obo-file /path_to_obo_file/go-basic.obo \
    --subject-folder-complete \
    --query-folder-complete
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Note:** Ensure all external tools are accessible from your `PATH` environment variable or provide the full paths to their executables when running the script.

If you encounter any issues or have questions, please open an issue on the repository or contact the maintainer.
