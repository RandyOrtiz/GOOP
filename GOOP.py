#!/usr/bin/env python3

import argparse
import os
import subprocess
import csv
import pandas as pd
import shutil
import sys
import gzip
from pathlib import Path
import re
from collections import defaultdict
from goatools.obo_parser import GODag
import seaborn as sns
import matplotlib.pyplot as plt
import requests


# First Part: Download and Rename Genomes
def download_and_rename(csv_file, prodigal_path):
    # Check if EDirect is installed
    if not shutil.which('esearch'):
        print("Error: NCBI Entrez Direct (EDirect) is not installed.")
        print("Please install it before running this script.")
        sys.exit(1)

    # Create output directory
    output_dir = 'genomes'
    os.makedirs(output_dir, exist_ok=True)

    # Open the CSV file and extract the relevant columns
    with open(csv_file, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        if 'Assembly Accession' not in reader.fieldnames or 'Genome Name' not in reader.fieldnames:
            print("Error: 'Assembly Accession' or 'Genome Name' column not found in the CSV file.")
            sys.exit(1)

        for row in reader:
            accession = row['Assembly Accession'].strip()
            genome_name = row['Genome Name'].strip()
            if not accession or not genome_name:
                continue

            print(f"Processing {accession} for genome {genome_name}")

            # Get the FTP path using EDirect
            ftp_path = get_ftp_path(accession)
            if not ftp_path:
                print(f"No FTP path found for {accession}")
                continue

            # Construct the download URL
            file_name = os.path.basename(ftp_path)
            file_url = f"{ftp_path}/{file_name}_genomic.fna.gz"

            # Download the file
            output_file = os.path.join(output_dir, f"{accession}_genomic.fna.gz")
            if os.path.exists(output_file):
                print(f"File {output_file} already exists. Skipping download.")
                continue

            print(f"Downloading {file_url}")
            download_success = download_file(file_url, output_file)
            if not download_success:
                print(f"Failed to download {file_url}")
                continue

            # Rename the file using the format 'Genome Name__Assembly Accession_genomic.fna.gz'
            new_file_name = f"{genome_name}__{accession}_genomic.fna.gz"
            new_file_path = os.path.join(output_dir, new_file_name)
            os.rename(output_file, new_file_path)
            print(f"Renamed {output_file} to {new_file_path}")

    print(f"Download and renaming completed. Genomes are saved in the '{output_dir}' directory.")

def get_ftp_path(accession):
    try:
        result = subprocess.run(
            ['esearch', '-db', 'assembly', '-query', accession],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        esummary = subprocess.run(
            ['esummary'],
            input=result.stdout,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        ftp_path = subprocess.run(
            ['xtract', '-pattern', 'DocumentSummary', '-element', 'FtpPath_RefSeq'],
            input=esummary.stdout,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout.strip()

        # If RefSeq path is empty, try GenBank FTP path
        if not ftp_path:
            ftp_path = subprocess.run(
                ['xtract', '-pattern', 'DocumentSummary', '-element', 'FtpPath_GenBank'],
                input=esummary.stdout,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout.strip()

        return ftp_path
    except Exception as e:
        print(f"Error retrieving FTP path for {accession}: {e}")
        return None

def download_file(url, output_path):
    try:
        subprocess.run(['wget', '-q', url, '-O', output_path], check=True)
        return True
    except subprocess.CalledProcessError:
        return False

# Second Part: Process with Prodigal
def gunzip_files_in_folder(folder_path):
    for filename in os.listdir(folder_path):
        if filename.endswith('.gz'):
            gz_file_path = os.path.join(folder_path, filename)
            output_file_path = os.path.join(folder_path, filename[:-3])
            try:
                with gzip.open(gz_file_path, 'rb') as f_in:
                    with open(output_file_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(gz_file_path)
                print(f'Unzipped and removed {gz_file_path}')
            except Exception as e:
                print(f'Failed to unzip {gz_file_path}: {e}')

def linearize_and_process_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    linearized_content = ""
    for line in lines:
        if line.strip():
            if line.startswith('>'):
                if linearized_content:
                    linearized_content += '\n'
                linearized_content += line.strip() + '\n'
            else:
                linearized_content += line.strip()
    new_file_name = file_path.stem.replace('_genomic', '_linear') + '.fasta'
    new_file_path = file_path.with_name(new_file_name)
    with open(new_file_path, 'w') as file:
        file.write(linearized_content)
    print(f"Linearized file created: {new_file_path}")
    with open(new_file_path, 'r') as file:
        lines = file.readlines()
    processed_content = ""
    for line in lines:
        if line.startswith('>'):
            header = re.sub(r'[^A-Za-z0-9>]+', '_', line.strip())
            processed_content += header + '\n'
        else:
            processed_content += line.strip() + '\n'
    with open(new_file_path, 'w') as file:
        file.write(processed_content)
    print(f"Processed file created: {new_file_path}")
    return new_file_path

def run_prodigal(input_file, prodigal):
    prodigal_output_name = input_file.stem.replace('_linear', '_prodigal') + '.fasta'
    prodigal_output_path = input_file.with_name(prodigal_output_name)
    subprocess.run([prodigal, '-i', str(input_file), '-d', str(prodigal_output_path)])
    print(f"Prodigal output created: {prodigal_output_path}")

def process_fasta_files(input_folder, prodigal):
    gunzip_files_in_folder(input_folder)
    for file_path in Path(input_folder).iterdir():
        if file_path.is_file() and file_path.suffix in ['.fna', '.fasta', '.fa']:
            print(f"Processing file: {file_path}")
            linearized_file = linearize_and_process_fasta(file_path)
            run_prodigal(linearized_file, prodigal)
    print("Processing complete.")

# Additional Functions from the Second Script
def download_uniprot_data(taxid, go_term=None, output_dir='.'):
    """
    Downloads UniProt FASTA and TSV data based on the provided Taxonomy ID and optional GO term.
    If the files already exist and are non-empty, it skips downloading.
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Define output file paths
    fasta_output_file = os.path.join(output_dir, f"uniprot_{taxid}_proteins.fasta")
    if not go_term:
        tsv_output_file = os.path.join(output_dir, f"uniprot_{taxid}_proteins.tsv")
    else:
        tsv_output_file = os.path.join(output_dir, f"uniprot_{taxid}_proteins_filtered_by_go_{go_term}.tsv")

    # Check if FASTA file exists and is non-empty
    if os.path.exists(fasta_output_file) and os.path.getsize(fasta_output_file) > 0:
        print(f"FASTA file {fasta_output_file} already exists and is non-empty. Skipping download.")
    else:
        # Download the FASTA file
        fasta_url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A{taxid}%29%29"
        fasta_wget_command = f"wget \"{fasta_url}\" -O \"{fasta_output_file}\""
        print(f"Downloading FASTA file with Taxonomy ID {taxid}...")
        try:
            subprocess.run(fasta_wget_command, shell=True, check=True)
            print(f"FASTA file downloaded as {fasta_output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading FASTA file: {e}")
            sys.exit(1)

    # Check if TSV file exists and is non-empty
    if os.path.exists(tsv_output_file) and os.path.getsize(tsv_output_file) > 0:
        print(f"TSV file {tsv_output_file} already exists and is non-empty. Skipping download.")
    else:
        # Download the TSV file
        if not go_term:
            tsv_url = f"https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cgo_id%2Corganism_id%2Cgo_p%2Cgo_f%2Cxref_unipathway&format=tsv&query=%28%28taxonomy_id%3A{taxid}%29%29"
        else:
            tsv_url = f"https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cgo_id%2Corganism_id%2Cgo_p%2Cgo_f%2Cxref_unipathway&format=tsv&query=%28%28taxonomy_id%3A{taxid}%29+AND+%28go%3A{go_term}%29%29"
        tsv_wget_command = f"wget \"{tsv_url}\" -O \"{tsv_output_file}\""
        print(f"Downloading TSV file{' with GO term ' + go_term if go_term else ''}...")
        try:
            subprocess.run(tsv_wget_command, shell=True, check=True)
            print(f"TSV file downloaded as {tsv_output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading TSV file: {e}")
            sys.exit(1)

    return fasta_output_file, tsv_output_file


def run_blast(fna_folder, protein_db_fasta, makeblastdb_path, blastx_path, num_threads, output_dir='.'):
    """
    Creates a BLAST database and runs BLASTX on all relevant .fasta files in the specified folder.
    """
    # Ensure the folder exists
    if not os.path.isdir(fna_folder):
        print(f"Error: Folder '{fna_folder}' does not exist.")
        sys.exit(1)

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Check if the protein database FASTA file exists and is non-empty
    if not os.path.exists(protein_db_fasta) or os.path.getsize(protein_db_fasta) == 0:
        print(f"Error: The protein database FASTA file '{protein_db_fasta}' does not exist or is empty.")
        sys.exit(1)

    # Create the protein database
    print(f"Creating BLAST protein database from {protein_db_fasta}...")
    try:
        db_output_path = os.path.join(output_dir, 'protein_db')
        makeblastdb_command = [
            makeblastdb_path, "-in", protein_db_fasta, "-dbtype", "prot", "-out", db_output_path
        ]
        subprocess.run(makeblastdb_command, check=True)
        print(f"Protein database created successfully at {db_output_path}.")
    except subprocess.CalledProcessError as e:
        print(f"Error creating protein database: {e}")
        sys.exit(1)

    # **Define fna_files before using it**
    # List all .fasta files ending with _prodigal.fasta in the folder
    fna_files = [f for f in os.listdir(fna_folder) if f.endswith('_prodigal.fasta')]
    if not fna_files:
        print(f"No .fasta files ending with '_prodigal.fasta' found in folder '{fna_folder}'.")
        sys.exit(1)

    for fna_file in fna_files:
        input_file = os.path.join(fna_folder, fna_file)
        output_file = os.path.join(fna_folder, fna_file.replace('_prodigal.fasta', '_blast_results.txt'))
        # Construct the BLASTX command
        blastx_command = [
            blastx_path, "-query", input_file, "-db", db_output_path, "-out", output_file,
            "-evalue", "1e-10", "-max_target_seqs", "1", "-qcov_hsp_perc", "80", "-outfmt", "6",
            "-num_threads", str(num_threads)
        ]
        print(f"Running BLASTX for {fna_file} with {num_threads} threads...")
        try:
            subprocess.run(blastx_command, check=True)
            print(f"Results saved to {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error running BLASTX for {fna_file}: {e}")

def extract_first_two_columns(input_file, output_file):
    """
    Extracts the first two columns from the TSV input file and writes them to the output file.
    """
    print(f"Extracting first two columns from {input_file} to {output_file}...")
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        for row in reader:
            writer.writerow(row[:2])
    print(f"Extraction complete. Saved to {output_file}")

def process_files(input_folder, reference_file, output_folder):
    """
    Processes BLAST result files by matching GO terms with the reference dictionary.
    """
    print(f"Processing BLAST results in {input_folder} using reference file {reference_file}...")
    with open(reference_file, 'r') as file:
        reference_dict = {line.split('\t')[0]: line.split('\t')[1] for line in file if len(line.split('\t')) > 1}

    os.makedirs(output_folder, exist_ok=True)
    for filename in os.listdir(input_folder):
        if filename.endswith('_blast_results.txt'):
            with open(os.path.join(input_folder, filename), 'r') as infile:
                lines = infile.readlines()
            output_filepath = os.path.join(output_folder, filename + '_matched_GOs.txt')
            with open(output_filepath, 'w') as outfile:
                for line in lines:
                    parts = line.strip().split('\t')
                    if len(parts) > 1:
                        second_col = parts[1]
                        match = second_col.split('|')[1] if '|' in second_col else None
                        if match and match in reference_dict:
                            outfile.write(reference_dict[match] + '\n')
    print(f"GO matching complete. Matched files are saved in {output_folder}")

def replace_semicolon_and_remove_spaces(folder_path):
    """
    Replaces semicolons with newlines and removes spaces in all files within the specified folder.
    """
    print(f"Replacing semicolons and removing spaces in files within {folder_path}...")
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            with open(file_path, 'r') as file:
                file_contents = file.readlines()
            updated_contents = [line.replace(';', '\n').replace(' ', '').lstrip() for line in file_contents]
            with open(file_path, 'w') as file:
                file.writelines(updated_contents)
    print(f"Semicolon replacement and space removal complete for files in {folder_path}")

def process_folder(folder_path, output_csv):
    """
    Processes the folder containing matched GO files and generates a summary CSV.
    """
    print(f"Processing matched GO files in {folder_path} to generate {output_csv}...")
    data = {}
    for filename in sorted(os.listdir(folder_path)):
        file_path = os.path.join(folder_path, filename)
        if not filename.endswith('_matched_GOs.txt'):
            continue
        row_header = filename.split('_matched_GOs.txt')[0]
        with open(file_path, 'r') as file:
            strings = file.read().splitlines()
            string_counts = {string: strings.count(string) for string in strings}
            data[row_header] = string_counts
    all_unique_strings = sorted(set().union(*(counts.keys() for counts in data.values())))
    df = pd.DataFrame(columns=all_unique_strings)
    for row_header, string_counts in data.items():
        df.loc[row_header] = [string_counts.get(string, 0) for string in all_unique_strings]
    df.index = df.index.str.replace('_blast_results.txt', '', regex=False)
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'organism'}, inplace=True)
    df.to_csv(output_csv, index=False)
    print(f"Summary CSV generated at {output_csv}")
    shutil.rmtree(folder_path)
    print(f"Temporary folder {folder_path} has been removed.")

def process_file(input_csv, output_csv, output_csv_uniq):
    """
    Sorts the input CSV by the sum of rows, removes duplicates, and saves both sorted and unique CSVs.
    """
    print(f"Processing file {input_csv} to generate sorted and unique CSVs...")
    if not os.path.exists(input_csv):
        print(f"Error: The file '{input_csv}' does not exist.")
        sys.exit(1)
    data = pd.read_csv(input_csv, index_col=0)
    sorted_data = data.loc[data.sum(axis=1).sort_values(ascending=False).index]
    sorted_data.to_csv(output_csv, index=True)
    print(f"Sorted CSV saved as {output_csv}")
    unique_data = sorted_data.drop_duplicates()
    unique_data.to_csv(output_csv_uniq, index=True)
    print(f"Unique CSV saved as {output_csv_uniq}")

def process_blast_results(tsv_file, fna_folder, output_dir):
    """
    Processes the BLAST results using the provided TSV file.
    """
    # Use the full path of the TSV file provided
    input_file = tsv_file
    output_file = os.path.join(output_dir, "T1_GO_KEY_T1.txt")
    output_folder = os.path.join(output_dir, "T1_GO_KEY_T1")

    # Step 1: Extract first two columns from TSV
    extract_first_two_columns(input_file, output_file)

    # Step 2: Process BLAST results and match GOs
    process_files(fna_folder, output_file, output_folder)

    # Step 3: Replace semicolons and remove spaces in matched GO files
    replace_semicolon_and_remove_spaces(output_folder)

    # Step 4: Process the folder to create a summary CSV
    final_csv = os.path.join(output_dir, "T2_QUERY_GO_TABLE_T2.csv")
    process_folder(output_folder, final_csv)

    # Step 5: Sort and deduplicate the final CSV
    final_csv_uniq = os.path.join(output_dir, "T2_QUERY_GO_TABLE_T2_uniq.csv")
    process_file(final_csv, final_csv, final_csv_uniq)

    # Return the final CSV filenames
    return final_csv, final_csv_uniq

def comparison(unique_file_results_uniq, unique_file_results_uniq_process, filtering_go, obo_file=None):
    """
    A function to compare the results using the provided unique files and filtering GO term.
    """
    print(f"Comparing results using:")
    print(f"Unique File Results 1: {unique_file_results_uniq}")
    print(f"Unique File Results 2: {unique_file_results_uniq_process}")
    print(f"Filtering GO Term: {filtering_go}")

    # Set the directory to the script's location
    script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(script_dir)

    # Always download the OBO file unless provided
    if obo_file and os.path.exists(obo_file):
        print(f"Using supplied OBO file: {obo_file}")
    else:
        if not obo_file:
            obo_file = 'go-basic.obo'
        if not os.path.exists(obo_file) or os.path.getsize(obo_file) == 0:
            print("Downloading go-basic.obo file...")
            url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
            response = requests.get(url)
            with open(obo_file, 'wb') as f:
                f.write(response.content)
            print("go-basic.obo downloaded.")
        else:
            print(f"Using existing OBO file: {obo_file}")

    # Step 0: Import goatools and load the Gene Ontology
    # Load the Gene Ontology
    print("Loading Gene Ontology...")
    godag = GODag(obo_file)
    print("Gene Ontology loaded.")

    # Function to get all descendants of a GO term
    def get_all_descendants(go_id, godag):
        descendants = set()
        def recurse(term):
            if term in descendants:
                return
            descendants.add(term)
            for child in godag[term].children:
                recurse(child.id)
        recurse(go_id)
        return descendants

    # Get all descendants of the parent GO term
    print(f"Getting all descendants of {filtering_go}...")
    descendants = get_all_descendants(filtering_go, godag)
    descendants = set(descendants)  # Convert to set for faster lookup
    print(f"Total descendants of {filtering_go}: {len(descendants)}")

    # Function to standardize GO IDs
    def standardize_go_id(go_id):
        # Replace 'GO_' with 'GO:'
        go_id = go_id.replace('GO_', 'GO:')
        # Remove any spaces
        go_id = go_id.replace(' ', '')
        # Extract 'GO:' followed by digits
        match = re.search(r'GO[:_]?(\d+)', go_id)
        if match:
            return f"GO:{match.group(1)}"
        else:
            return go_id.strip()

    # ------------------------ Processing Function ------------------------

    # Function to process GO terms
    def process_go_terms(df1_filename, df2_filename, godag, filtering_go, output_suffix):
        # Step 1: Read the CSV files
        df1 = pd.read_csv(df1_filename)
        df2 = pd.read_csv(df2_filename)

        # Standardize column names in df1 and df2
        df1.columns = [standardize_go_id(col) for col in df1.columns]
        df2.columns = [standardize_go_id(col) for col in df2.columns]

        # Get the first column name from df1 and df2 (assuming it's the species/organism name)
        first_column_name_df1 = df1.columns[0]
        first_column_name_df2 = df2.columns[0]

        # Step 2: Get GO IDs from the columns of both dataframes and filter to descendants
        go_terms_df1 = set(df1.columns[1:])  # Exclude the first column (species name)
        go_terms_df2 = set(df2.columns[1:])  # Exclude the first column (species name)

        # Filter GO terms to only those that are descendants of filtering_go
        go_terms_df1 = go_terms_df1.intersection(descendants)
        go_terms_df2 = go_terms_df2.intersection(descendants)

        # Step 3: Identify shared and unique GO terms among the filtered terms
        shared_go_terms = go_terms_df1.intersection(go_terms_df2)
        unique_go_terms_df1 = go_terms_df1 - go_terms_df2
        unique_go_terms_df2 = go_terms_df2 - go_terms_df1

        # Adjust unique_go_terms_df1 and unique_go_terms_df2 based on child terms in shared_go_terms and the other organism's GO terms
        def filter_unique_go_terms(unique_go_terms, other_go_terms, shared_go_terms, godag):
            adjusted_unique_go_terms = set()
            for go_id in unique_go_terms:
                descendants_of_go_id = get_all_descendants(go_id, godag)
                # If any of the descendants are in shared_go_terms or in other_go_terms, skip this go_id
                if not descendants_of_go_id.intersection(shared_go_terms.union(other_go_terms)):
                    adjusted_unique_go_terms.add(go_id)
            return adjusted_unique_go_terms

        unique_go_terms_df1 = filter_unique_go_terms(unique_go_terms_df1, go_terms_df2, shared_go_terms, godag)
        unique_go_terms_df2 = filter_unique_go_terms(unique_go_terms_df2, go_terms_df1, shared_go_terms, godag)

        print(f"Total GO terms in df1 (descendants of {filtering_go}): {len(go_terms_df1)}")
        print(f"Total GO terms in df2 (descendants of {filtering_go}): {len(go_terms_df2)}")
        print(f"Number of shared GO terms: {len(shared_go_terms)}")
        print(f"Number of unique GO terms in df1 after adjustment: {len(unique_go_terms_df1)}")
        print(f"Number of unique GO terms in df2 after adjustment: {len(unique_go_terms_df2)}")

        # Step 4: Combine all filtered GO terms for processing
        all_go_terms = go_terms_df1.union(go_terms_df2)

        # Step 5: Create a mapping of GO IDs to descriptions using goatools
        go_descriptions = {}
        for go_id in all_go_terms:
            if go_id in godag:
                go_descriptions[go_id] = godag[go_id].name.lower()
            else:
                go_descriptions[go_id] = ''

        # Function to get full name for a GO ID
        def get_full_name(go_id):
            if go_id in godag:
                return f"{go_id} {godag[go_id].name}"
            else:
                return go_id

        # List of common words to exclude (used later)
        common_words = {'process', 'metabolic', 'catabolic', 'biosynthetic', 'acid', 'via', 'cycle'}

        # Lists to keep track of the output data
        output_a = []  # GO terms unique to df1
        output_b = []  # GO terms unique to df2
        output_c = []  # GO terms shared between df1 and df2

        # Dictionary to store per-organism maximum values for each GO term
        organism_max_values = defaultdict(lambda: defaultdict(float))

        # Step 6: Process each GO ID separately
        for go_id in all_go_terms:
            in_df1 = go_id in df1.columns
            in_df2 = go_id in df2.columns
            output_data = []
            # Compare each row of the first CSV with all rows of the second CSV
            for i, row1 in df1.iterrows():
                species1 = row1[first_column_name_df1]  # Use the first column name from df1
                species1_val = row1.get(go_id, 0) if in_df1 else 0
                for j, row2 in df2.iterrows():
                    species2 = row2[first_column_name_df2]  # Use the first column name from df2
                    species2_val = row2.get(go_id, 0) if in_df2 else 0
                    comparison_val = species1_val - species2_val
                    output_data.append([f"{species1} v. {species2}", comparison_val])
                    # Update the organism's max value for this GO term
                    current_max = organism_max_values[species2].get(go_id, None)
                    if current_max is None or comparison_val > current_max:
                        organism_max_values[species2][go_id] = comparison_val
                # Add a newline of space between comparisons
                output_data.append(["", ""])
            # Append the output data to the appropriate list
            full_name = get_full_name(go_id)  # Get the full GO ID and name
            if output_data:
                if in_df1 and in_df2:
                    output_c.append((full_name, output_data))
                elif in_df1:
                    if go_id in unique_go_terms_df1:
                        output_a.append((full_name, output_data))
                elif in_df2:
                    if go_id in unique_go_terms_df2:
                        output_b.append((full_name, output_data))

        # Function to save output data
        def save_output(output_list, filename, organism_max_values=None, create_counts_report=False, other_go_terms=None):
            merged_df = None
            for full_name, data in output_list:
                df = pd.DataFrame(data, columns=['comparison', full_name])
                if merged_df is None:
                    merged_df = df
                else:
                    # To avoid alignment issues, reset index before concatenation
                    merged_df = pd.concat([merged_df.reset_index(drop=True), df.iloc[:, 1].reset_index(drop=True)], axis=1)
            if merged_df is not None:
                # If organism_max_values is provided, append the max values block
                if organism_max_values:
                    # Add an empty line to separate the blocks
                    empty_row = pd.DataFrame([[""] + [""] * (merged_df.shape[1] - 1)], columns=merged_df.columns)
                    merged_df = pd.concat([merged_df, empty_row], ignore_index=True)
                    # Prepare the max values DataFrame
                    max_values_data = []
                    organisms = sorted(organism_max_values.keys())
                    for organism in organisms:
                        go_values = organism_max_values[organism]
                        row = [organism]
                        for full_name in merged_df.columns[1:]:
                            # Extract GO ID from the full_name using regex
                            match = re.search(r'GO:\d+', full_name)
                            if match:
                                go_id = match.group()
                                value = go_values.get(go_id, '')
                            else:
                                value = ''
                            row.append(value)
                        max_values_data.append(row)
                    max_values_df = pd.DataFrame(max_values_data, columns=merged_df.columns)
                    # Append the max values DataFrame to merged_df
                    merged_df = pd.concat([merged_df, max_values_df], ignore_index=True)
                else:
                    max_values_df = pd.DataFrame()  # Empty DataFrame if organism_max_values is None
                # Now, create counts of negative values per GO term in the organism max values
                if create_counts_report:
                    # For each GO term, count how many organisms have negative values
                    counts = []
                    go_terms = merged_df.columns[1:]
                    for go_term in go_terms:
                        if not max_values_df.empty and go_term in max_values_df.columns:
                            # Get the column from max_values_df
                            values = max_values_df[go_term].dropna()
                            # Convert values to numeric, coerce errors to NaN
                            values = pd.to_numeric(values, errors='coerce')
                            # Count number of negative values
                            count_negative = (values < 0).sum()
                            counts.append(count_negative)
                        else:
                            counts.append(0)  # If max_values_df is empty or column not found, append 0
                    # Create counts DataFrame
                    # First line: column headers (go_terms)
                    headers_row = [""] + list(go_terms)
                    # Second line: counts
                    counts_row = ["Count"] + counts
                    counts_df = pd.DataFrame([headers_row, counts_row], columns=merged_df.columns)
                    # Add an empty line to separate the blocks
                    empty_row = pd.DataFrame([[""] + [""] * (merged_df.shape[1] - 1)], columns=merged_df.columns)
                    # Append empty row and counts_df to merged_df
                    merged_df = pd.concat([merged_df, empty_row, counts_df], ignore_index=True)
                    # Save the merged DataFrame to CSV
                    merged_df.to_csv(filename, index=False)
                    print(f"Saved {filename}")
                    # Now, create the final report (transpose the last two rows)
                    # Transpose counts_df
                    transposed_df = counts_df.transpose()
                    # The first row is empty (since the first column is empty), so we drop it
                    transposed_df = transposed_df.iloc[1:]
                    # Rename columns
                    transposed_df.columns = ['GO Term', 'Count']
                    # Reset index
                    transposed_df = transposed_df.reset_index(drop=True)
                    # Convert 'Count' column to numeric
                    transposed_df['Count'] = pd.to_numeric(transposed_df['Count'], errors='coerce')
                    # Calculate the total number of organisms
                    total_organisms = len(max_values_df)
                    if total_organisms == 0:
                        total_organisms = 1  # To avoid division by zero
                    # Add 'Percentage' column
                    transposed_df['Percentage'] = (transposed_df['Count'] / total_organisms) * 100
                    # Round to 2 decimal places
                    transposed_df['Percentage'] = transposed_df['Percentage'].round(2)
                    # Sort by 'Count', high to low
                    transposed_df = transposed_df.sort_values(by='Count', ascending=False).reset_index(drop=True)
                    # Add 'Similar GOs' and 'Similar Strings' columns
                    similar_gos_list = []
                    similar_strings_list = []
                    for index, row in transposed_df.iterrows():
                        current_go_term = row['GO Term']
                        # Extract GO ID
                        match = re.search(r'GO:\d+', current_go_term)
                        if match:
                            current_go_id = match.group()
                            current_description = go_descriptions.get(current_go_id, '')
                            # Exclude common words
                            current_words = set(current_description.split()) - common_words
                            similarities = []
                            similar_strings = []
                            # Compare with other GO terms
                            for other_go_id in other_go_terms:
                                other_description = go_descriptions.get(other_go_id, '')
                                # Exclude common words
                                other_words = set(other_description.split()) - common_words
                                # Compute Jaccard similarity
                                intersection = current_words & other_words
                                union = current_words | other_words
                                if intersection:
                                    # Similar Strings: any GO term that shares at least one word
                                    similar_strings.append(get_full_name(other_go_id))
                                if union:
                                    similarity = len(intersection) / len(union)
                                    similarities.append((similarity, other_go_id))
                            # Get top 3 similar GO IDs for 'Similar GOs' column
                            top_similar = sorted(similarities, key=lambda x: x[0], reverse=True)[:3]
                            similar_go_terms = [get_full_name(go_id) for _, go_id in top_similar]
                            similar_gos_list.append("; ".join(similar_go_terms))
                            # Add all similar GO terms for 'Similar Strings' column
                            similar_strings_list.append("; ".join(set(similar_strings)))
                        else:
                            similar_gos_list.append("")
                            similar_strings_list.append("")
                    transposed_df['Similar GOs'] = similar_gos_list
                    transposed_df['Similar Strings'] = similar_strings_list
                    # Save to TSV file
                    report_filename = filename.replace('.csv', '_counts_report.tsv')
                    transposed_df.to_csv(report_filename, index=False, sep='\t')
                    print(f"Saved {report_filename}")
                else:
                    # Save the merged DataFrame to CSV without counts report
                    merged_df.to_csv(filename, index=False)
                    print(f"Saved {filename}")

        # Prepare the list of GO IDs from outputs
        def get_go_ids_from_output(output_list):
            go_ids = []
            for full_name, _ in output_list:
                # Extract GO ID from full_name
                match = re.search(r'GO:\d+', full_name)
                if match:
                    go_ids.append(standardize_go_id(match.group()))
            return go_ids

        # Get GO IDs from output_a and output_c for 'other_go_terms'
        other_go_ids = get_go_ids_from_output(output_a) + get_go_ids_from_output(output_c)

        # Save the outputs
        # Define output filenames dynamically based on the input filenames
        base_name1 = os.path.splitext(os.path.basename(df1_filename))[0]
        base_name2 = os.path.splitext(os.path.basename(df2_filename))[0]

        output_a_filename = f'GOs_present_in_{base_name1}_but_missing_in_{base_name2}_{output_suffix}.csv'
        output_b_filename = f'GOs_present_in_{base_name2}_but_missing_in_{base_name1}_{output_suffix}.csv'
        output_c_filename = f'GOs_present_in_both_{base_name1}_{base_name2}_{output_suffix}.csv'

        save_output(output_a, output_a_filename, organism_max_values, create_counts_report=False)
        save_output(output_b, output_b_filename, organism_max_values, create_counts_report=True, other_go_terms=other_go_ids)
        save_output(output_c, output_c_filename, organism_max_values, create_counts_report=False)

        print(f"Processing and merging complete for {output_suffix}.")

        # Return the filenames of the counts reports for use in the heatmap generation
        counts_report_filename = output_b_filename.replace('.csv', '_counts_report.tsv')
        return counts_report_filename

    # ------------------------ Main Processing ------------------------

    print("Starting GO term processing...")

    # First run: df1 as df1 (Query), df2 as df2 (Subject)
    counts_report1 = process_go_terms(unique_file_results_uniq, unique_file_results_uniq_process, godag, filtering_go, 'original')

    # Second run: df1 as df2 (Subject), df2 as df1 (Query)
    counts_report2 = process_go_terms(unique_file_results_uniq_process, unique_file_results_uniq, godag, filtering_go, 'swapped')

    print("All processing complete.")

    # ------------------------ Heatmap Generation ------------------------

    print("Starting heatmap generation...")

    # Read the TSV files generated from the processing function
    df_counts1 = pd.read_csv(counts_report1, sep='\t')
    df_counts2 = pd.read_csv(counts_report2, sep='\t')

    # Prepare and sort data for the first DataFrame
    data1 = df_counts1[['GO Term', 'Percentage']].copy()
    data1['Percentage'] = pd.to_numeric(data1['Percentage'], errors='coerce')
    data1.dropna(subset=['Percentage'], inplace=True)

    # **Apply cutoff to exclude values under 15**
    data1 = data1[data1['Percentage'] >= 15]

    data1.sort_values(by='Percentage', ascending=False, inplace=True)
    data1.reset_index(drop=True, inplace=True)
    data1_pivot = data1.set_index('GO Term')

    # Prepare and sort data for the second DataFrame
    data2 = df_counts2[['GO Term', 'Percentage']].copy()
    data2['Percentage'] = pd.to_numeric(data2['Percentage'], errors='coerce')
    data2.dropna(subset=['Percentage'], inplace=True)

    # **Apply cutoff to exclude values under 15**
    data2 = data2[data2['Percentage'] >= 15]

    data2.sort_values(by='Percentage', ascending=False, inplace=True)
    data2.reset_index(drop=True, inplace=True)
    data2_pivot = data2.set_index('GO Term')

    # Calculate global min and max percentages (after cutoff)
    all_percentages = pd.concat([data1['Percentage'], data2['Percentage']])
    min_percentage = all_percentages.min()
    max_percentage = all_percentages.max()

    # Number of rows in each heatmap
    num_rows1 = data1_pivot.shape[0]
    num_rows2 = data2_pivot.shape[0]

    # Define cell sizes
    cell_height = 0.5  # Inches per cell
    cell_width = 2     # Inches per cell

    # Calculate heights of each heatmap
    height1 = num_rows1 * cell_height
    height2 = num_rows2 * cell_height

    # Determine the maximum height for the figure
    fig_height = max(height1, height2) + 1  # Add 1 inch padding

    # Figure width
    fig_width = cell_width * 2 + 2  # Two heatmaps and padding

    # Create the figure
    fig = plt.figure(figsize=(fig_width, fig_height))

    # Left heatmap axes (ax0)
    left_ax0 = 0.1  # 10% from the left
    bottom_ax0 = (fig_height - height1) / fig_height / 2
    width_ax0 = 0.35  # 35% width
    height_ax0 = height1 / fig_height

    ax0 = fig.add_axes([left_ax0, bottom_ax0, width_ax0, height_ax0])

    # Right heatmap axes (ax1)
    left_ax1 = 0.55  # 55% from the left
    bottom_ax1 = (fig_height - height2) / fig_height / 2
    width_ax1 = 0.35  # 35% width
    height_ax1 = height2 / fig_height

    ax1 = fig.add_axes([left_ax1, bottom_ax1, width_ax1, height_ax1])

    # First heatmap
    sns.heatmap(
        data1_pivot,
        ax=ax0,
        annot=True,
        fmt=".2f",
        cmap="Blues",
        cbar=False,
        vmin=min_percentage,
        vmax=max_percentage,
        yticklabels=True,
        xticklabels=False,
        linewidths=0.5,
        linecolor='white'
    )

    ax0.set_title('Query', fontsize=12)
    ax0.set_xlabel('')
    ax0.set_ylabel('')
    ax0.set_xticks([])
    ax0.set_yticklabels(ax0.get_yticklabels(), rotation=0, fontsize=10)

    # Second heatmap
    sns.heatmap(
        data2_pivot,
        ax=ax1,
        annot=True,
        fmt=".2f",
        cmap="Blues",
        cbar=False,
        vmin=min_percentage,
        vmax=max_percentage,
        yticklabels=True,
        xticklabels=False,
        linewidths=0.5,
        linecolor='white'
    )
    ax1.set_title('Subject', fontsize=12)
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    ax1.set_xticks([])

    # Move y-axis labels to the right for the second heatmap
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0, fontsize=10)

    # Remove unnecessary labels
    for ax in [ax0, ax1]:
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])

    # Save and show the figure
    plt.savefig('combined_heatmaps.png', dpi=300, bbox_inches='tight')
    plt.show()

    print("Heatmap generation complete.")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Download and process genomic data with Prodigal and BLAST")
    # Arguments from the first script
    parser.add_argument('-i', '--input', required=True, help='Path to the input CSV file')
    parser.add_argument('-p', '--prodigal', required=True, help='Path to the prodigal executable')
    # Arguments for the genomes folder
    parser.add_argument("-d", "--taxid", required=True, help="Taxonomy ID for the organism (e.g., 194 for Campylobacter)")
    parser.add_argument("-g", "--go", help="Gene Ontology (GO) term (without 'GO:') to filter the TSV by (e.g., 0005975)")
    parser.add_argument('--filtering-go', help='GO term to be used later (without "GO:" prefix).')
    # Arguments for the process folder
    parser.add_argument('--process-folder', help='Path to a folder containing genome files to process starting from linearize_and_process_fasta')
    parser.add_argument('--process-folder-taxid', help='Taxonomy ID for the organism in the process folder')
    parser.add_argument('--process-folder-go', help='GO term for the process folder')
    # Shared arguments
    parser.add_argument("-m", "--makeblastdb", required=True, help="Path to the makeblastdb command")
    parser.add_argument("-b", "--blastx", required=True, help="Path to the blastx command")
    parser.add_argument("-t", "--threads", type=int, default=24, help="Number of threads to use for blastx (default: 24)")
    # Local file usage
    parser.add_argument('--uniprot-fasta', help='Path to a local UniProt FASTA file to use instead of downloading.')
    parser.add_argument('--uniprot-tsv', help='Path to a local UniProt TSV file to use instead of downloading.')
    parser.add_argument('--obo-file', help='Path to a local OBO file to use instead of downloading.')
    parser.add_argument('--process-folder-uniprot-fasta', help='Path to a local UniProt FASTA file for the process folder to use instead of downloading.')
    parser.add_argument('--process-folder-uniprot-tsv', help='Path to a local UniProt TSV file for the process folder to use instead of downloading.')

    args = parser.parse_args()

    # Validate threads
    if args.threads < 1:
        print("Error: Number of threads must be at least 1.")
        sys.exit(1)

    # Ensure filtering_go starts with 'GO:' if provided
    if args.filtering_go:
        if not args.filtering_go.startswith('GO:'):
            args.filtering_go = 'GO:' + args.filtering_go

    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

    # Step 1: Run download_and_rename from the first script
    download_and_rename(args.input, args.prodigal)

    # Step 2: Process the fasta files (unzipping, linearizing, running Prodigal) in 'genomes' folder
    process_fasta_files('genomes', args.prodigal)

    # If --process-folder is specified, process that folder separately
    if args.process_folder_uniprot_fasta and args.process_folder_uniprot_tsv:
        if not os.path.exists(args.process_folder_uniprot_fasta):
            print(f"Error: The UniProt FASTA file '{args.process_folder_uniprot_fasta}' does not exist.")
            sys.exit(1)
        if not os.path.exists(args.process_folder_uniprot_tsv):
            print(f"Error: The UniProt TSV file '{args.process_folder_uniprot_tsv}' does not exist.")
            sys.exit(1)
        fasta_file_process = args.process_folder_uniprot_fasta
        tsv_file_process = args.process_folder_uniprot_tsv
    else:
        fasta_file_process, tsv_file_process = download_uniprot_data(
            args.process_folder_taxid, args.process_folder_go, output_dir=args.process_folder
        )

    # Step 3: Download UniProt data for 'genomes' folder or use local files
    if args.uniprot_fasta and args.uniprot_tsv:
        if not os.path.exists(args.uniprot_fasta):
            print(f"Error: The UniProt FASTA file '{args.uniprot_fasta}' does not exist.")
            sys.exit(1)
        if not os.path.exists(args.uniprot_tsv):
            print(f"Error: The UniProt TSV file '{args.uniprot_tsv}' does not exist.")
            sys.exit(1)
        fasta_file_genomes = args.uniprot_fasta
        tsv_file_genomes = args.uniprot_tsv
    else:
        fasta_file_genomes, tsv_file_genomes = download_uniprot_data(args.taxid, args.go, output_dir='genomes')

    # Step 4: Run BLAST on 'genomes' folder
    protein_db_fasta_genomes = os.path.abspath(fasta_file_genomes)  # Use the downloaded FASTA as the protein DB
    run_blast('genomes', protein_db_fasta_genomes, args.makeblastdb, args.blastx, args.threads, output_dir='genomes')

    # Step 5: Process BLAST results for 'genomes' folder
    final_csv_genomes, final_csv_uniq_genomes = process_blast_results(tsv_file_genomes, 'genomes', output_dir='genomes')

    # Save the unique file results at the end as two separate variables
    unique_file_results = final_csv_genomes
    unique_file_results_uniq = final_csv_uniq_genomes

    # If --process-folder was specified, process it separately
    if args.process_folder:
        # Validate that taxid is provided for process-folder
        if not args.process_folder_taxid:
            print("Error: --process-folder-taxid is required when --process-folder is specified.")
            sys.exit(1)
        # Step 2b: Process the fasta files (unzipping, linearizing, running Prodigal) in the process folder
        process_fasta_files(args.process_folder, args.prodigal)
        # Step 3b: Download UniProt data for process-folder or use local files
        if args.process_folder_uniprot_fasta and args.process_folder_uniprot_tsv:
            fasta_file_process = args.process_folder_uniprot_fasta
            tsv_file_process = args.process_folder_uniprot_tsv
        else:
            fasta_file_process, tsv_file_process = download_uniprot_data(
                args.process_folder_taxid, args.process_folder_go, output_dir=args.process_folder
            )

        # Step 4b: Run BLAST on process-folder
        protein_db_fasta_process = os.path.abspath(fasta_file_process)
        run_blast(args.process_folder, protein_db_fasta_process, args.makeblastdb, args.blastx, args.threads, output_dir=args.process_folder)

        # Step 5b: Process BLAST results for process-folder
        final_csv_process, final_csv_uniq_process = process_blast_results(tsv_file_process, args.process_folder, output_dir=args.process_folder)

        # Save the unique file results at the end as two separate variables
        unique_file_results_process = final_csv_process
        unique_file_results_uniq_process = final_csv_uniq_process

    comparison(unique_file_results_uniq, unique_file_results_uniq_process, args.filtering_go, obo_file=args.obo_file)

    print("All tasks completed successfully.")

if __name__ == "__main__":
    main()

