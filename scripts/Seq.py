import csv
import time
import urllib.request
import gzip
import shutil
import os
from Bio import Entrez

# Always tell NCBI who you are
Entrez.email = "lewistem8@gmail.com"

def download_assembly_sequences(accession_list):
    for assembly_id in accession_list:
        print(f"--- Processing Assembly: {assembly_id} ---")
        try:
            # 1. Get the Assembly Summary to find the FTP link
            handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
            record = Entrez.read(handle)
            handle.close()

            if not record['DocumentSummarySet']['DocumentSummary']:
                print(f"No summary found for {assembly_id}")
                continue

            # Get the FTP URL for the GenBank or RefSeq assembly
            summary = record['DocumentSummarySet']['DocumentSummary'][0]
            ftp_url = summary.get('FtpPath_RefSeq', '') or summary.get('FtpPath_GenBank', '')

            if not ftp_url:
                print(f"No FTP path found for {assembly_id}")
                continue

            # Construct the filename for the genomic FASTA
            # Format usually: [FTP_URL]/[Assembly_Name]_genomic.fna.gz
            assembly_name = ftp_url.split('/')[-1]
            fasta_filename = f"{assembly_name}_genomic.fna.gz"
            fasta_url = f"{ftp_url}/{fasta_filename}"
            
            output_filename = f"{assembly_id}.fasta"
            
            print(f"Downloading genome from: {fasta_url}")
            
            # 2. Download the file
            # We use urllib because Entrez efetch doesn't support large binary assembly downloads well
            urllib.request.urlretrieve(fasta_url, fasta_filename)
            
            # 3. Unzip the file
            print(f"Unzipping to {output_filename}...")
            with gzip.open(fasta_filename, 'rb') as f_in:
                with open(output_filename, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            # Cleanup: remove the .gz file
            os.remove(fasta_filename)

            print(f"Successfully saved genome: {output_filename}")

        except Exception as e:
            print(f"Error downloading {assembly_id}: {e}")

if __name__ == "__main__":
    # Example: List of Assembly Accessions (GCF or GCA)
    my_assemblies = ["GCF_001642285.1", "GCA_000219515.3", "GCA_000219515.2"] # Replace with your list
    download_assembly_sequences(my_assemblies)
