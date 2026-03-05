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
            # 1. Search to get the internal Assembly UID
            # esummary requires an internal integer ID (UID), not the string Accession (GCF_...)
            search_handle = Entrez.esearch(db="assembly", term=assembly_id)
            search_record = Entrez.read(search_handle)
            search_handle.close()

            if not search_record['IdList']:
                print(f"No results found for {assembly_id}")
                continue

            assembly_uid = search_record['IdList'][0]

            # 2. Get the Assembly Summary using the UID
            handle = Entrez.esummary(db="assembly", id=assembly_uid, report="full")
            record = Entrez.read(handle, validate=False)
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

            # Switch to HTTPS for better reliability
            ftp_url = ftp_url.replace("ftp://", "https://")

            # Construct the filename for the genomic FASTA
            # Format usually: [FTP_URL]/[Assembly_Name]_genomic.fna.gz
            assembly_name = ftp_url.split('/')[-1]
            fasta_filename = f"{assembly_name}_genomic.fna.gz"
            fasta_url = f"{ftp_url}/{fasta_filename}"
            
            output_filename = f"{assembly_id}.fasta"
            
            print(f"Downloading genome from: {fasta_url}")
            
            # 2. Download the file with retries
            # We use urllib because Entrez efetch doesn't support large binary assembly downloads well
            max_retries = 3
            for attempt in range(max_retries):
                try:
                    urllib.request.urlretrieve(fasta_url, fasta_filename)
                    break # Success!
                except Exception as download_error:
                    if attempt == max_retries - 1:
                        raise download_error
                    print(f"Download incomplete, retrying ({attempt + 1}/{max_retries})...")
                    time.sleep(2)  # Wait a bit before retrying
            
            # 3. Unzip the file
            print(f"Unzipping to {output_filename}...")
            with gzip.open(fasta_filename, 'rb') as f_in:
                with open(f'data/fasta_files/{output_filename}', 'wb') as f_out: #TODO check this line
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
