import csv
import time
import urllib.request
import gzip
import shutil
import os
from Bio import Entrez

# Always tell NCBI who you are
Entrez.email = "lewistem8@gmail.com"

def download_sra_linked_assemblies(sra_list):
    # Ensure the output directory exists
    os.makedirs('data/fasta_files', exist_ok=True)
    
    for sra_id in sra_list:
        print(f"--- Processing SRA Accession: {sra_id} ---")
        try:
            # 1. Search to get the internal SRA UID
            search_handle = Entrez.esearch(db="sra", term=sra_id)
            search_record = Entrez.read(search_handle)
            search_handle.close()

            if not search_record['IdList']:
                print(f"No SRA record found for {sra_id}")
                continue

            sra_uid = search_record['IdList'][0]

            # 2. Route SRA -> BioSample -> Assembly
            # a. SRA to BioSample
            bs_link_handle = Entrez.elink(dbfrom="sra", db="biosample", id=sra_uid)
            bs_link_record = Entrez.read(bs_link_handle)
            bs_link_handle.close()

            biosample_uids = []
            if bs_link_record and bs_link_record[0].get('LinkSetDb'):
                for linkset in bs_link_record[0]['LinkSetDb']:
                    if linkset['LinkName'] == 'sra_biosample' or linkset['DbTo'] == 'biosample':
                        biosample_uids.extend([link['Id'] for link in linkset['Link']])
            
            if not biosample_uids:
                print(f"No linked BioSample found for SRA {sra_id}")
                continue

            # b. BioSample to Assembly
            assembly_uids = []
            for bs_uid in biosample_uids:
                asm_link_handle = Entrez.elink(dbfrom="biosample", db="assembly", id=bs_uid)
                asm_link_record = Entrez.read(asm_link_handle)
                asm_link_handle.close()

                if asm_link_record and asm_link_record[0].get('LinkSetDb'):
                    for linkset in asm_link_record[0]['LinkSetDb']:
                        if linkset['LinkName'] == 'biosample_assembly' or linkset['DbTo'] == 'assembly':
                            assembly_uids.extend([link['Id'] for link in linkset['Link']])
            
            if not assembly_uids:
                print(f"No linked Assembly found for BioSample of SRA {sra_id}")
                continue
                
            assembly_uid = assembly_uids[0] # Take the first linked assembly if multiple exist
            print(f"Found linked Assembly UID: {assembly_uid}")

            # 3. Get the Assembly Summary using the UID
            handle = Entrez.esummary(db="assembly", id=assembly_uid, report="full")
            record = Entrez.read(handle, validate=False)
            handle.close()

            if not record['DocumentSummarySet']['DocumentSummary']:
                print(f"No assembly summary found for linked UID {assembly_uid}")
                continue

            summary = record['DocumentSummarySet']['DocumentSummary'][0]
            ftp_url = summary.get('FtpPath_RefSeq', '') or summary.get('FtpPath_GenBank', '')

            if not ftp_url:
                print(f"No FTP path found for linked Assembly of {sra_id}")
                continue

            # Switch to HTTPS for better reliability
            ftp_url = ftp_url.replace("ftp://", "https://")

            # Construct the filename for the genomic FASTA
            assembly_name = ftp_url.split('/')[-1]
            fasta_filename = f"{assembly_name}_genomic.fna.gz"
            fasta_url = f"{ftp_url}/{fasta_filename}"
            
            output_filename = f"{sra_id}_linked.fasta"
            
            print(f"Downloading genome from: {fasta_url}")
            
            # 4. Download the file with retries
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
            
            # 5. Unzip the file
            print(f"Unzipping to {output_filename}...")
            with gzip.open(fasta_filename, 'rb') as f_in:
                with open(f'data/fasta_files/{output_filename}', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            # Cleanup
            os.remove(fasta_filename)

            print(f"Successfully saved linked assembly to: {output_filename}")

        except Exception as e:
            print(f"Error downloading for SRA {sra_id}: {e}")

if __name__ == "__main__":
    # Example: List of SRA Accessions (SRR, ERR, DRR)
    my_sra_accessions = ['SRS4578163'] # Replace with your actual SRA accessions
    download_sra_linked_assemblies(my_sra_accessions)
