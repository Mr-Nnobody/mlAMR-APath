import csv
import time
from http.client import IncompleteRead
from Bio import Entrez, SeqIO

# Always tell NCBI who you are
Entrez.email = "lewistem8@gmail.com"

def download_assembly_sequences(accession_list):
    for assembly_id in accession_list:
        print(f"--- Processing Assembly: {assembly_id} ---")
        try:
            # 1. Search for the assembly to get the UID
            search_handle = Entrez.esearch(db="assembly", term=assembly_id)
            search_record = Entrez.read(search_handle)
            search_handle.close()

            if not search_record['IdList']:
                print(f"No results found for {assembly_id}")
                continue

            assembly_uid = search_record['IdList'][0]

            # 2. Get links to the Nucleotide (nuccore) database
            link_handle = Entrez.elink(dbfrom="assembly", db="nuccore", id=assembly_uid)
            link_record = Entrez.read(link_handle)
            link_handle.close()

            # Extract list of IDs for all components (chromosomes/scaffolds)
            nucleotide_ids = [
                link['Id'] for link in link_record[0]['LinkSetDb'][0]['Link']
            ]

            # 3. Fetch and save sequences to a FASTA file
            output_filename = f"{assembly_id}_full.fasta"
            print(f"Downloading {len(nucleotide_ids)} sequences to {output_filename}...")

            with open(output_filename, "w") as out_f:
                # Fetch in batches to prevent server timeouts
                for i in range(0, len(nucleotide_ids), 10):
                    batch = nucleotide_ids[i : i + 10]
                    fetch_handle = Entrez.efetch(
                        db="nuccore", 
                        id=batch, 
                        rettype="fasta", 
                        retmode="text"
                    )
                    out_f.write(fetch_handle.read())
                    fetch_handle.close()
                    time.sleep(1) # Be kind to NCBI servers

            print(f"Successfully saved {assembly_id}")

        except Exception as e:
            print(f"Error downloading {assembly_id}: {e}")

if __name__ == "__main__":
    # Example: List of Assembly Accessions (GCF or GCA)
    my_assemblies = ["GCF_001642285.1", "GCA_000219515.3", "GCA_000219515.2"] # Replace with your list
    download_assembly_sequences(my_assemblies)
