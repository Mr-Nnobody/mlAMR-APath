import os
import time
import requests


def download_bvbrc_genomes(genome_ids):
    os.makedirs('data/fasta_files', exist_ok=True)
    
    for gid in genome_ids:
        print(f"Downloading Genome ID: {gid}...")
        url = f"https://www.bv-brc.org/api/genome_sequence/?eq(genome_id,{gid})&http_accept=application/dna+fasta"
        
        response = requests.get(url)
        if response.status_code == 200 and len(response.text) > 100:
            with open(f"data/fasta_files/{gid}.fasta", "w") as f:
                f.write(response.text)
            print(f"Success: {gid}.fasta")
        else:
            # log failed genome IDs for review
            with open("data/failed_genome_ids.txt", "a") as f:
                    f.write(gid + "\n")
            print(f"Failed to find sequence for {gid}")
        time.sleep(0.5) # respect the API 
        

# Example testing with a Genome ID you know doesn't have an NCBI Assembly
with open("data/genome_ids.txt", "r") as f:
    genome_ids = [line.strip() for line in f]
gene= ['562.97996', '562.99422', '562.7338', '562.100175', '562.99496', '562.99937', '562.99739', '562.99389', '562.98309', '562.9977', '562.10006', '562.144216', '562.98868', '562.22727']
download_bvbrc_genomes(gene)
print("Done downloading genomes from BV-BRC. Check data/fasta_files/ for results and data/failed_genome_ids.txt for any missing sequences.")