import os
import subprocess
import tempfile
import pandas as pd

def check_database_exists(db_path):
    """Verifies that BLAST database index files exist before running."""
    extensions = [".nhr", ".nin", ".nsq"]
    for ext in extensions:
        if os.path.exists(db_path + ext):
            return True
    return False

def array_to_temp_fasta(kmers_list, fasta_path):
    """Converts a Python array/list of k-mers into FASTA format on the fly."""
    with open(fasta_path, "w") as f:
        for idx, kmer in enumerate(kmers_list, 1):
            kmer_clean = str(kmer).strip()
            if kmer_clean and kmer_clean != "nan":
                f.write(f">kmer_{idx}\n{kmer_clean}\n")

def run_blast_short(query_fasta, db_path, output_tsv):
    """Runs BLASTN for short k-mers."""
    blast_cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", db_path,
        "-task", "blastn-short",
        "-word_size", "7",
        "-evalue", "1000",
        "-perc_identity", "100",
        "-out", output_tsv,
        "-outfmt", "6 qseqid sseqid pident length mismatch evalue bitscore stitle"
    ]
    subprocess.run(blast_cmd, check=True)

def categorize_kmer(in_card, in_ecoli):
    """Applies the biological interpretation decision matrix."""
    if in_card and in_ecoli:
        return "Known E. coli Resistance", "Keep"
    elif in_card and not in_ecoli:
        return "Possible Horizontal Gene Transfer (HGT)", "Flag"
    elif not in_card and in_ecoli:
        return "Lineage / Non-canonical / Novel Candidate", "Keep"
    else:
        return "Unmapped Noise / Technical Artifact", "Discard"

def main():
    # Database Paths
    project_root = os.path.abspath(os.curdir)
    card_db = os.path.join(project_root, "card_database", "card_db")
    ecoli_db = os.path.join(project_root, "ecoli_database", "ecoli_ref_db")
    feature_dir = "output/feature_importances/piperacillin_tazobactam_feature_importance.txt"
    
    # Check that BLAST databases exist
    if not check_database_exists(card_db):
        raise FileNotFoundError(f"CARD database index missing at: {card_db}. Ensure Step 2 setup was completed.")
    if not check_database_exists(ecoli_db):
        raise FileNotFoundError(f"E. coli database index missing at: {ecoli_db}. Run Step 1 PowerShell setup first.")

    # Input: You can pass an array of k-mers directly or read from your feature importance CSV
    # Example input array (Replace or load your model's actual top k-mers here):
   
   # Initialize as an empty list to match your desired structure
    kmers_input: list[dict[str, any]] = []

    with open(feature_dir, "r", encoding="utf8", errors="ignore") as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            kmer, imp = parts
            try:
                # Create a dictionary for the current line and append it to the list
                kmers_input.append({
                    "Kmer": kmer, 
                    "Importance": round(float(imp),2)  # Convert the string to a decimal float
                })
            except ValueError:
                continue


    
    # If reading from a CSV file instead, uncomment the line below:
    # kmers_input = pd.read_csv("top_kmers.csv").to_dict('records')

    # CREATE A SELF-CLEANING TEMPORARY DIRECTORY
    # All intermediate FASTA and raw BLAST TSVs will be created inside here
    # and automatically DELETED when this block finishes!
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_query_fasta = os.path.join(temp_dir, "query_kmers.fasta")
        card_out_tsv = os.path.join(temp_dir, "card_hits.tsv")
        ecoli_out_tsv = os.path.join(temp_dir, "ecoli_hits.tsv")

        print("[1/4] Converting k-mers array to temporary FASTA buffer...")
        kmer_sequences = [item["Kmer"] for item in kmers_input]
        array_to_temp_fasta(kmer_sequences, temp_query_fasta)

        print("[2/4] Running BLAST against CARD database...")
        run_blast_short(temp_query_fasta, card_db, card_out_tsv)

        print("[3/4] Running BLAST against E. coli reference database...")
        run_blast_short(temp_query_fasta, ecoli_db, ecoli_out_tsv)

        print("[4/4] Parsing hits and categorizing features...")
        cols = ["kmer_id", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "stitle"]
        
        card_hits = pd.read_csv(card_out_tsv, sep="\t", names=cols) if os.path.exists(card_out_tsv) and os.path.getsize(card_out_tsv) > 0 else pd.DataFrame(columns=cols)
        ecoli_hits = pd.read_csv(ecoli_out_tsv, sep="\t", names=cols) if os.path.exists(ecoli_out_tsv) and os.path.getsize(ecoli_out_tsv) > 0 else pd.DataFrame(columns=cols)

        # Extract matching sequences (by comparing raw sequence matches)
        # Match k-mer string to BLAST subject descriptions
        results = []
        for item in kmers_input:
            kmer_seq = item["Kmer"]
            importance = item["Importance"]

            # Check presence in CARD hits
            card_matches = card_hits[card_hits["stitle"].str.contains(kmer_seq, case=False, na=False) | (card_hits["kmer_id"].notna())]
            in_card = len(card_hits) > 0 and any(kmer_seq in str(row) for row in card_hits["stitle"].tolist() + card_hits["sseqid"].tolist())
            
            # Simplified check: sequence exists in output
            in_card = len(card_hits) > 0
            in_ecoli = len(ecoli_hits) > 0

            # Match specifically for this k-mer based on query sequence
            # Since query IDs are kmer_1, kmer_2...
            kmer_idx = kmer_sequences.index(kmer_seq) + 1
            q_id = f"kmer_{kmer_idx}"

            kmer_card_rows = card_hits[card_hits["kmer_id"] == q_id]
            kmer_ecoli_rows = ecoli_hits[ecoli_hits["kmer_id"] == q_id]

            has_card_hit = len(kmer_card_rows) > 0
            has_ecoli_hit = len(kmer_ecoli_rows) > 0

            category, action = categorize_kmer(has_card_hit, has_ecoli_hit)

            gene_title = "Unmapped Region"
            if has_card_hit:
                raw_title = kmer_card_rows.iloc[0]["stitle"]
                gene_title = raw_title.split("|")[-1].strip() if "|" in raw_title else raw_title

            results.append({
                "Kmer_Feature": kmer_seq,
                "Importance_Score": importance,
                "In_CARD": has_card_hit,
                "In_Ecoli_Genome": has_ecoli_hit,
                "Mapped_Target": gene_title,
                "Biological_Category": category,
                "Pipeline_Action": action
            })

        # Save ONLY the final clean CSV output file
        final_df = pd.DataFrame(results)
        final_output_file = "output/biological_annotations/piperacillin_tazobactam_kmer_biological_interpretation_results.csv"
        final_df.to_csv(final_output_file, index=False)
        print(f"[+] Complete! Only the final results file was saved: '{final_output_file}'")

if __name__ == "__main__":
    main()