import csv


def write_kmers_long_format(genome_ids_path: str,
                            output_csv_path: str,
                            not_counted_path: str = "data/not_counted_ids.txt") -> None:
    """Stream all k-mers for all genomes into a single long-format CSV.

    Output columns: GenomeID, kmer, count
    This avoids building a huge in-memory wide DataFrame that can crash the PC.
    """

    # Read genome IDs to process
    with open(genome_ids_path, "r") as f:
        genome_ids = [line.strip() for line in f if line.strip()]

    # Prepare not-counted log
    not_counted_file = open(not_counted_path, "w")

    # Open output CSV once and stream rows
    with open(output_csv_path, "w", newline="") as out_f:
        writer = csv.writer(out_f)
        writer.writerow(["GenomeID", "kmer", "count"])  # header

        for gid in genome_ids[:1000]:
            dump_file = f"output/{gid}_db_kmers.txt"
            try:
                with open(dump_file, "r") as df:
                    for line in df:
                        parts = line.strip().split()
                        if len(parts) != 2:
                            continue
                        kmer, cnt = parts
                        writer.writerow([gid, kmer, int(cnt)])
            except FileNotFoundError:
                not_counted_file.write(f"{gid}\n")
                continue
    not_counted_file.close()
    print(f"Done. Long-format k-mer data written to {output_csv_path}")


if __name__ == "__main__":
    write_kmers_long_format(
        genome_ids_path="data/genome_ids.txt",
        output_csv_path="kmer-data-long.csv",
    )