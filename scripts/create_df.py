import csv
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def write_kmers_long_format(genome_ids_path: str,
                            output_csv_path: str,
                            not_counted_path: str = "data/not-in-df_ids.txt") -> None:
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

        total = len(genome_ids)
        for idx, gid in enumerate(genome_ids, start=1):
            print(f"Processing {idx}/{total} genome IDs", end="\r", flush=True)
            dump_file = f"output/counted_kmers/{gid}_db_kmers.txt"
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
            print(" " * 60, end="\r", flush=True)
    not_counted_file.close()
    print(f"Done. Long-format k-mer data written to {output_csv_path}")


def write_kmers_long_parquet(
    genome_ids_path: str,
    output_parquet_path: str,
    not_counted_path: str = "data/not-in-df_ids.txt",
    chunk_rows: int = 500_000,
) -> None:
    """Stream all k-mers into a long-format Parquet file.

    Output columns: GenomeID, kmer, count
    """

    with open(genome_ids_path, "r") as f:
        genome_ids = [line.strip() for line in f if line.strip()]

    not_counted_file = open(not_counted_path, "w")

    out_path = Path(output_parquet_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    schema = pa.schema(
        [
            ("GenomeID", pa.string()),
            ("kmer", pa.string()),
            ("count", pa.int32()),
        ]
    )

    writer = pq.ParquetWriter(out_path, schema=schema, compression="zstd")

    total = len(genome_ids)
    rows_buffer = {"GenomeID": [], "kmer": [], "count": []}

    def flush_buffer() -> None:
        if not rows_buffer["GenomeID"]:
            return
        table = pa.Table.from_pydict(rows_buffer, schema=schema)
        writer.write_table(table)
        rows_buffer["GenomeID"].clear()
        rows_buffer["kmer"].clear()
        rows_buffer["count"].clear()

    try:
        for idx, gid in enumerate(genome_ids, start=1):
            print(f"Processing {idx}/{total} genome IDs", end="\r", flush=True)
            dump_file = f"output/counted_kmers/{gid}_db_kmers.txt"
            try:
                with open(dump_file, "r") as df:
                    for line in df:
                        parts = line.strip().split()
                        if len(parts) != 2:
                            continue
                        kmer, cnt = parts
                        rows_buffer["GenomeID"].append(gid)
                        rows_buffer["kmer"].append(kmer)
                        rows_buffer["count"].append(int(cnt))
                        if len(rows_buffer["GenomeID"]) >= chunk_rows:
                            flush_buffer()
            except FileNotFoundError:
                not_counted_file.write(f"{gid}\n")
                continue

        flush_buffer()
        print(" " * 60, end="\r", flush=True)
    finally:
        writer.close()
        not_counted_file.close()

    print(f"Done. Long-format k-mer data written to {output_parquet_path}")


if __name__ == "__main__":
    write_kmers_long_parquet(
        genome_ids_path="data/genome_ids.txt",
        output_parquet_path="kmer-data-long.parquet",
    )