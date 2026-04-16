"""Build a wide, memory-efficient sparse k-mer matrix from the long-format CSV.

This avoids using pandas.pivot_table (dense) which can exhaust RAM.

Output:
  - kmer_matrix_sparse.npz : SciPy CSR matrix (rows = genomes, cols = selected k-mers)
  - kmer_matrix_genome_index.csv : mapping from row index -> GenomeID
  - kmer_matrix_kmer_index.csv   : mapping from column index -> kmer

Requires: scipy, pandas, numpy
Run from project root:
  python scripts/wide_format_sparse.py
"""

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, save_npz


def build_kmer_vocab(
    long_csv_path: Path,
    min_genome_freq: int = 5,
    top_k: int | None = None,
    chunksize: int = 1_000_000,
) -> dict:
    """Build a k-mer vocabulary from the long-format CSV.

    We count in how many rows each k-mer appears. Since each row is one
    (GenomeID, kmer) pair, this is approximately the number of genomes
    that contain that k-mer.

    Parameters
    ----------
    long_csv_path : Path
        Path to kmer-data-long.csv (GenomeID, kmer, count).
    min_genome_freq : int
        Keep only k-mers that appear in at least this many rows (genomes).
    top_k : int or None
        If not None, keep only the top_k most frequent k-mers after filtering.
    chunksize : int
        Number of rows per chunk when streaming the CSV.
    """

    print("Building k-mer vocabulary (first pass over long CSV)...")
    freq = Counter()

    for chunk in pd.read_csv(long_csv_path, chunksize=chunksize):
        vc = chunk["kmer"].value_counts()
        for k, v in vc.items():
            freq[k] += int(v)

    # Apply minimum genome frequency filter
    items = [(k, c) for k, c in freq.items() if c >= min_genome_freq]
    if not items:
        raise ValueError("No k-mers passed the min_genome_freq filter; try lowering it.")

    # Sort by frequency (descending)
    items.sort(key=lambda x: x[1], reverse=True)

    # Optionally keep only top_k
    if top_k is not None:
        items = items[:top_k]

    vocab = {kmer: idx for idx, (kmer, _) in enumerate(items)}
    print(f"Selected {len(vocab)} k-mers for the vocabulary.")
    return vocab


def build_genome_index(
    long_csv_path: Path,
    chunksize: int = 1_000_000,
) -> dict:
    """Assign a row index to each GenomeID found in the long CSV."""

    print("Building genome index (second pass over long CSV)...")
    genome_to_row: dict[str, int] = {}

    for chunk in pd.read_csv(long_csv_path, chunksize=chunksize):
        for gid in chunk["GenomeID"].astype(str).unique():
            if gid not in genome_to_row:
                genome_to_row[gid] = len(genome_to_row)

    print(f"Found {len(genome_to_row)} genomes.")
    return genome_to_row


def build_sparse_matrix(
    long_csv_path: Path,
    vocab: dict,
    genome_to_row: dict,
    chunksize: int = 1_000_000,
) -> coo_matrix:
    """Build a sparse genome x k-mer matrix from the long CSV.

    Only k-mers in `vocab` are kept.
    """

    print("Building sparse matrix (third pass over long CSV)...")
    rows: list[int] = []
    cols: list[int] = []
    data: list[int] = []

    for chunk in pd.read_csv(long_csv_path, chunksize=chunksize):
        # Ensure correct dtypes
        chunk["GenomeID"] = chunk["GenomeID"].astype(str)

        # Filter to k-mers we care about
        mask = chunk["kmer"].isin(vocab)
        if not mask.any():
            continue
        sub = chunk.loc[mask, ["GenomeID", "kmer", "count"]]

        # Map to indices
        row_idx = sub["GenomeID"].map(genome_to_row).to_numpy()
        col_idx = sub["kmer"].map(vocab).to_numpy()
        counts = sub["count"].to_numpy()

        rows.extend(row_idx.tolist())
        cols.extend(col_idx.tolist())
        data.extend(counts.tolist())

    n_rows = len(genome_to_row)
    n_cols = len(vocab)
    print(f"Constructing COO matrix with shape ({n_rows}, {n_cols}) and {len(data)} nonzeros...")

    coo = coo_matrix((np.array(data, dtype=np.int32), (np.array(rows), np.array(cols))),
                     shape=(n_rows, n_cols))
    return coo.tocsr()


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    long_csv_path = project_root / "kmer-data-long.csv"

    if not long_csv_path.exists():
        raise FileNotFoundError(f"Long-format file not found: {long_csv_path}")

    # Hyperparameters: tune these based on available RAM and desired feature size
    min_genome_freq = 5   # keep k-mers that appear in at least this many genomes
    top_k = 50_000        # or None to keep all after filtering

    # 1) Build k-mer vocabulary
    vocab = build_kmer_vocab(
        long_csv_path=long_csv_path,
        min_genome_freq=min_genome_freq,
        top_k=top_k,
    )

    # 2) Build genome index
    genome_to_row = build_genome_index(long_csv_path=long_csv_path)

    # 3) Build sparse matrix
    X_csr = build_sparse_matrix(
        long_csv_path=long_csv_path,
        vocab=vocab,
        genome_to_row=genome_to_row,
    )

    # 4) Save matrix and index mappings
    out_matrix = project_root / "kmer_matrix_sparse.npz"
    out_genome_index = project_root / "kmer_matrix_genome_index.csv"
    out_kmer_index = project_root / "kmer_matrix_kmer_index.csv"

    print(f"Saving sparse matrix to {out_matrix}...")
    save_npz(out_matrix, X_csr)

    # Save genome index: row -> GenomeID
    genome_index_df = (
        pd.Series(genome_to_row)
        .reset_index()
        .rename(columns={"index": "GenomeID", 0: "row"})
    )
    genome_index_df.to_csv(out_genome_index, index=False)

    # Save k-mer index: col -> kmer
    kmer_index_df = (
        pd.Series(vocab)
        .reset_index()
        .rename(columns={"index": "kmer", 0: "col"})
    )
    kmer_index_df.to_csv(out_kmer_index, index=False)

    print("Done. Outputs:")
    print(f"  Sparse matrix: {out_matrix}")
    print(f"  Genome index: {out_genome_index}")
    print(f"  K-mer index:  {out_kmer_index}")


if __name__ == "__main__":
    main()
