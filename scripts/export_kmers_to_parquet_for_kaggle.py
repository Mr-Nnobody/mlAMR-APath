#!/usr/bin/env python3
"""Export k-mer features + phenotype to a wide Parquet for Kaggle/AutoScientist.

This script builds a prevalence-filtered vocabulary from per-genome k-mer dump
files, constructs a sparse presence matrix, applies chi2 selection to pick
the top `--select-k` k-mers, densifies only the selected features and writes
a wide Parquet file with a `GenomeID` column and one column per selected k-mer.

Designed for use on Kaggle where memory may be limited; keep `--select-k`
moderate (e.g., 5000) to avoid OOM when densifying.

Example:
  python scripts/export_kmers_to_parquet_for_kaggle.py \
    --dump-dir /kaggle/input/counted_kmers \
    --phenotype /kaggle/input/ampicillin_phenotype.csv \
    --out-parquet /kaggle/working/ampicillin_kmers_5k.parquet \
    --select-k 5000
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import sys
import json

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.feature_selection import chi2, SelectKBest


def iter_genome_kmers(dump_path: Path) -> dict[str, int]:
    kmers: dict[str, int] = {}
    with dump_path.open("r", encoding="utf8", errors="ignore") as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            kmer, cnt = parts
            try:
                kmers[kmer] = int(cnt)
            except ValueError:
                continue
    return kmers


def build_prevalence(dump_dir: Path, genome_ids: list[str]) -> Counter:
    prevalence: Counter = Counter()
    for gid in genome_ids:
        dump_path = dump_dir / f"{gid}_db_kmers.txt"
        if not dump_path.exists():
            continue
        kmers = iter_genome_kmers(dump_path)
        prevalence.update(kmers.keys())
    return prevalence


def select_vocab_by_prevalence(prevalence: Counter, n_genomes: int, min_frac: float, max_frac: float) -> list[str]:
    min_count = int(np.ceil(min_frac * n_genomes))
    max_count = int(np.floor(max_frac * n_genomes))
    vocab = [k for k, c in prevalence.items() if min_count <= c <= max_count]
    return vocab


def build_sparse_matrix(dump_dir: Path, genome_ids: list[str], vocab: list[str], binary: bool = True, chunk_size: int = 100) -> sparse.csr_matrix:
    vocab_index = {kmer: idx for idx, kmer in enumerate(vocab)}
    n_features = len(vocab)
    blocks = []
    dtype = np.int8 if binary else np.int32

    for start in range(0, len(genome_ids), chunk_size):
        chunk_ids = genome_ids[start:start + chunk_size]
        rows: list[int] = []
        cols: list[int] = []
        data: list[int] = []
        for row_idx, gid in enumerate(chunk_ids):
            dump_path = dump_dir / f"{gid}_db_kmers.txt"
            if not dump_path.exists():
                continue
            kmers = iter_genome_kmers(dump_path)
            for kmer in kmers.keys():
                col_idx = vocab_index.get(kmer)
                if col_idx is None:
                    continue
                rows.append(row_idx)
                cols.append(col_idx)
                data.append(1 if binary else int(kmers.get(kmer, 0)))
        if rows:
            block = sparse.csr_matrix((data, (rows, cols)), shape=(len(chunk_ids), n_features), dtype=dtype)
        else:
            block = sparse.csr_matrix((len(chunk_ids), n_features), dtype=dtype)
        blocks.append(block)
    if not blocks:
        return sparse.csr_matrix((0, n_features), dtype=dtype)
    return sparse.vstack(blocks, format='csr')


def load_labels(labels_path: Path, gid_col_candidates=('Genome ID', 'GenomeID', 'genome id', 'genomeid'), pheno_col_candidates=('phenotype', 'phenotype_label', 'pheno')) -> pd.Series:
    df = pd.read_csv(labels_path)
    cols_lower = {c.strip().lower(): c for c in df.columns}
    gid_col = None
    pheno_col = None
    for c in cols_lower:
        if c in [g.lower() for g in gid_col_candidates]:
            gid_col = cols_lower[c]
        if c in [p.lower() for p in pheno_col_candidates]:
            pheno_col = cols_lower[c]
    if gid_col is None or pheno_col is None:
        raise ValueError(f"Could not find GenomeID/phenotype columns in {labels_path}. Found: {list(df.columns)}")
    pheno = df.set_index(gid_col)[pheno_col]
    pheno = pd.to_numeric(pheno, errors='coerce')
    return pheno


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Export selected k-mer features + phenotype to a wide Parquet for Kaggle")
    p.add_argument('--dump-dir', required=True, help='Directory containing {GenomeID}_db_kmers.txt files')
    p.add_argument('--phenotype', required=True, help='CSV with genome phenotype labels')
    p.add_argument('--out-parquet', required=True, help='Output Parquet path (wide format)')
    p.add_argument('--out-vocab', default=None, help='Path to save selected k-mers (txt)')
    p.add_argument('--min-frac', type=float, default=0.02, help='Min prevalence fraction')
    p.add_argument('--max-frac', type=float, default=0.95, help='Max prevalence fraction')
    p.add_argument('--max-vocab', type=int, default=200000, help='Cap raw vocab to this many top-by-prevalence k-mers before selection')
    p.add_argument('--select-k', type=int, default=5000, help='Number of top k-mers to select with chi2 (final feature columns)')
    p.add_argument('--chunk-size', type=int, default=100, help='Chunk size when building sparse matrix')
    p.add_argument('--compression', default='zstd', help='Parquet compression (zstd/snappy/gzip)')
    args = p.parse_args(argv)

    dump_dir = Path(args.dump_dir)
    phen_path = Path(args.phenotype)
    out_parquet = Path(args.out_parquet)
    out_parquet.parent.mkdir(parents=True, exist_ok=True)

    labels = load_labels(phen_path)
    labels.index = labels.index.astype(str)
    # genomes that have labels and a dump file
    available = [p.stem.replace('_db_kmers', '') for p in dump_dir.glob('*_db_kmers.txt')]
    genome_ids = [g for g in labels.index.astype(str) if g in available]
    if not genome_ids:
        raise RuntimeError('No genomes found with both labels and k-mer dumps')
    print(f'Found {len(genome_ids)} genomes with dumps and labels')

    # prevalence and raw vocab
    prevalence = build_prevalence(dump_dir, genome_ids)
    print(f'Unique k-mers observed: {len(prevalence):,d}')

    vocab = select_vocab_by_prevalence(prevalence, n_genomes=len(genome_ids), min_frac=args.min_frac, max_frac=args.max_frac)
    print(f'Vocab after prevalence filter: {len(vocab):,d}')

    if len(vocab) > args.max_vocab:
        # keep top by prevalence
        vocab = [k for k, _ in prevalence.most_common(args.max_vocab) if k in set(vocab)]
        print(f'Vocab capped to top {args.max_vocab:,d} -> {len(vocab):,d}')

    # Build sparse matrix (rows in genome_ids order)
    X = build_sparse_matrix(dump_dir, genome_ids, vocab, binary=True, chunk_size=args.chunk_size)
    print('Built sparse matrix:', X.shape)

    # Build label array aligned to genome_ids
    y = np.array([int(labels.loc[g]) for g in genome_ids], dtype=int)

    # SelectKBest with chi2 (supports sparse input)
    k = min(args.select_k, X.shape[1])
    if k <= 0:
        raise ValueError('select-k must be > 0')
    print(f'Selecting top {k} k-mers with chi2')
    selector = SelectKBest(score_func=chi2, k=k)
    selector.fit(X, y)
    X_sel = selector.transform(X)
    selected_idx = selector.get_support(indices=True)
    selected_kmers = [vocab[i] for i in selected_idx]
    print('Selected k-mers:', len(selected_kmers))

    if args.out_vocab:
        Path(args.out_vocab).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out_vocab, 'w', encoding='utf8') as fh:
            fh.write('\n'.join(selected_kmers))
        print('Wrote selected vocab to', args.out_vocab)

    # Safety: estimate memory to densify
    n_rows = X_sel.shape[0]
    n_cols = X_sel.shape[1]
    est_bytes = n_rows * n_cols * 8  # float64
    est_gb = est_bytes / (1024**3)
    print(f'Estimate densified array size: {est_gb:.2f} GB')
    if est_gb > 16:
        print('Warning: densifying >16GB will likely OOM on Kaggle. Consider reducing --select-k')

    # Densify selected features and write wide parquet
    print('Densifying selected features...')
    X_dense = X_sel.toarray()
    df = pd.DataFrame(X_dense, columns=selected_kmers)
    df.insert(0, 'GenomeID', genome_ids)
    # attach phenotype column name 'phenotype'
    df['phenotype'] = [int(labels.loc[g]) for g in genome_ids]

    print('Writing Parquet to', out_parquet)
    df.to_parquet(out_parquet, index=False, compression=args.compression)
    # save a small metadata file
    meta = {
        'n_genomes': int(n_rows),
        'n_features': int(n_cols),
        'selected_k': int(k),
        'parquet_file': str(out_parquet),
    }
    with open(out_parquet.with_suffix('.metadata.json'), 'w', encoding='utf8') as fh:
        json.dump(meta, fh, indent=2)
    print('Wrote metadata to', out_parquet.with_suffix('.metadata.json'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
