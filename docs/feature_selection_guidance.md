# Feature Selection Guidance for E. coli 10-mer AMR (Binary)

This note outlines a scientifically motivated, reproducible feature-selection pipeline for
E. coli 10-mer k-mers with binary AMR labels (R/I/S), based on ~7,800 genomes. Each step
includes a rationale that is biological or statistical, not just technical.

## 1) Decide how to treat "Intermediate" (I)

**Option A: Treat I as Resistant**

- Rationale: Intermediate often reflects reduced susceptibility and can indicate partial
  resistance mechanisms. Grouping I with R can capture clinically conservative signals.

**Option B: Exclude I from training**

- Rationale: Intermediate can be a heterogeneous category and may blur the phenotype signal.
  Removing it can yield a cleaner separation between R and S for feature discovery.

Pick one and document it. For model comparisons, test both strategies if time allows.

## 2) Prevalence filter (remove rare k-mers)

**Recommendation**: keep k-mers present in at least 1% of genomes

- With 7,800 genomes, 1% corresponds to about 78 genomes.
- Rationale: extremely rare k-mers are more likely to be sequencing/assembly artifacts orprivate variants that do not generalize. Removing them reduces noise and improves reproducibility across cohorts.

min_frac: minimum fraction of genomes that must contain a k‑mer to keep it.
If a k‑mer is too rare (only in a handful of genomes), it’s often noise or too specific to be useful.

Optional: evaluate 0.5% and 2% as sensitivity analyses.

## 3) Upper prevalence filter (remove nearly fixed k-mers)

**Recommendation**: remove k-mers present in more than 98-99% of genomes

- With 7,800 genomes, 98% corresponds to about 7,644 genomes.
- Rationale: near-constant k-mers represent core-genome background and add little to
  discrimination between resistant and susceptible phenotypes.

max_frac: maximum fraction of genomes that may contain a k‑mer to keep it.
If a k‑mer is in almost every genome, it carries little discriminative signal (it’s nearly constant)'

## 4) Association-based selection (label-linked signal)

**Recommendation**: rank remaining k-mers by association with the label and keep top N

- Use chi-square or mutual information (MI) with the binary label.
- Rationale: AMR is often mediated by mobile elements and specific genes. Association tests enrich for k-mers that track with these biological factors.

**Starting N**: 50,000 to 200,000 k-mers

- The exact number should be tuned based on downstream model performance and stability.

## 5) Stability selection (reproducible signals)

**Recommendation**: bootstrap the data and keep k-mers selected repeatedly by an L1 model

- Example: run L1-logistic regression on 50-100 bootstraps, keep features that appear in
  at least 50% of runs.
- Rationale: features that repeat across resamples are more likely to be biologically
  stable and less likely to be dataset-specific noise.

## 6) Suggested practical pipeline

1. Choose I handling (merge with R or exclude)
2. Build a k-mer presence matrix (binary presence is usually enough for selection)
3. Apply prevalence filter (>= 1% and <= 98-99%)
4. Rank by chi-square or MI and keep top 50k-200k
5. Optional stability selection for a high-confidence set
6. Train models and validate on held-out data

## 7) Why this is scientifically defensible

- It balances sensitivity (detecting true AMR signals) with specificity (avoiding noise).
- It reduces the multiple-testing burden inherent to millions of possible k-mers.
- It aligns with known AMR biology: resistance is driven by specific genomic elements,
  not random rare variants or universal core-genome k-mers.

## 8) Reporting and reproducibility

Always record:

- The exact prevalence thresholds
- The association metric used (chi-square or MI)
- The number of retained k-mers
- How "Intermediate" was handled

This makes results interpretable and comparable across studies.

## 9) Python implementation (end-to-end example)

Below is a concrete, reproducible implementation that:

- handles the I label
- computes k-mer prevalence directly from per-genome dumps
- applies prevalence filters
- selects top k-mers by chi-square
- builds a sparse matrix for modeling

This assumes each genome has a dump file like:
`output/counted_kmers/{GenomeID}_db_kmers.txt` with lines: `kmer count`.

```python
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.feature_selection import chi2


def load_labels(labels_path: str, treat_intermediate_as_resistant: bool) -> pd.Series:
  df = pd.read_csv(labels_path)
  # Expect columns: GenomeID, phenotype (R/I/S)
  pheno = df.set_index("GenomeID")["phenotype"].str.upper()

  if treat_intermediate_as_resistant:
    pheno = pheno.replace({"I": "R"})
    pheno = pheno[pheno.isin(["R", "S"])]
  else:
    pheno = pheno[pheno.isin(["R", "S"])]

  y = pheno.map({"R": 1, "S": 0})
  return y


def iter_genome_kmers(dump_path: Path) -> dict[str, int]:
  kmers: dict[str, int] = {}
  with dump_path.open("r") as f:
    for line in f:
      parts = line.strip().split()
      if len(parts) != 2:
        continue
      kmer, cnt = parts
      kmers[kmer] = int(cnt)
  return kmers


def build_prevalence(dump_dir: str, genome_ids: list[str]) -> Counter:
  prevalence = Counter()
  for gid in genome_ids:
    dump_path = Path(dump_dir) / f"{gid}_db_kmers.txt"
    if not dump_path.exists():
      continue
    kmers = iter_genome_kmers(dump_path)
    prevalence.update(kmers.keys())
  return prevalence


def select_vocab_by_prevalence(prevalence: Counter, n_genomes: int,
                 min_frac: float = 0.01, max_frac: float = 0.99) -> list[str]:
  min_count = int(np.ceil(min_frac * n_genomes))
  max_count = int(np.floor(max_frac * n_genomes))
  vocab = [k for k, c in prevalence.items() if min_count <= c <= max_count]
  return vocab


def build_sparse_matrix(dump_dir: str, genome_ids: list[str], vocab: list[str]) -> sparse.csr_matrix:
  vocab_index = {kmer: idx for idx, kmer in enumerate(vocab)}
  rows = []
  cols = []
  data = []

  for row_idx, gid in enumerate(genome_ids):
    dump_path = Path(dump_dir) / f"{gid}_db_kmers.txt"
    if not dump_path.exists():
      continue
    kmers = iter_genome_kmers(dump_path)
    for kmer, cnt in kmers.items():
      col_idx = vocab_index.get(kmer)
      if col_idx is None:
        continue
      rows.append(row_idx)
      cols.append(col_idx)
      data.append(cnt)

  mat = sparse.csr_matrix((data, (rows, cols)),
              shape=(len(genome_ids), len(vocab)),
              dtype=np.int32)
  return mat


def select_top_kmers_chi2(X: sparse.csr_matrix, y: np.ndarray, k: int) -> np.ndarray:
  scores, _ = chi2(X, y)
  top_idx = np.argsort(scores)[::-1][:k]
  return top_idx


def main():
  dump_dir = "output/counted_kmers"
  labels_path = "data/isolates.csv"  # adjust if needed
  genome_ids_path = "data/genome_ids.txt"

  with open(genome_ids_path, "r") as f:
    genome_ids = [line.strip() for line in f if line.strip()]

  y = load_labels(labels_path, treat_intermediate_as_resistant=True)
  genome_ids = [gid for gid in genome_ids if gid in y.index]

  prevalence = build_prevalence(dump_dir, genome_ids)
  vocab = select_vocab_by_prevalence(prevalence, n_genomes=len(genome_ids),
                     min_frac=0.01, max_frac=0.99)

  X = build_sparse_matrix(dump_dir, genome_ids, vocab)
  y_arr = y.loc[genome_ids].to_numpy()

  top_idx = select_top_kmers_chi2(X, y_arr, k=100000)
  X_sel = X[:, top_idx]
  vocab_sel = [vocab[i] for i in top_idx]

  print(f"Matrix shape after selection: {X_sel.shape}")
  print(f"Selected vocab size: {len(vocab_sel)}")


if __name__ == "__main__":
  main()
```

Notes:

- This uses count values; for binary presence/absence, replace `data.append(cnt)` with `1`.
- If memory becomes a bottleneck during prevalence counting, consider writing the prevalence Counter to disk periodically (pickle) and merging in batches.
