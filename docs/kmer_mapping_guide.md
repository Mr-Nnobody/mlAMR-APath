# K-mer Mapping Note

This project uses sparse k-mer features for AMR prediction. The matrix itself is only part of the story: to interpret model outputs biologically, each feature column must remain traceable back to its original k-mer sequence.

## Purpose

The goal of this note is to describe how to preserve that traceability, how to map important model features back to genomic regions, and how to strengthen the overall project so it is both technically rigorous and biologically useful.

## Core idea

A sparse wide matrix is appropriate for modeling, but it does not need to carry human-readable column names in memory. Instead, the feature identity should be stored separately in a column index file.

Recommended project outputs:

- `kmer_matrix_sparse.npz` for the sparse feature matrix
- `kmer_matrix_genome_index.csv` for row-to-GenomeID mapping
- `kmer_matrix_kmer_index.csv` for column-to-k-mer mapping

Together, these files preserve the exact relationship between model features and biological sequences.

## Workflow overview

1. Generate the long-format table with rows of `GenomeID, kmer, count`.
2. Filter or select a manageable k-mer vocabulary.
3. Convert the long table into a sparse wide matrix for modeling.
4. Save the row and column index files alongside the matrix.
5. Train a model on the sparse representation.
6. Extract the most important feature indices from the fitted model.
7. Translate those indices back into k-mer strings using the saved index file.
8. Search the k-mers against genome FASTA files to identify likely genomic regions.
9. Group overlapping or nearby signals and check whether they align with known genes, operons, or mobile elements.

## Why the index file is essential

Sparse matrices are efficient because they store values by numeric position rather than by large dense labels. That is ideal for memory use, but it means the biological meaning of each column must be preserved externally. Without the k-mer index file, the model can still run, but the resulting feature importance values cannot be interpreted reliably.

The index file solves that problem by making the mapping explicit:

- matrix column `0` corresponds to a specific k-mer string
- matrix column `1` corresponds to another k-mer string
- and so on

## How to map features back to genomic regions

### 1. Identify informative features

After training, rank features using model coefficients, permutation importance, or SHAP values when available.

### 2. Recover the k-mer sequence

Use the saved k-mer index file to convert each important feature index into its actual k-mer string.

### 3. Locate the k-mer in the genome

Search the k-mer sequence against the genome FASTA file for the relevant GenomeID. If the k-mer appears in multiple places, keep all candidate locations.

### 4. Consolidate adjacent signals

If multiple significant k-mers cluster within the same gene or genomic neighborhood, treat that cluster as stronger evidence than an isolated feature.

### 5. Annotate the candidate region

Compare the region against known genes, resistance islands, plasmids, transposons, or other mobile elements to support a biological interpretation.

## Recommendations to strengthen the project

- Use a sparse-aware classifier such as logistic regression with `saga`, linear SVM, or SGD-based learning.
- Train one model per antibiotic so the labels stay clear and the results remain interpretable.
- Filter k-mers by minimum genome frequency or top-k selection before modeling to control dimensionality.
- Keep a fixed train, validation, and test split so results are reproducible.
- Report AUROC, AUPRC, F1, precision, and recall rather than relying on a single metric.
- Include a simple baseline model to show that the k-mer approach adds value.
- Report how many genomes and k-mers remain after each filtering step.
- Highlight memory savings from the sparse pipeline compared with a dense pivot.
- Map the strongest k-mers back to genes and discuss biological relevance, not only classification performance.
- Add at least one interpretation figure that links significant k-mers to candidate resistance regions.

## Strong final deliverables

The project will be stronger if it ends with a small but complete package of outputs:

- a reproducible sparse feature pipeline,
- a per-antibiotic prediction model,
- a table of top features and their mapped k-mers,
- a genomic region mapping figure,
- and a short discussion of known resistance mechanisms that match the top signals.

That combination makes the work look technically solid, biologically grounded, and easy to defend.
