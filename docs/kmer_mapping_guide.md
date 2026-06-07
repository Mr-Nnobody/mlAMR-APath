# K-mer Mapping Note

This project uses sparse k-mer features for AMR prediction. The matrix itself is only part of the story: to interpret model outputs biologically, each feature column must remain traceable back to its original k-mer sequence.

## Purpose

The goal of this note is to describe how to preserve that traceability, how to map important model features back to genomic regions, and how to strengthen the overall project so it is both technically rigorous and biologically useful.

## Core idea

A sparse wide matrix is appropriate for modeling, but it does not need to carry human-readable column names in memory. Instead, the feature identity should be stored separately in a column index file.

You do not need to save the full sparse matrix to do the interpretation phase later. What you need is:

- the selected k-mer vocabulary in a fixed order,
- the trained model,
- and the original genome FASTA or per-genome k-mer dump used at prediction time.

At inference time, rebuild the sample's sparse vector from the FASTA/k-mer dump, score it with the model, and then map only the important active k-mers back to genomic coordinates or known resistance markers.

Recommended project outputs:

- `kmer_matrix_sparse.npz` for the sparse feature matrix
- `kmer_matrix_genome_index.csv` for row-to-GenomeID mapping
- `kmer_matrix_kmer_index.csv` for column-to-k-mer mapping

Together, these files preserve the exact relationship between model features and biological sequences.

## Workflow overview

1. Generate the long-format table with rows of `GenomeID, kmer, count`.
2. Filter or select a manageable k-mer vocabulary.
3. Save the vocabulary order as the canonical feature index.
4. Train a model on the sparse representation.
5. At prediction time, rebuild the sparse vector for the new genome from its FASTA or k-mer dump.
6. Score the genome and identify the active k-mers that contributed most to the prediction.
7. Translate those indices back into k-mer strings using the saved index file.
8. Search the k-mers against the genome FASTA file to identify likely genomic regions.
9. Group overlapping or nearby signals and check whether they align with known genes, operons, or mobile elements.

## Why the index file is essential

Sparse matrices are efficient because they store values by numeric position rather than by large dense labels. That is ideal for memory use, but it means the biological meaning of each column must be preserved externally. Without the k-mer index file, the model can still run, but the resulting feature importance values cannot be interpreted reliably.

The index file solves that problem by making the mapping explicit:

- matrix column `0` corresponds to a specific k-mer string
- matrix column `1` corresponds to another k-mer string
- and so on

If you later want to map a prediction back to biology, the index file is the bridge from model feature number to actual sequence. The sparse matrix itself is only the temporary numeric container used for training or scoring.

## Do you need to save the sparse matrix?

Usually no.

Save the matrix only if you want a cache of the exact training dataset for reruns or debugging. For the prediction and mapping phase, the useful artifacts are:

- `vocab_selected_*.pkl` or a CSV index file with ordered k-mers
- the fitted model file
- optionally a scaler or transformer if one was used
- the genome FASTA file or per-genome k-mer dump for the sample being analyzed

That is enough to:

1. rebuild the sample's feature vector,
2. score the model,
3. inspect the non-zero k-mers,
4. map those k-mers back to the genome sequence,
5. and annotate nearby genes or resistance markers.

## How to map features back to genomic regions

### 1. Identify informative features

After training, rank features using model coefficients, permutation importance, or SHAP values when available.

### 2. Recover the k-mer sequence

Use the saved k-mer index file to convert each important feature index into its actual k-mer string.

### 2b. Rebuild the sample feature vector

For a new genome, create the same sparse vector using the saved vocabulary order. This gives you the exact active k-mers for the prediction without needing the original training matrix.

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

## Prediction and mapping workflow

Use this flow after you have trained a model and want to explain a new genome prediction.

### Save these artifacts from training

- the selected k-mer vocabulary in order, for example `vocab_selected_binary.pkl` or `vocab_selected_counts.pkl`
- the fitted model, for example `logistic_binary.joblib` or `rf_counts.joblib`
- the scaler or transformer, if the model used one
- the antibiotic-specific label file or metadata used for training
- optionally the training genome IDs, if you want reproducible splits and traceability

### Rebuild the prediction vector

For a new genome:

1. read the FASTA file or k-mer dump,
2. count the k-mers that appear in the saved vocabulary,
3. build a sparse vector using the exact saved vocabulary order,
4. apply the same preprocessing steps used during training,
5. run the model to get the resistance prediction.

### Map prediction back to biology

Once you have the prediction, take the non-zero or highly weighted k-mers and:

1. translate the feature indices back into k-mer strings using the saved vocabulary,
2. locate each k-mer in the genome sequence,
3. cluster nearby hits into candidate genomic regions,
4. annotate those regions against known genes, operons, plasmids, transposons, or resistance islands,
5. summarize the strongest matches as likely resistance markers.

### Practical rule

If you did not save the training sparse matrix, that is fine. You can still do the full mapping phase as long as you saved the feature vocabulary and can rebuild the sample vector from the genome FASTA or k-mer dump.
