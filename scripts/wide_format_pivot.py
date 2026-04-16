import pandas as pd
kmers_long = pd.read_csv("kmer-data-long.csv")  # GenomeID, kmer, count

kmer_matrix = kmers_long.pivot_table(
    index="GenomeID",
    columns="kmer",
    values="count",
    aggfunc="sum",
    fill_value=0,
).reset_index()  # <- makes GenomeID a column
print('Pivoted to wide format, now saving to CSV...')

print(kmer_matrix.head())
kmer_matrix.to_csv("kmer-data-wide.csv", index=False)
# Optionally make it sparse for memory
# kmer_matrix.iloc[:, 1:] = kmer_matrix.iloc[:, 1:].astype("Sparse[int]")