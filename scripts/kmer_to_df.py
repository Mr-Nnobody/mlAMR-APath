import pandas as pd

def kmc_dump_to_series(genome_id: str, dump_path: str) -> pd.Series:
    data = {}
    with open(dump_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            kmer, cnt = parts
            data[kmer] = int(cnt)
    return pd.Series(data, name=genome_id)

with open ('data/genome_ids.txt', 'r') as f:
    genome_ids = [line.strip() for line in f]
rows = []
for gid in genome_ids:
    try: 
        dump_file = f"output/{gid}_db_kmers.txt"
        rows.append(kmc_dump_to_series(gid, dump_file))
    except:
        with open ('data/not_counted_ids.txt', 'a') as p:
            p.write(f'{gid}\n')
        continue

kmer_matrix = pd.DataFrame(rows).fillna(0)
kmer_matrix.to_csv('kmer-data.csv')
print('Done')
# kmer_matrix.index.name = "GenomeID"