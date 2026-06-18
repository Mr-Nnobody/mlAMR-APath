# Copy sampled 1000 k-mer dump files into an export folder for Colab/Kaggle
import shutil
from pathlib import Path

source_dir = Path('../data/counted_kmers')
export_dir = Path('../data/export_kaggle_colab/counted_kmers_1000')
export_dir.mkdir(parents=True, exist_ok=True)

if 'sampled_ids' not in globals() or not sampled_ids:
    sampled_ids_path = Path('../data/phenotype/ampicillin_1000_ids.txt')
    if not sampled_ids_path.exists():
        raise FileNotFoundError('sampled_ids not found in memory and IDs file is missing: ../data/phenotype/ampicillin_1000_ids.txt')
    sampled_ids = [line.strip() for line in sampled_ids_path.read_text(encoding='utf8').splitlines() if line.strip()]

copied = 0
missing = []
for gid in sampled_ids:
    src = source_dir / f'{gid}_db_kmers.txt'
    dst = export_dir / src.name
    if src.exists():
        shutil.copy2(src, dst)
        copied += 1
    else:
        missing.append(str(gid))

print(f'Copied {copied} files to {export_dir}')
if missing:
    print(f'Missing files: {len(missing)}')
    print('First few missing GenomeIDs:', missing[:10])

# Optional: create a single zip file for easy upload
zip_base = export_dir.parent / 'counted_kmers_1000'
zip_path = shutil.make_archive(str(zip_base), 'zip', root_dir=export_dir)
print(f'Created zip: {zip_path}')