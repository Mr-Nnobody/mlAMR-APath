# Copy sampled 1000 k-mer dump files into an export folder for Colab/Kaggle
import shutil
from pathlib import Path

# change the path depending on which samples you are compiling
source_dir = Path('data/counted_kmers')
export_dir = Path('D:/MSc/Thesis/compressed_kmers')
export_dir.mkdir(parents=True, exist_ok=True)

if 'sampled_ids' not in globals() or not sampled_ids:
    sampled_ids_path = Path('data/sampled_data/cefuroxime_1444_ids.txt')
    if not sampled_ids_path.exists():
        raise FileNotFoundError(f'sampled_ids not found in memory and IDs file is missing: {sampled_ids_path}')
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

# # Optional: create a single zip file for easy upload
# zip_base = export_dir.parent
zip_base = 'D:/MSc/Thesis/cefuroxime_1364'
zip_path = shutil.make_archive(str(zip_base), 'zip', root_dir=export_dir)
print(f'Created zip: {zip_path}')