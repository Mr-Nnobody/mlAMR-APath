import os
import glob
import shutil
import subprocess


def count_kmers(fasta_file, output_prefix, kmer_size=10):
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Paths to KMC binaries (bundled under scripts/kmer/bin)
    kmc_bin_dir = os.path.join(script_dir, "kmc", "bin")
    kmc_exe = os.path.join(kmc_bin_dir, "kmc.exe")
    kmc_tools_exe = os.path.join(kmc_bin_dir, "kmc_tools.exe")

    # KMC requires a temporary directory for RAM overflow (relative to project root / CWD)
    tmp_dir = "kmc_tmp"
    os.makedirs(tmp_dir, exist_ok=True)

    # STEP 1: Count K-mers   
    kmc_cmd = [
        kmc_exe,
        f"-k{kmer_size}",
        "-fm",      # -fm: FASTA format
        "-ci1",     # -ci1: minimum count to report (1 means show all)
        "-m4",   # use at most 4 GB RAM; try 2 if needed: "-m2"
        "-t2",   # (optional) 2 threads to reduce peak memory
        fasta_file,
        output_prefix,
        tmp_dir,
    ]

    print("Running KMC Count:", " ".join(kmc_cmd))
    subprocess.run(kmc_cmd, check=True)

    # STEP 2: Dump to Text
    dump_file = f"{output_prefix}_kmers.txt"
    dump_cmd = [
        kmc_tools_exe,
        "transform",
        output_prefix,
        "dump",
        dump_file,
    ]

    print("Running KMC Dump:", " ".join(dump_cmd))
    subprocess.run(dump_cmd, check=True)

    # STEP 3: Cleanup KMC database and temporary files
    for kmc_db_file in glob.glob(f"{output_prefix}.kmc_*"):
        try:
            os.remove(kmc_db_file)
        except OSError:
            pass

    if os.path.isdir(tmp_dir):
        try:
            shutil.rmtree(tmp_dir)
        except OSError:
            pass

    print(f"Finished! File saved to {dump_file} (temporary KMC files removed)")


with open ('data/genome_ids.txt', 'r') as f:
    genome_ids = [line.strip() for line in f]
rows = []
for gid in genome_ids:
    try:
       count_kmers(f"data/fasta_files/{gid}.fasta", f"output/{gid}_db")
       
        
    except Exception as e:
        with open('data/kmer_not_counted.txt', 'a') as p:
           p.write(f"{gid}\n")
        continue