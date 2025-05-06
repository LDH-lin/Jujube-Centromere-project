#!/data05/lindonghui/software/mamba/mambaforge/bin/python
import sys
import subprocess
import argparse

def runblast6(type, db, query, num, out, threads, evalue, mode, mt_mode):
    outfmt = "6 qseqid sseqid pident mismatch gapopen length qlen qstart qend slen sstart send"
    run_blast = [type,
                 "-db", db,
                 "-query", query,
                 "-outfmt", outfmt,
                 "-mt_mode", str(mt_mode),  # Use provided mt_mode value
                 "-out", f'{out}.blast',
                 "-num_threads", str(threads),
                 "-evalue", evalue]  # Set number of threads

    if num:
        run_blast.extend(["-max_target_seqs", str(num)])  # Add max_target_seqs if num provided

    run_blastall = [type,
                    "-db", db,
                    "-query", query,
                    "-out", f'{out}.blastall',
                    "-num_threads", str(threads),
                    "-mt_mode", str(mt_mode),  # Use provided mt_mode value
                    "-evalue", evalue]  # Include mt_mode 1

    if num:
        run_blastall.extend(["-max_target_seqs", str(num)])  # Add max_target_seqs if num provided

    try:
        if mode == "blast":
            subprocess.run(run_blast, check=True)
            print('Generated blast results')
        elif mode == "blastall":
            subprocess.run(run_blastall, check=True)
            print('Generated blastall results')
        elif mode == "both":
            subprocess.run(run_blast, check=True)
            subprocess.run(run_blastall, check=True)
            print('Generated both blast and blastall results')
        else:
            print(f'Unknown mode: {mode}')
            sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f'Execution error: {e}')
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Run BLAST with specified database, query, and output files.")
    parser.add_argument('-t', '--type', required=True)
    parser.add_argument('-d', '--db', required=True)
    parser.add_argument('-q', '--query', required=True)
    parser.add_argument('-n', '--num', required=False, type=int, help="Max target sequences (optional)")
    parser.add_argument('-o', '--out', required=True)
    parser.add_argument('-th', '--threads', required=False, default=1, type=int, help="Number of threads (default: 1)")
    parser.add_argument('-e', '--evalue', required=True)
    parser.add_argument('-m', '--mode', required=True, choices=["blast", "blastall", "both"],
                        help="Select run mode: 'blast', 'blastall' or 'both'")
    parser.add_argument('-mt', '--mt_mode', required=False, default=1, type=int,
                        help="Set mt_mode parameter (default: 1)")  # New mt_mode parameter
    args = parser.parse_args()

    runblast6(args.type, args.db, args.query, args.num, args.out, args.threads, args.evalue, args.mode, args.mt_mode)

if __name__ == "__main__":
    main()