#!/data05/lindonghui/software/mamba/mambaforge/bin/python
import sys

def read_bed_file(bed_file):
    """Read BED file and return a set of all regions"""
    regions = set()
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            regions.add((chrom, start, end))
    return regions

def parse_alignment_file(alignment_file):
    """Parse alignment file and return a list of alignment information"""
    alignments = []
    with open(alignment_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            seq_id = parts[0]
            chrom = parts[1]
            alignment_rate = float(parts[2])
            alignment_length = int(parts[5])
            start = int(parts[10])
            end = int(parts[11])
            alignments.append((seq_id, chrom, alignment_rate, alignment_length, start, end))
    return alignments

def calculate_abundance(alignment_rate, alignment_length):
    """Calculate alignment abundance"""
    return alignment_rate * alignment_length * 0.01

def main(alignment_file, bed_file, min_alignment_count):
    # Read BED file
    regions = read_bed_file(bed_file)

    # Parse alignment file
    alignments = parse_alignment_file(alignment_file)

    # Initialize statistics
    abundance_in_bed = {}
    abundance_out_of_bed = {}
    alignment_count = {}  # Track alignment counts per sequence

    # Process each alignment
    for alignment in alignments:
        seq_id, chrom, alignment_rate, alignment_length, start, end = alignment

        # Calculate abundance
        abundance = calculate_abundance(alignment_rate, alignment_length)

        # Update alignment count
        if seq_id in alignment_count:
            alignment_count[seq_id] += 1
        else:
            alignment_count[seq_id] = 1

        # Check region containment
        in_bed = False
        for region in regions:
            if chrom == region[0] and start >= region[1] and end <= region[2]:
                in_bed = True
                break

        # Update abundance statistics
        if in_bed:
            target_dict = abundance_in_bed
        else:
            target_dict = abundance_out_of_bed
            
        if seq_id in target_dict:
            target_dict[seq_id] += abundance
        else:
            target_dict[seq_id] = abundance

    # Output results with filtering
    for seq_id in set(abundance_in_bed.keys()).union(set(abundance_out_of_bed.keys())):
        if alignment_count[seq_id] < min_alignment_count:
            continue  # Skip sequences with insufficient alignments

        in_bed_abundance = abundance_in_bed.get(seq_id, 0)
        out_of_bed_abundance = abundance_out_of_bed.get(seq_id, 0)
        print(f"{seq_id}\t{in_bed_abundance}\t{out_of_bed_abundance}\t{alignment_count[seq_id]}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <alignment_file> <bed_file> <min_alignment_count>")
        print("Example: python script.py alignments.txt regions.bed 5")
        sys.exit(1)

    alignment_file = sys.argv[1]
    bed_file = sys.argv[2]
    try:
        min_alignment_count = int(sys.argv[3])
        if min_alignment_count < 0:
            raise ValueError
    except ValueError:
        print("Error: min_alignment_count must be a non-negative integer.")
        sys.exit(1)

    main(alignment_file, bed_file, min_alignment_count)