#!/data05/lindonghui/software/mamba/mambaforge/bin/python
import sys
from collections import defaultdict

# Check argument count (requires 3 arguments + script name)
if len(sys.argv) != 4:
    print("Usage: python script.py <input.blast> <output.txt> <bin_size>")
    sys.exit(1)

input_file = sys.argv[1]   # First argument: input file
output_file = sys.argv[2]  # Second argument: output file
bin_size = int(sys.argv[3])# Third argument: bin size (converted to integer)

# ------------------ Part 1: Determine maximum query length ------------------
query_length = 0

# First pass: read file to find maximum length
with open(input_file, 'r') as f:
    for line in f:
        # Skip comment lines and invalid lines
        if line.startswith('#') or len(line.strip().split()) < 9:
            continue

        # Extract columns 8-9 (Python uses 0-based indexing)
        fields = line.strip().split()
        s = int(fields[7])  # Column 8: alignment start
        e = int(fields[8])  # Column 9: alignment end

        # Update maximum length (consider both forward and reverse alignments)
        current_max = max(s, e)
        if current_max > query_length:
            query_length = current_max

# ------------------ Part 2: Initialize coverage array ------------------
# Calculate number of bins needed (rounded up)
num_bins = (query_length + bin_size - 1) // bin_size
coverage = [0] * num_bins  # Create array initialized with zeros

# ------------------ Part 3: Calculate coverage counts ------------------
# Second pass: process alignments
with open(input_file, 'r') as f:
    for line in f:
        # Filter invalid lines
        if line.startswith('#') or len(line.strip().split()) < 9:
            continue

        fields = line.strip().split()
        s = int(fields[7])
        e = int(fields[8])

        # Determine actual start/end positions (handles reverse alignments)
        start = min(s, e)
        end = max(s, e)

        # Calculate affected bin range
        start_bin = (start - 1) // bin_size  # e.g., position 100 belongs to bin 0 (1-100) when bin_size=100
        end_bin = (end - 1) // bin_size

        # Iterate through each affected bin
        for bin_idx in range(start_bin, end_bin + 1):
            if bin_idx >= num_bins:
                continue  # Prevent index out of bounds

            # Calculate exact bin range
            bin_start = bin_idx * bin_size + 1
            bin_end = (bin_idx + 1) * bin_size
            if bin_idx == num_bins - 1:  # Handle last bin specially
                bin_end = query_length

            # Calculate overlap region
            overlap_start = max(start, bin_start)
            overlap_end = min(end, bin_end)

            # Accumulate coverage length
            if overlap_start <= overlap_end:
                coverage[bin_idx] += (overlap_end - overlap_start + 1)

# ------------------ Part 4: Output results ------------------
with open(output_file, 'w') as out:
    for i in range(num_bins):
        bin_start = i * bin_size + 1
        bin_end = (i + 1) * bin_size
        if i == num_bins - 1:
            bin_end = query_length
        out.write(f"{bin_start}\t{bin_end}\t{coverage[i]}\n")