#!/data05/lindonghui/software/mamba/mambaforge/bin/python
import sys

def parse_passlist(passlist_path):
    """Parse pass.list file and build chromosome region dictionary (auto-handles reverse sequences)"""
    ltr_regions = {}

    with open(passlist_path, 'r') as f:
        # Skip header line
        header = next(f)

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            # Parse LTR_loc field (first column)
            ltr_loc = fields[0]
            chr_part, coord = ltr_loc.split(':')
            start, end = map(int, coord.split('..'))

            # New: Auto-handle reverse sequences to ensure start < end
            if start > end:  # Handle reverse regions
                start, end = end, start  # Swap positions

            # Parse insertion time (last column)
            ins_time = fields[-1]

            # Store regions by chromosome
            if chr_part not in ltr_regions:
                ltr_regions[chr_part] = []
            ltr_regions[chr_part].append( (start, end, ins_time) )

    return ltr_regions  # Subsequent functions remain unchanged

# The following functions remain completely unchanged (same as original script)
def process_fasta(fasta_path, ltr_regions):
    """Process FASTA file to generate ID mapping table"""
    id_map = {}

    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                old_id = line[1:].strip()
                chr_part, pos = old_id.split(':')
                a_start, a_end = map(int, pos.split('-'))

                # Check if within LTR regions
                matched = False
                if chr_part in ltr_regions:
                    for (b_start, b_end, ins_time) in ltr_regions[chr_part]:
                        # Check region overlap
                        if (a_start <= b_end) and (a_end >= b_start):
                            new_id = f"{chr_part}_{a_start}_{a_end}_{ins_time}"
                            id_map[old_id] = new_id
                            matched = True
                            break  # Match first overlapping region

                # Keep original ID if no match
                if not matched:
                    id_map[old_id] = old_id
    return id_map

def rewrite_fasta(input_path, output_path, id_map):
    """Rewrite FASTA file"""
    with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                old_id = line[1:].strip()
                new_id = id_map[old_id]
                fout.write(f'>{new_id}\n')
            else:
                fout.write(line)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <input.fa> <pass.list> <output.fa>")
        sys.exit(1)

    input_fa = sys.argv[1]
    passlist_file = sys.argv[2]
    output_fa = sys.argv[3]

    # 1. Parse pass.list (with reverse handling added)
    ltr_regions = parse_passlist(passlist_file)

    # 2. Process FASTA to generate ID mapping (unchanged)
    id_mapping = process_fasta(input_fa, ltr_regions)

    # 3. Output new FASTA (unchanged)
    rewrite_fasta(input_fa, output_fa, id_mapping)
    print(f"Processing complete! Output file: {output_fa}")