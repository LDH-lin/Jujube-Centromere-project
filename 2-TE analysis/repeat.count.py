#!/usr/bin/env python
import sys

# Define merging rules for element types
def merge_element_type(element_type):
    if element_type.startswith("DNA"):
        return "DNA"
    elif element_type == "LTR/Copia":
        return "LTR/Copia"
    elif element_type == "LTR/Gypsy":
        return "LTR/Gypsy"
    elif element_type in ["Satellite", "Low_complexity"]:
        return "Satellite"
    elif element_type == "Unknown":
        return "Unknown"
    else:
        return "Other"

# Parse TRF file and merge overlapping intervals
def parse_trf_file(trf_file):
    intervals = {}
    with open(trf_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            parts = line.split('\t')
            if len(parts) < 6:  # Ensure sufficient fields
                continue  # Skip malformed lines
            # Extract period information
            info_field = parts[8]
            info_pairs = info_field.split(';')
            period = None
            for pair in info_pairs:
                if pair.startswith('period='):
                    period_str = pair.split('=')[1]
                    try:
                        period = int(period_str)
                    except ValueError:
                        period = None
                        break  # Skip if non-integer period
                    break  # Stop after finding period
            if period is None or period < 10:
                continue  # Skip lines with missing/small period
            chromosome = parts[0]
            try:
                start = int(parts[3])
                end = int(parts[4])
            except (IndexError, ValueError):
                continue  # Skip invalid positions
            if chromosome not in intervals:
                intervals[chromosome] = []
            intervals[chromosome].append((start, end))
    return intervals

# Merge overlapping intervals algorithm
def merge_intervals(intervals_list):
    if not intervals_list:
        return []
    # Sort intervals by start position
    sorted_intervals = sorted(intervals_list, key=lambda x: x[0])
    merged = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            # Merge overlapping/adjacent intervals
            new_start = last[0]
            new_end = max(last[1], current[1])
            merged[-1] = (new_start, new_end)
        else:
            merged.append(current)
    return merged

# Calculate total length of merged intervals
def calculate_total_length(merged_intervals):
    total = 0
    for start, end in merged_intervals:
        total += (end - start + 1)
    return total

# Process TRF GFF file and calculate total TRF length
def trf_gff_(trf_file):
    intervals = parse_trf_file(trf_file)
    trf_total_length = 0
    trf_intervals = {}  # Store merged intervals per chromosome
    for chrom, chrom_intervals in intervals.items():
        merged = merge_intervals(chrom_intervals)
        chrom_length = calculate_total_length(merged)
        trf_total_length += chrom_length
        trf_intervals[chrom] = merged
    return trf_total_length, trf_intervals

# Parse BED file to get chromosomal regions
def parse_bed_file(bed_file):
    chromosome_info = {}  # Track min_start and max_end per chromosome
    with open(bed_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            parts = line.split('\t')
            if len(parts) < 3:
                continue  # Skip malformed lines
            chromosome = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue  # Skip invalid positions
            if chromosome not in chromosome_info:
                chromosome_info[chromosome] = {
                    'min_start': start,
                    'max_end': end
                }
            else:
                if start < chromosome_info[chromosome]['min_start']:
                    chromosome_info[chromosome]['min_start'] = start
                if end > chromosome_info[chromosome]['max_end']:
                    chromosome_info[chromosome]['max_end'] = end
    return chromosome_info

# Process RepeatMasker .out file while excluding TRF overlaps
def parse_repeatmasker_out(file_path, trf_intervals):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith("\s*SW") and not line.startswith("\s*score"):  # Skip headers
                parts = line.strip().split()
                chromosome = parts[4]  # Chromosome name
                query_start = int(parts[5])
                query_end = int(parts[6])
                element_type = parts[10]
                length = query_end - query_start + 1
                merged_type = merge_element_type(element_type)

                # Check TRF overlap
                if chromosome in trf_intervals:
                    overlap = False
                    for (start, end) in trf_intervals[chromosome]:
                        if not (query_end < start or query_start > end):
                            overlap = True
                            break
                    if overlap:
                        continue  # Skip overlapping records
                data.append((merged_type, length))
    return data

# Calculate total genome length from BED regions
def calculate_genome_length(chromosome_info):
    genome_length = 0
    for chrom, info in chromosome_info.items():
        chrom_length = info['max_end'] - info['min_start'] + 1
        genome_length += chrom_length
    return genome_length

# Calculate lengths and percentages for each element type
def calculate_length_and_percentage(data, trf_total_length, genome_length):
    # Calculate total repeatmasker length
    repeatmasker_total_length = sum(length for _, length in data)
    element_total_length = repeatmasker_total_length + trf_total_length

    # Calculate non-element regions
    non_element_length = genome_length - element_total_length

    # Initialize length dictionary
    length_dict = {"TRF": trf_total_length}
    for element_type, length in data:
        length_dict[element_type] = length_dict.get(element_type, 0) + length
    length_dict["Non-element"] = non_element_length

    # Calculate percentages
    result_dict = {}
    for element_type, length in length_dict.items():
        percentage = (length / genome_length) * 100
        result_dict[element_type] = {"Length": length, "Percentage": percentage}
    return result_dict

# Main workflow
def main(repeatmasker_file, trf_file, bed_file, output_file):
    # Process TRF data
    trf_total_length, trf_intervals = trf_gff_(trf_file)

    # Get chromosomal regions
    chromosome_info = parse_bed_file(bed_file)

    # Process RepeatMasker data
    data = parse_repeatmasker_out(repeatmasker_file, trf_intervals)

    # Calculate genomic metrics
    genome_length = calculate_genome_length(chromosome_info)
    result_dict = calculate_length_and_percentage(data, trf_total_length, genome_length)

    # Save results
    with open(output_file, 'w') as file:
        file.write("ElementType\tLength\tPercentage\n")
        for element_type, values in result_dict.items():
            file.write(f"{element_type}\t{values['Length']}\t{values['Percentage']:.2f}\n")
        # Write summary stats
        all_elements_length = sum(values['Length'] for values in result_dict.values() if element_type != "Non-element")
        file.write(f"\nAll Elements Total Length:\t{all_elements_length}\n")
        file.write(f"RepeatMasker Total Length:\t{sum(length for _, length in data)}\n")
    print(f"Results saved to {output_file}")

# File paths from command-line arguments
repeatmasker_file = sys.argv[1]  # RepeatMasker .out file
trf_file = sys.argv[2]          # TRF file
bed_file = sys.argv[3]           # BED file
output_file = sys.argv[4]        # Output file

# Execute main program
main(repeatmasker_file, trf_file, bed_file, output_file)