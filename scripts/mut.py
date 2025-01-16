import pandas as pd
import re

# Function to load and validate the input files
def load_data(file_path, sep, col_names):
    try:
        # Load data with appropriate column names
        data = pd.read_csv(file_path, sep=sep, header=None, names=col_names, on_bad_lines="skip")
        print(f"Successfully loaded {file_path}.")
        return data
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        exit(1)
    except pd.errors.ParserError as e:
        print(f"Error parsing {file_path}: {e}")
        exit(1)

# Function to extract counts from the `{ALLELE:COUNT}` format
def extract_mutation_count(allele_count_str):
    print(f"Processing: {allele_count_str}")  # Debugging print
    # This regular expression matches "A:1" or "GTACCTCCT:2" or similar patterns
    counts = re.findall(r'(\S+):(\d+)', allele_count_str.replace(" ", ""))
    print(f"Extracted Counts: {counts}")  # Debugging print
    total_count = sum(int(count) for allele, count in counts)  # Sum all mutation counts
    return total_count

# Main function to calculate median background mutation level
def main():
    # File paths
    depth_file = 'depth_info.INFO'
    variant_file = 'normal_counts.frq.count'

    # Load depth information and variant counts
    depth = load_data(depth_file, sep='\t', col_names=['CHROM', 'POS', 'DP'])
    variant_counts = load_data(variant_file, sep='\t', col_names=['CHROM', 'POS', '{ALLELE:COUNT}'])

    # Ensure columns are numeric
    try:
        depth['DP'] = pd.to_numeric(depth['DP'], errors='coerce')
    except ValueError as e:
        print(f"Error converting columns to numeric: {e}")
        exit(1)

    # Extract mutation counts from '{ALLELE:COUNT}' column
    variant_counts['COUNT'] = variant_counts['{ALLELE:COUNT}'].apply(extract_mutation_count)

    # Debugging: Print the first few entries of the COUNT column
    print("\nVariant Counts with Mutation Count:")
    print(variant_counts[['CHROM', 'POS', '{ALLELE:COUNT}', 'COUNT']].head(20))  # Checking first 20 entries

    # Drop rows with invalid numeric data
    depth = depth.dropna(subset=['DP'])
    variant_counts = variant_counts.dropna(subset=['COUNT'])

    # Check for missing mutation counts
    missing_counts = variant_counts['COUNT'].isna().sum()
    if missing_counts > 0:
        print(f"Warning: {missing_counts} rows have missing mutation counts (NaN) in the 'COUNT' column.")

    # Check total mutations and depth for debugging
    total_mutations = variant_counts['COUNT'].sum()
    total_depth = depth['DP'].sum()

    # Debugging prints
    print(f"Total Mutations (COUNT sum): {total_mutations}")
    print(f"Total Depth (DP sum): {total_depth}")
    
    # Show basic statistics of the data
    print("\nDepth Info Sample:")
    print(depth.describe())
    print("\nVariant Counts Sample:")
    print(variant_counts.describe())

    # Handle division by zero
    if total_depth == 0:
        print("Error: Total depth is zero, cannot calculate background mutation level.")
        exit(1)

    # Median background mutation level
    background_mutation_level = (total_mutations / total_depth) * 1e6
    print(f"Median Background Mutation Level: {background_mutation_level:.2f} mutations per million reads")

# Entry point
if __name__ == "__main__":
    main()

