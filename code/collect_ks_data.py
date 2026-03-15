#!/usr/bin/env python3
"""
Collect KS test data from output/metrics folder and combine into a single table.

This script reads all KS test TSV files matching the pattern:
    *_ks_test.tsv
from the output/metrics folder and combines them into a single summary table.

Output columns:
    - Group: Cancer type | Gene | Mutation type (e.g., BRCA|BRCA1|snv)
    - Driver: Driver number (extracted from 'driver' column)
    - D: KS test statistic
    - 95% D_low: Lower bound of 95% CI for D
    - 95% D_high: Upper bound of 95% CI for D
    - p: P-value label
"""

import pandas as pd
import glob
import os
from pathlib import Path


def parse_driver_number(driver_str):
    """
    Convert driver string to number.

    Examples:
        'ratio' -> 1
        'two' -> 2
        'three' -> 3
        'four' -> 4
        etc.
    """
    driver_map = {
        "ratio": 1,
        "two": 2,
        "three": 3,
        "four": 4,
        "five": 5,
        "six": 6,
        "seven": 7,
        "eight": 8,
        "nine": 9,
        "ten": 10,
        "eleven": 11,
        "twelve": 12,
        "thirteen": 13,
        "fourteen": 14,
        "fifteen": 15,
        "sixteen": 16,
    }
    return driver_map.get(driver_str.lower(), driver_str)


def collect_ks_data(
    metrics_dir="output/metrics", output_file="output/combined_ks_test.tsv"
):
    """
    Collect all KS test data from TSV files and combine into a single table.

    Parameters:
        metrics_dir (str): Directory containing KS test TSV files
        output_file (str): Path to save the combined output table

    Returns:
        pd.DataFrame: Combined KS test data
    """
    # Find all KS test files
    pattern = os.path.join(metrics_dir, "*_ks_test.tsv")
    ks_files = glob.glob(pattern)

    if not ks_files:
        print(f"No KS test files found matching pattern: {pattern}")
        return None

    print(f"Found {len(ks_files)} KS test files")

    # List to store data from all files
    all_data = []

    # Process each file
    for filepath in sorted(ks_files):
        filename = os.path.basename(filepath)
        print(f"Processing: {filename}")

        try:
            # Read the TSV file
            df = pd.read_csv(filepath, sep="\t")

            # Check if required columns exist
            required_cols = ["driver", "D", "D_lo", "D_hi", "P_label", "group"]
            if not all(col in df.columns for col in required_cols):
                print(f"  Warning: Missing required columns in {filename}")
                continue

            # Extract relevant columns and rename
            df_subset = df[["driver", "D", "D_lo", "D_hi", "P_label", "group"]].copy()

            # Convert driver names to numbers
            df_subset["driver_num"] = df_subset["driver"].apply(parse_driver_number)

            # Rename columns to match output format
            df_subset = df_subset.rename(
                columns={
                    "group": "Group",
                    "driver_num": "Driver",
                    "D_lo": "95% D_low",
                    "D_hi": "95% D_high",
                    "P_label": "p",
                }
            )

            # Select and reorder columns
            df_subset = df_subset[
                ["Group", "Driver", "D", "95% D_low", "95% D_high", "p"]
            ]

            # Add to list
            all_data.append(df_subset)

        except Exception as e:
            print(f"  Error processing {filename}: {e}")
            continue

    if not all_data:
        print("No data collected from any files")
        return None

    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)

    # Sort by Group and Driver
    combined_df = combined_df.sort_values(["Group", "Driver"]).reset_index(drop=True)

    # Round numeric columns for cleaner output
    combined_df["D"] = combined_df["D"].round(4)
    combined_df["95% D_low"] = combined_df["95% D_low"].round(4)
    combined_df["95% D_high"] = combined_df["95% D_high"].round(4)

    # Save to file
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined_df.to_csv(output_file, sep="\t", index=False)
    print(f"\nCombined data saved to: {output_file}")
    print(f"Total rows: {len(combined_df)}")

    # Print summary statistics
    print(f"\nSummary:")
    print(f"  Unique groups: {combined_df['Group'].nunique()}")
    print(f"  Groups: {sorted(combined_df['Group'].unique())}")

    return combined_df


def main():
    """Main function to run the script."""
    print("=" * 80)
    print("Collecting KS Test Data")
    print("=" * 80)

    # Collect and combine data
    df = collect_ks_data()

    if df is not None:
        print("\n" + "=" * 80)
        print("First 10 rows of combined data:")
        print("=" * 80)
        print(df.head(10).to_string(index=False))

        print("\n" + "=" * 80)
        print("Last 10 rows of combined data:")
        print("=" * 80)
        print(df.tail(10).to_string(index=False))


if __name__ == "__main__":
    main()
