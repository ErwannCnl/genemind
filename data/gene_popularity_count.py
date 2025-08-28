#!/usr/bin/env python3
import csv
import sys

def get_gene_value(file_path: str, gene_name: str):
    """
    Search a TSV file for a gene name and return the value in the last column
    of the first matching row. If not found, return None.
    """
    try:
        with open(file_path, newline='') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                # match only the seventh column (gene symbol)
                if row and row[6] == gene_name:
                    return row[-1]
    except FileNotFoundError:
        print(f"❌ File '{file_path}' not found.", file=sys.stderr)
        sys.exit(1)

    return None


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <GENE_NAME>", file=sys.stderr)
        sys.exit(1)

    file = "all_gene_counts.tsv"
    gene = sys.argv[1]
    value = get_gene_value(file, gene)

    print(value if value is not None else "NA")

