import sys
import csv

def load_gene_ids(csv_file):
    gene_ids = set()
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for idx,row in enumerate(reader):
            if row and idx > 0:  # Skip header
                gene_ids.add(row[0].strip())
    return gene_ids

def filter_fasta_by_genes(input_fasta, output_fasta, gene_ids):
    write_seq = False
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Extract gene id (remove '>')
                gene_id = line[1:].strip().split(" ")[0]
                if gene_id in gene_ids:
                    write_seq = True
                    outfile.write(line)
                else:
                    write_seq = False
            else:
                if write_seq:
                    outfile.write(line)

def main():
    if len(sys.argv) != 4:
        print("Usage: python filter_fasta_by_expression.py input_cleaned.fasta gene_ids.csv output_filtered.fasta")
        sys.exit(1)
    input_fasta = sys.argv[1]
    csv_file = sys.argv[2]
    output_fasta = sys.argv[3]
    gene_ids = load_gene_ids(csv_file)
    filter_fasta_by_genes(input_fasta, output_fasta, gene_ids)
    print(f"Filtered FASTA written to {output_fasta}")

if __name__ == "__main__":
    main()
