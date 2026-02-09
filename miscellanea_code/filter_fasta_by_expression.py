import sys
import csv

def load_gene_ids(txt_file, sep='\t'):
    gene_ids = set()
    with open(txt_file, 'r') as f:
        reader = csv.reader(f, delimiter=sep)
        next(reader, None)  # Skip header row
        for row in reader:
            if row:
                gene_ids.add(row[0].strip())
    return gene_ids

def filter_fasta_by_genes(input_fasta, output_fasta, gene_ids):
    write_seq = False
    found_genes = set()
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                gene_id = line[1:].strip().split()[0]
                if gene_id in gene_ids:
                    write_seq = True
                    outfile.write(line)
                    found_genes.add(gene_id)
                else:
                    write_seq = False
            else:
                if write_seq:
                    outfile.write(line)
    return found_genes

def main():
    if len(sys.argv) != 5:
        print("Usage: python filter_fasta_by_expression.py input_cleaned.fasta gene_ids.txt output_filtered.fasta")
        sys.exit(1)
    input_fasta = sys.argv[1]
    txt_file = sys.argv[2]
    sep = sys.argv[3]
    output_fasta = sys.argv[4]
    gene_ids = load_gene_ids(txt_file, sep)
    found_genes = filter_fasta_by_genes(input_fasta, output_fasta, gene_ids)
    print(f"Filtered FASTA written to {output_fasta}")
    print(f"Found {len(found_genes)} genes in the FASTA file that matched the gene IDs from the text file.")

if __name__ == "__main__":
    main()
