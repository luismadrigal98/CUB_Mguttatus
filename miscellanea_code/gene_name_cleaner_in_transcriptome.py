import sys

def clean_fasta_headers(input_fasta, output_fasta):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Keep only the first field (gene name)
                header = line.strip().split()[0]
                # Remove .1 only if it is at the end of the header
                if header.endswith('.1'):
                    header = header[:-2]
                outfile.write(header + '\n')
            else:
                outfile.write(line)

def main():
    if len(sys.argv) != 3:
        print("Usage: python gene_name_cleaner_in_transcriptome.py input.fasta output.fasta")
        sys.exit(1)
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    clean_fasta_headers(input_fasta, output_fasta)
    print(f"Cleaned FASTA written to {output_fasta}")

if __name__ == "__main__":
    main()
