import sys
import gzip
import csv
import re

def get_complement(nuc):
    """Return the complement of a single nucleotide."""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return mapping.get(nuc, 'N')

def parse_gff_introns(gff_file):
    """
    Parses GFF3 to find gene coordinates/strand and subtracts exons to get introns.
    Returns: {chrom: [(start, end, strand), ...]}
    """
    print(f"Loading GFF: {gff_file}")
    genes = {}
    exons = {}
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            chrom = parts[0]
            feature = parts[2]
            start = int(parts[3]) - 1 # 0-based
            end = int(parts[4])
            strand = parts[6] # '+' or '-'
            
            if feature == 'gene':
                # Parse ID
                gid_match = re.search(r'ID=([^;]+)', parts[8])
                if gid_match:
                    gid = gid_match.group(1)
                    genes[gid] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand}
            elif feature == 'exon':
                # Parse Parent
                pid_match = re.search(r'Parent=([^;]+)', parts[8])
                if pid_match:
                    parent = pid_match.group(1)
                    if parent not in exons: exons[parent] = []
                    exons[parent].append((start, end))

    # Calculate Introns
    introns = {}
    print("Calculating intron coordinates...")
    
    for gid, gene_info in genes.items():
        if gid not in exons: continue
        
        chrom = gene_info['chrom']
        strand = gene_info['strand']
        
        # Sort exons by genomic position
        sorted_exons = sorted(exons[gid], key=lambda x: x[0])
        
        current_pos = sorted_exons[0][1] # End of first exon
        
        if chrom not in introns: introns[chrom] = []
        
        for i in range(1, len(sorted_exons)):
            exon_start = sorted_exons[i][0]
            
            # Gap validation
            if exon_start > current_pos + 10:
                # Apply 30bp trimming as per your methods
                intron_start = current_pos + 30
                intron_end = exon_start - 30
                
                if intron_end > intron_start:
                    # Store strand with the interval
                    introns[chrom].append((intron_start, intron_end, strand))
            
            current_pos = sorted_exons[i][1]
            
    # Sort intervals for binary search
    for chrom in introns:
        introns[chrom].sort(key=lambda x: x[0])
        
    print(f"Introns parsed for {len(introns)} chromosomes.")
    return introns

def find_intron_info(chrom, pos, intron_dict):
    """
    Binary search to check if pos is in intron list for chrom.
    Returns (True, strand) or (False, None).
    """
    if chrom not in intron_dict: return False, None
    regions = intron_dict[chrom]
    
    low = 0
    high = len(regions) - 1
    
    while low <= high:
        mid = (low + high) // 2
        start, end, strand = regions[mid]
        
        if start <= pos < end:
            return True, strand
        elif pos < start:
            high = mid - 1
        else:
            low = mid + 1
    return False, None

def process_vcf(vcf_path, intron_dict):
    """
    Reads VCF, filters for introns, corrects for strand, counts G and C alleles.
    """
    print(f"Processing VCF: {vcf_path}")
    
    # Output files now represent "Coding Strand Target"
    out_g = open('intron_counts_coding_G.csv', 'w')
    out_c = open('intron_counts_coding_C.csv', 'w')
    
    writer_g = csv.writer(out_g)
    writer_c = csv.writer(out_c)
    
    header = ['chrom', 'pos', 'k', 'n', 'strand']
    writer_g.writerow(header)
    writer_c.writerow(header)
    
    opener = gzip.open if vcf_path.endswith('.gz') else open
    
    count_sites = 0
    
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1]) - 1 # 0-based
            
            # 1. Spatial Filter + Get Strand
            in_intron, strand = find_intron_info(chrom, pos, intron_dict)
            if not in_intron:
                continue
            
            ref_genomic = parts[3]
            alt_genomic = parts[4]
            
            # 2. Quality Filter
            if len(ref_genomic) > 1 or len(alt_genomic) > 1 or ',' in alt_genomic:
                continue
                
            # 3. Get Counts (Genomic Strand)
            try:
                format_keys = parts[8].split(':')
                gt_idx = format_keys.index('GT')
            except ValueError:
                continue
            
            ref_count = 0
            alt_count = 0
            
            for s in parts[9:]:
                gt_str = s.split(':')[gt_idx]
                # Fast parsing for 0/0, 0/1, 1/1
                if '0' in gt_str: ref_count += gt_str.count('0')
                if '1' in gt_str: alt_count += gt_str.count('1')
            
            total_n = ref_count + alt_count
            if total_n < 10: continue
            
            # 4. Strand Correction
            # We want to know: What are the alleles on the CODING strand?
            if strand == '+':
                ref_coding = ref_genomic
                alt_coding = alt_genomic
                # Counts map directly
                count_of_ref_allele = ref_count
                count_of_alt_allele = alt_count
            else: # Negative strand
                # The "Ref" on coding strand is complement of genomic Ref
                ref_coding = get_complement(ref_genomic)
                alt_coding = get_complement(alt_genomic)
                # Counts still map to Ref/Alt indices, but the physical bases swapped
                count_of_ref_allele = ref_count
                count_of_alt_allele = alt_count

            # 5. Assign to Targets (Coding Strand Perspective)
            
            # -- Target G (Coding) --
            k_g = 0
            is_relevant_g = False
            
            if ref_coding == 'G':
                k_g = count_of_ref_allele
                is_relevant_g = True
            elif alt_coding == 'G':
                k_g = count_of_alt_allele
                is_relevant_g = True
            
            if is_relevant_g:
                writer_g.writerow([chrom, pos, k_g, total_n, strand])
                
            # -- Target C (Coding) --
            k_c = 0
            is_relevant_c = False
            
            if ref_coding == 'C':
                k_c = count_of_ref_allele
                is_relevant_c = True
            elif alt_coding == 'C':
                k_c = count_of_alt_allele
                is_relevant_c = True
                
            if is_relevant_c:
                writer_c.writerow([chrom, pos, k_c, total_n, strand])
                
            count_sites += 1
            if count_sites % 100000 == 0:
                print(f"Processed {count_sites} intronic sites...")

    out_g.close()
    out_c.close()
    print("Done.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_intron_sfs_stranded.py <gff_file> <vcf_file>")
        sys.exit(1)
        
    gff_file = sys.argv[1]
    vcf_file = sys.argv[2]
    
    introns = parse_gff_introns(gff_file)
    process_vcf(vcf_file, introns)