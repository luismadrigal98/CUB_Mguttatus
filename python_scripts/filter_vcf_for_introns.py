import sys
import csv
import re
import argparse
import multiprocessing
from collections import defaultdict
import time
import gzip

# --- Global Storage for Workers ---
# This dictionary will be shared with workers via inheritance (Linux) 
# or initializer (Windows/Mac) to avoid pickling overhead.
INTRON_INDEX = {}

# --- Helper Functions ---

def get_complement(nuc):
    """Return complement of a nucleotide."""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '.': '.'}
    return mapping.get(nuc, 'N')

def parse_gff_introns(gff_file):
    """
    Parses GFF3 to find gene coordinates/strand and calculates intron intervals.
    Returns: {chrom: [(start, end, strand), ...]}
    """
    print(f"Loading GFF: {gff_file}")
    genes = {}
    exons = {}
    
    # 1. Parse Genes and mRNA/CDS features
    # Note: Some GFF3 files use 'exon', others use 'CDS' - we handle both
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            chrom, feature = parts[0], parts[2]
            if "Chr" not in chrom:
                continue # Skip scaffolds/non-chromosomal
            start, end = int(parts[3]) - 1, int(parts[4]) # 0-based
            strand = parts[6]
            attributes = parts[8]
            
            if feature == 'gene':
                gid_match = re.search(r'ID=([^;]+)', attributes)
                if gid_match:
                    genes[gid_match.group(1)] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand}
            elif feature == 'mRNA':
                # Store mRNA as gene-level feature (some GFFs only have mRNA, not gene)
                mid_match = re.search(r'ID=([^;]+)', attributes)
                pid_match = re.search(r'Parent=([^;]+)', attributes)
                if mid_match:
                    mrna_id = mid_match.group(1)
                    # Link mRNA to gene parent if exists
                    if pid_match and pid_match.group(1) not in genes:
                        genes[pid_match.group(1)] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand}
            elif feature in ['exon', 'CDS']:
                # CDS features define coding exons
                pid_match = re.search(r'Parent=([^;]+)', attributes)
                if pid_match:
                    parent = pid_match.group(1)
                    if parent not in exons: exons[parent] = []
                    exons[parent].append((start, end))

    # 2. Calculate Introns
    introns = {}
    print("Calculating intron coordinates...")
    count_introns = 0
    count_genes_processed = 0
    
    # For GFF3 files with mRNA->CDS hierarchy, we need to:
    # 1. Find all mRNA IDs
    # 2. Get their parent genes
    # 3. Use CDS as exons
    
    # Build mRNA to gene mapping
    mrna_to_gene = {}
    for gid, gene_info in genes.items():
        # Gene IDs often have pattern: MgIM767.01G000100.v2.1
        # mRNA IDs often have pattern: MgIM767.01G000100.1.v2.1
        # We need to match them
        pass  # We'll handle this in the loop below
    
    for parent_id, exon_list in exons.items():
        # parent_id might be mRNA ID, we need to find the gene
        # Try to find matching gene
        gene_info = None
        
        # Strategy 1: Direct match (parent_id is gene ID)
        if parent_id in genes:
            gene_info = genes[parent_id]
            gene_id = parent_id
        else:
            # Strategy 2: Parent is mRNA, find its gene parent
            # Look for gene ID by removing transcript suffix
            # Example: MgIM767.01G000100.1.v2.1 -> MgIM767.01G000100.v2.1
            for gid in genes:
                if parent_id.startswith(gid.rsplit('.', 2)[0]):  # Match base name
                    gene_info = genes[gid]
                    gene_id = gid
                    break
        
        if not gene_info:
            continue  # Skip if we can't find gene info
        
        count_genes_processed += 1
        chrom, strand = gene_info['chrom'], gene_info['strand']
        
        # Sort exons by position
        sorted_exons = sorted(exon_list, key=lambda x: x[0])
        
        if len(sorted_exons) < 2:
            continue  # Need at least 2 exons to have an intron
        
        # Introns are gaps between exons
        current_pos = sorted_exons[0][1]
        
        if chrom not in introns: introns[chrom] = []
        
        for i in range(1, len(sorted_exons)):
            exon_start = sorted_exons[i][0]
            
            # 10bp gap check (minimal intron size)
            if exon_start > current_pos + 10:
                # 30bp trimming from exon boundaries
                intron_start = current_pos + 30
                intron_end = exon_start - 30
                
                if intron_end > intron_start:
                    introns[chrom].append((intron_start, intron_end, strand))
                    count_introns += 1
            
            current_pos = sorted_exons[i][1]
            
    # Sort for binary search
    for chrom in introns:
        introns[chrom].sort(key=lambda x: x[0])
        
    print(f"Loaded {count_introns} introns across {len(introns)} chromosomes.")
    print(f"Processed {count_genes_processed} genes with CDS/exon annotations.")
    
    if count_introns == 0:
        print("\n⚠️  WARNING: No introns found!")
        print("Possible issues:")
        print("  1. GFF3 file may not contain multi-exon genes")
        print("  2. CDS features may not be properly formatted")
        print("  3. Parent-child relationships may be incorrect")
        print("\nShowing first few genes and their exon counts:")
        for i, (parent_id, exon_list) in enumerate(list(exons.items())[:5]):
            print(f"  {parent_id}: {len(exon_list)} exons")
    
    return introns

def binary_search_intron(chrom, pos, intron_dict):
    """Returns (True, strand) if pos is in an intron, else (False, None)."""
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

# --- Worker Logic ---

def worker_init(intron_dict_shared):
    """Initializer to set the global intron index in each worker process."""
    global INTRON_INDEX
    INTRON_INDEX = intron_dict_shared

def process_batch(lines):
    """
    Process a batch of VCF lines.
    Returns two SFS dictionaries: sfs_G, sfs_C
    Key: (n, k), Value: count
    """
    local_sfs_G = defaultdict(int)
    local_sfs_C = defaultdict(int)
    
    for line in lines:
        if line.startswith('#'): continue
        
        parts = line.strip().split('\t')
        if len(parts) < 10: continue
        
        chrom = parts[0]
        try:
            pos = int(parts[1]) - 1
        except ValueError:
            continue
            
        in_intron, strand = binary_search_intron(chrom, pos, INTRON_INDEX)
        if not in_intron:
            continue
            
        ref_genomic = parts[3]
        alt_genomic = parts[4]
        is_invariant = (alt_genomic == '.' or alt_genomic == '<NON_REF>' or alt_genomic == '*')
        
        if not is_invariant and (len(ref_genomic) > 1 or len(alt_genomic) > 1 or ',' in alt_genomic):
            continue

        # --- CORRECTED LOGIC START ---
        try:
            fmt = parts[8]
            fmt_fields = fmt.split(':')
            
            try:
                gt_idx = fmt_fields.index('GT')
                ad_idx = fmt_fields.index('AD')
            except ValueError:
                # If GT or AD is missing, we can't perform this specific filter
                continue
            
            c0 = 0 # Ref allele count (Haplotypes)
            c1 = 0 # Alt allele count (Haplotypes)
            
            for sample_str in parts[9:]:
                sample_fields = sample_str.split(':')
                if len(sample_fields) <= max(gt_idx, ad_idx):
                    continue
                
                # 1. CHECK DEPTH (AD)
                ad_val = sample_fields[ad_idx]
                if ad_val == '.' or ad_val == '': continue
                
                try:
                    # Handle "0,0" or "0,0,0" etc
                    depths = [int(x) for x in ad_val.split(',') if x != '.']
                    total_depth = sum(depths)
                except ValueError:
                    continue
                
                # FILTER: If total depth is 0, this is missing data
                if total_depth == 0:
                    continue
                
                # 2. COUNT GENOTYPES (GT)
                # Only if depth > 0 do we trust the call
                gt_val = sample_fields[gt_idx]
                
                # Count alleles (0 or 1)
                # This ensures n is sample size, not read depth
                c0 += gt_val.count('0')
                c1 += gt_val.count('1')
                
        except Exception:
            continue
        # --- CORRECTED LOGIC END ---

        total_n = c0 + c1
        if total_n < 10: continue 

        # Strand Correction
        if strand == '+':
            ref_coding = ref_genomic
            alt_coding = alt_genomic
            count_ref = c0
            count_alt = c1
        else: # Negative Strand
            ref_coding = get_complement(ref_genomic)
            alt_coding = get_complement(alt_genomic)
            count_ref = c0 
            count_alt = c1

        # Populate SFS G
        k_g = -1
        if ref_coding == 'G':
            k_g = count_ref
        elif not is_invariant and alt_coding == 'G':
            k_g = count_alt
        elif is_invariant and ref_coding != 'G':
            k_g = 0
        elif not is_invariant and ref_coding != 'G' and alt_coding != 'G':
            k_g = 0
            
        if k_g != -1:
            local_sfs_G[(total_n, k_g)] += 1

        # Populate SFS C
        k_c = -1
        if ref_coding == 'C':
            k_c = count_ref
        elif not is_invariant and alt_coding == 'C':
            k_c = count_alt
        elif is_invariant and ref_coding != 'C':
            k_c = 0
        elif not is_invariant and ref_coding != 'C' and alt_coding != 'C':
            k_c = 0
            
        if k_c != -1:
            local_sfs_C[(total_n, k_c)] += 1
            
    return local_sfs_G, local_sfs_C

# --- Main Driver ---

def main():
    parser = argparse.ArgumentParser(description="Parallel SFS generation for Introns from massive VCF.")
    parser.add_argument('--vcf', type=str, help='Path to input VCF (gzip ok, but piping via --stream is faster)')
    parser.add_argument('--stream', action='store_true', help='Read VCF from stdin (Recommended: zcat file.vcf.gz | python script.py --stream)')
    parser.add_argument('--gff', type=str, required=True, help='GFF3 annotation file')
    parser.add_argument('--workers', type=int, default=max(1, multiprocessing.cpu_count() - 1), help='Number of worker processes')
    parser.add_argument('--batch_size', type=int, default=20000, help='Lines per batch to send to workers')
    
    args = parser.parse_args()
    
    # 1. Load Introns (Main Process)
    introns = parse_gff_introns(args.gff)
    
    # 2. Setup Parallel Pool
    print(f"Starting {args.workers} workers...")
    pool = multiprocessing.Pool(
        processes=args.workers,
        initializer=worker_init,
        initargs=(introns,)
    )
    
    # 3. Process Stream
    if args.stream:
        input_stream = sys.stdin
        print("Reading VCF from STDIN...")
    elif args.vcf:
        if args.vcf.endswith('.gz'):
            input_stream = gzip.open(args.vcf, 'rt')
        else:
            input_stream = open(args.vcf, 'r')
        print(f"Reading VCF from {args.vcf}...")
    else:
        print("Error: Must provide --vcf or --stream")
        sys.exit(1)
        
    # Master Aggregators
    GLOBAL_SFS_G = defaultdict(int)
    GLOBAL_SFS_C = defaultdict(int)
    
    batch_lines = []
    active_jobs = []
    total_lines_read = 0
    total_batches_sent = 0
    
    start_time = time.time()
    
    # Callback to merge results as they finish
    def update_result(result):
        sfs_g_part, sfs_c_part = result
        for k, v in sfs_g_part.items():
            GLOBAL_SFS_G[k] += v
        for k, v in sfs_c_part.items():
            GLOBAL_SFS_C[k] += v

    for line in input_stream:
        if line.startswith('##'): continue # Skip meta headers
        
        batch_lines.append(line)
        total_lines_read += 1
        
        # When batch is full, submit to pool
        if len(batch_lines) >= args.batch_size:
            pool.apply_async(process_batch, (batch_lines,), callback=update_result)
            batch_lines = []
            total_batches_sent += 1
            
            if total_batches_sent % 100 == 0:
                elapsed = time.time() - start_time
                print(f"Read {total_lines_read:,} lines... ({total_lines_read/elapsed:.0f} lines/sec)")
                
                # Prevent memory overflow if workers are slow
                # Simple backpressure: check pending jobs (optional but good for stability)
                # For simplicity in this script, we trust the OS scheduler + memory
    
    # Submit final batch
    if batch_lines:
        pool.apply_async(process_batch, (batch_lines,), callback=update_result)
    
    input_stream.close()
    
    print("Finished reading file. Waiting for workers to finish...")
    pool.close()
    pool.join()
    
    # 4. Write Results
    print("Writing output files...")
    
    with open('sfs_introns_G.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['n', 'k', 'count'])
        for (n, k), count in GLOBAL_SFS_G.items():
            w.writerow([n, k, count])
            
    with open('sfs_introns_C.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['n', 'k', 'count'])
        for (n, k), count in GLOBAL_SFS_C.items():
            w.writerow([n, k, count])
            
    print("Done!")

if __name__ == "__main__":
    main()