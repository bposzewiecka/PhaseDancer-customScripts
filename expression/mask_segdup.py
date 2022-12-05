from Bio import SeqIO
import csv

DIR = '/home/basia/express/data/'

def get_segdups(fn, regions):
    
    segdups = []
    
    with open(fn) as f:
        reader = csv.DictReader(f, delimiter='\t')

        for segdup in reader:
            segdups.append(segdup)          

    return segdups

def get_regions(fn):
    
    regions = []
    
    with open(fn) as f:
        for line in f:
            line = line.split('\t')
            regions.append((line[0], int(line[1]), int(line[2])))
            
    return regions  

def mask_reference(ref_fn, regions_fn, segdups_fn, ref_masked_fn):

    ref = SeqIO.to_dict(SeqIO.parse(ref_fn, "fasta"))

    regions = get_regions(regions_fn)
    segdups = get_segdups(segdups_fn, regions)

    seq_regions = []

    for region in regions:
        seq_regions.append(ref[region[0]][region[1]:region[2]])

    ref_seqs = { name: list(str(contig.seq)) for name, contig in ref.items() }    

    for segdup in segdups:

        for region in regions:

            segdup_start = int(segdup['chromStart'])
            segdup_end = int(segdup['chromEnd'])

            region_start = region[1]
            region_end = region[2]

            if segdup['chrom'] == region[0] and (region_start < segdup_start < region_end or region_start < segdup_end < region_end):

                other_chrom = segdup['otherChrom']
                other_start = int(segdup['otherStart'])
                other_end = int(segdup['otherEnd'])

                print(other_chrom,  other_start, other_end)

                ref_seqs[other_chrom][other_start:other_end] = ['N'] * (other_end - other_start)

    for region, seq_region in zip(regions, seq_regions):
        ref_seqs[region[0]][region[1]:region[2]] = list(seq_region)

    with open(ref_masked_fn, "w") as f:

        for name, seq in ref_seqs.items():
            f.write(f'>{name}\n')
            f.write(''.join(seq))
            f.write('\n') 

mask_reference(DIR + 'refs/hg38/hg38.fa', DIR + 'frag.bed', DIR + 'segdups.tsv', DIR + 'refs/hg38_masked/hg38_masked.fa')            

