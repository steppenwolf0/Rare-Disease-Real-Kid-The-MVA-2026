import gzip

VCF_FILE = "WGS_EX2312012_HGWCNDSX7.vcf.gz"

    # "BUB1B": {"chrom": "chr15", "start": 40161069, "end": 40221123},
    # "CEP57": {"chrom": "chr11", "start": 95790014, "end": 95832693},
    # "TRIP13": {"chrom": "chr5", "start": 892884, "end": 919348},

START = 892884
END = 919348
TARGET_CHROMS = {"5", "chr5"}

with gzip.open(VCF_FILE, "rt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")
        vcf_chrom = fields[0]
        pos = int(fields[1])

        # Chromosome check
        if vcf_chrom not in TARGET_CHROMS:
            continue

        # Position check
        if pos < START or pos > END:
            continue

        ref = fields[3]
        alt = fields[4]

        # Calculate end position based on REF length
        start_pos = pos
        end_pos = pos + len(ref) - 1

        # Print formatted output: CHROM, START, END, REF/ALT
        print(f"{vcf_chrom:<10} {start_pos:<10} {end_pos:<10} {ref}/{alt}")