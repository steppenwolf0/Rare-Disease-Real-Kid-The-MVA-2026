import csv
import os
import gzip
import vcf

# 1. Define paths and metadata
VCF_FILE = "WGS_EX2312012_HGWCNDSX7.vcf.gz"
OUTPUT_CSV = "steppenwolf0_mva_candidates.csv"
PROBAND_ID = "PROBAND01"

# 2. Exact GRCh38 coordinates for our high-probability target genes
CANDIDATE_GENES = {
    "BUB1B": {"chrom": "chr15", "start": 40161069, "end": 40221123},
    "CEP57": {"chrom": "chr11", "start": 95790014, "end": 95832693},
    "TRIP13": {"chrom": "chr5", "start": 892884, "end": 919348},
    "MAD1L1": {"chrom": "chr7", "start": 1815787, "end": 2233243},


}


def scan_vcf():
    # Dictionary to pool found variants by gene to check for compound-hets
    variants_by_gene = {gene: [] for gene in CANDIDATE_GENES}

    print("🧬 Scanning target genomic regions for rare disease variants using PyVCF3...")

    # Open the file in binary mode ('rb') and decode the stream line-by-line
    with gzip.open(VCF_FILE, "rb") as gz_file:
        # Wrap the binary stream so PyVCF3 reads it as strings cleanly
        text_stream = (line.decode("utf-8") for line in gz_file)
        vcf_reader = vcf.Reader(text_stream)

        for record in vcf_reader:
            chrom = str(record.CHROM)
            pos = int(record.POS)

            # Match against target genes
            for gene, region in CANDIDATE_GENES.items():
                # Normalize chromosome names (handles both '15' and 'chr15')
                target_chrom = region["chrom"]
                target_chrom_alt = target_chrom.replace("chr", "")

                if chrom == target_chrom or chrom == target_chrom_alt:
                    # Check if the variant position falls inside the window
                    if region["start"] <= pos <= region["end"]:
                        
                        # Grab the target proband sample
                        try:
                            sample = record.genotype(PROBAND_ID)
                        except KeyError:
                            sample = None
                        
                        # Fallback to the first sample if ID name doesn't match column headers
                        if sample is None and record.samples:
                            sample = record.samples[0]

                        if not sample:
                            continue

                        # Extract indices and sequence names safely
                        gt_indices = sample.gt_alleles
                        gt_string = sample.gt_type  # 0=HOM_REF, 1=HET, 2=HOM_ALT

                        # Skip wildtype (0/0), missing, or empty records
                        if gt_string == 0 or not gt_indices or None in gt_indices:
                            continue

                        # Parse alt allele string formats cleanly
                        alt_allele = str(record.ALT[0]) if record.ALT else "N"

                        # Build the unified tracking entry
                        norm_chrom_name = chrom if chrom.startswith("chr") else f"chr{chrom}"
                        variants_by_gene[gene].append(
                            {
                                "chrom": norm_chrom_name,
                                "pos": pos,
                                "ref": str(record.REF),
                                "alt": alt_allele,
                                "gt": f"{gt_indices[0]}/{gt_indices[1]}",
                            }
                        )

    # 3. Format into Hackathon rows based on your EPCR criteria
    submission_rows = []

    for gene, vars_list in variants_by_gene.items():
        if not vars_list:
            continue

        print(f" Found {len(vars_list)} alternative variants in {gene}")

        # Check if it's a compound heterozygote (multiple different 0/1 variants)
        het_variants = [v for v in vars_list if v["gt"] in ("0/1", "1/0", "0|1", "1|0")]
        hom_variants = [v for v in vars_list if v["gt"] in ("1/1", "1|1")]

        # Scenario A: We found a compound-het pair (highly likely for MVA!)
        if len(het_variants) >= 2:
            # Take the first two heterozygous variants as a pair
            v1, v2 = het_variants[0], het_variants[1]
            submission_rows.append(
                {
                    "proband_id": PROBAND_ID,
                    "chrom_1": v1["chrom"],
                    "pos_1": v1["pos"],
                    "ref_1": v1["ref"],
                    "alt_1": v1["alt"],
                    "chrom_2": v2["chrom"],
                    "pos_2": v2["pos"],
                    "ref_2": v2["ref"],
                    "alt_2": v2["alt"],
                    "epcr": 0.95,  # Strong candidate pair
                    "finding_type": "primary",
                    "notes": f"Potential compound heterozygous causal pair in {gene} matching MVA symptoms.",
                }
            )
        # Scenario B: We found a homozygous mutation
        elif hom_variants:
            v = hom_variants[0]
            submission_rows.append(
                {
                    "proband_id": PROBAND_ID,
                    "chrom_1": v["chrom"],
                    "pos_1": v["pos"],
                    "ref_1": v["ref"],
                    "alt_1": v["alt"],
                    "chrom_2": "",
                    "pos_2": "",
                    "ref_2": "",
                    "alt_2": "",
                    "epcr": 0.90,
                    "finding_type": "primary",
                    "notes": f"Homozygous mutation in clinical candidate gene {gene}.",
                }
            )
        # Scenario C: Just single heterozygous variants
        else:
            for v in vars_list:
                submission_rows.append(
                    {
                        "proband_id": PROBAND_ID,
                        "chrom_1": v["chrom"],
                        "pos_1": v["pos"],
                        "ref_1": v["ref"],
                        "alt_1": v["alt"],
                        "chrom_2": "",
                        "pos_2": "",
                        "ref_2": "",
                        "alt_2": "",
                        "epcr": 0.50,
                        "finding_type": "secondary",
                        "notes": f"Single heterozygous variant in {gene}.",
                    }
                )

    # 4. Sort predictions by EPCR (Highest probability first)
    submission_rows = sorted(
        submission_rows, key=lambda x: x["epcr"], reverse=True
    )

    # Limit to top 10 rows per hackathon rules
    submission_rows = submission_rows[:10]

    # Write out the submission file
    fields = [
        "proband_id",
        "chrom_1",
        "pos_1",
        "ref_1",
        "alt_1",
        "chrom_2",
        "pos_2",
        "ref_2",
        "alt_2",
        "epcr",
        "finding_type",
        "notes",
    ]
    with open(OUTPUT_CSV, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(submission_rows)

    print(
        f"\n Submission file created successfully: {os.path.abspath(OUTPUT_CSV)}"
    )


if __name__ == "__main__":
    scan_vcf()
