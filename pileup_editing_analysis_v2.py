import pysam
from collections import defaultdict
import pandas as pd
import argparse
from tqdm import tqdm

"""
Script Description:
This script performs pileup analysis to calculate the ADAR editing level at specific positions in MIRROR sequencing data.
It uses BAM alignment files and reference FASTA sequences to compute base coverage and editing levels.

Key Features:
1. Extract flanking sequences around a target position.
2. Convert reference FASTA files to tabular format.
3. Analyze base composition (A, T, C, G) and calculate editing levels for specified positions.
4. Save results to a CSV file with detailed coverage and editing level information.

Dependencies:
- pysam
- pandas
- argparse
- tqdm
- Bio (from Biopython)

Usage:
python pileup_editing_analysis_v2.py -b sample.bam -f reference.fa -c 10000 -o ./output/ -F 10 -p 1234
example : for UAA site in our paper 
`python pileup_editing_analysis_v2.py -b ./5L2_E2E_STAR.Aligned.out.sorted.bam -f ./UAA_pool_reference.fa  -o ./output/ -F 10 -p 124`
"""

# Intercept a sequence based on position and flanking size
def intercept_seq(seq, flanking, pos):
    return seq[pos - flanking - 1:pos + flanking]

# Convert a FASTA file to a DataFrame
def fa2csv(input_fasta):
    from Bio import SeqIO
    sequences = [(record.id, str(record.seq)) for record in SeqIO.parse(input_fasta, "fasta")]
    return pd.DataFrame(sequences, columns=["IDs", "seq"])

# Calculate pileup levels for a given sample and reference
def pileup_level(sample_bam, ref_df, pos_sites, cover_max):
    def calculate_level(row):
        return round(row["G"] / row["coverage"], 3) if row["coverage"] > 0 else "NaN"

    output = defaultdict(dict)
    with pysam.AlignmentFile(sample_bam, "rb") as bam:
        for idx in tqdm(ref_df.index):
            contig = ref_df.at[idx, "IDs"]

            # Determine coverage
            # cover = sum(sum(base_counts) for base_counts in bam.count_coverage(contig=contig, start=options.pos - 1, stop=options.pos)) if cover_max.upper() == "MAX" else int(cover_max)

            ref_seq = ref_df.at[idx, "seq"]
            count = bam.count_coverage(contig=contig)

            for pos in pos_sites:
                ref_base = ref_seq[pos - 1]
                key = f"{contig}*{pos}*{ref_base}"
                output[key]["A"] = count[0][pos - 1]
                output[key]["C"] = count[1][pos - 1]
                output[key]["G"] = count[2][pos - 1]
                output[key]["T"] = count[3][pos - 1]

    # Convert to DataFrame
    output_df = pd.DataFrame.from_dict(output, orient="index")
    output_df["coverage"] = output_df.sum(axis=1)
    output_df["level"] = output_df.apply(calculate_level, axis=1)
    output_df["pos"] = output_df.index.str.split("*").str[1]
    output_df["ref"] = output_df.index.str.split("*").str[0]
    output_df["ref_base"] = output_df.index.str.split("*").str[2]

    # Save to CSV
    output_df.to_csv(f"{options.output_dir}/df_level_flanking{options.flanking}_cover_count{cover_max}.csv")
    print(f"Processed {len(output_df)} positions!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script for pileup analysis to calculate editing levels at specific positions.",
        fromfile_prefix_chars='@',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-b', '--bam', dest='sample_bam_file', required=True, help='Path to BAM file')
    parser.add_argument('-f', '--fa', dest='fa', required=True, help='Path to reference FASTA file')
    parser.add_argument('-c', '--cover', dest='cover', required=False, help='Pileup depth (MAX or int for pysam pileup)')
    parser.add_argument('-o', '--output', dest='output_dir', required=True, help='Output directory')
    parser.add_argument('-F', '--flanking', dest='flanking', required=True, help='Number of surrounding bases that need to be pileup')
    parser.add_argument('-p', '--pos', dest='pos', required=True, help='Editing site position in your reference')

    options = parser.parse_args()

    ref_df = fa2csv(input_fasta=options.fa)
    flanking = int(options.flanking)
    editing_site = int(options.pos)
    position_sites = list(range(editing_site - flanking, editing_site + flanking + 1))

    pileup_level(
        sample_bam=options.sample_bam_file,
        ref_df=ref_df,
        pos_sites=position_sites,
        cover_max=options.cover
    )
