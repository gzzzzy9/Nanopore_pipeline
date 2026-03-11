import argparse
import logging
import sys
import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    stream=sys.stdout  # 把日志输出到 stdout，而不是默认的 stderr
)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter reads by MiXCR results.")
    parser.add_argument('input_tsv', help="Path to input MiXCR alignment tsv file")
    parser.add_argument('input_fastq', help="Path to input fastq file")
    parser.add_argument('output_fastq', help="Path to output fastq file")

    args = parser.parse_args()
    input_tsv = args.input_tsv
    input_fastq = args.input_fastq
    output_fastq = args.output_fastq

    mixcr_align_df = pd.read_csv(input_tsv, sep='\t', engine='python',on_bad_lines='skip')
    read_id_set = set(mixcr_align_df['readId'].tolist())
        
    with open(output_fastq, 'w') as out_handle:
        for i, record in enumerate(SeqIO.parse(input_fastq, 'fastq')):
            if i in read_id_set:
                SeqIO.write(record, out_handle, 'fastq')