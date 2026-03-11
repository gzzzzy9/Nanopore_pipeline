# mixcr_parser.py
import pandas as pd # type: ignore
import re
import argparse
import glob
import warnings
from typing import List

def extract_data_from_slurms_out(batch:str, files:str|List[str], output:str='mixcr_report.csv'):
    """ Basic statistics extraction from MiXCR output

    Args:
        batch (str): Batch name to target to the result directory.
        files (List[str]): Files containing the MiXCR stdout logs.
        output (str, optional): Output CSV file name. Defaults to 'mixcr_report.csv'.

    Returns:
        pd.DataFrame: dataframe containing the extracted statistics, including number of reads at each step and the number of clones. 
                      Results will be saved under the current dir by default.
    """
    text = ''
    if isinstance(files,str):
        files = [files]
    file_not_found_flag = True
    for file in files:
        try:
            with open(file, 'r') as f:
                text += f.read()
            file_not_found_flag = False
        except FileNotFoundError:
            warnings.warn(f'File {file} not found. Skipping.')
            continue
    if file_not_found_flag:
        return

    sample_pattern = re.compile(r"Processing (\S+) - (\S{4}) - (\S{3})")
    groups = sample_pattern.findall(text)
    libraries = pd.Series([group[0] for group in groups], name='Library')
    barcodes = pd.Series([group[1] for group in groups], name='Barcode')
    isotypes = pd.Series([group[2] for group in groups], name='Isotype')

    raw_reads_pattern = re.compile(r"Raw reads: (\d+)")
    raw_reads = pd.Series(raw_reads_pattern.findall(text)).astype(int)

    reads_after_length_filter_pattern = re.compile(r"Length filtered reads: (\d+)")
    reads_after_length_filter = pd.Series(reads_after_length_filter_pattern.findall(text)).astype(int)

    reads_after_MiXCR_filter_pattern = re.compile(r"MiXCR filtered reads: (\d+)")
    reads_after_MiXCR_filter = pd.Series(reads_after_MiXCR_filter_pattern.findall(text)).astype(int)

    reads_after_primer_filter_pattern = re.compile(r"Primer filtered reads: (\d+)")
    reads_after_primer_filter = pd.Series(reads_after_primer_filter_pattern.findall(text)).astype(int)

    reads_aligned_pattern = re.compile(r"Successfully aligned reads: (\d+)")
    reads_aligned = pd.Series(reads_aligned_pattern.findall(text)[1::2]).astype(int)

    reads_tagrefined_pattern = re.compile(r"Number of output records: (\d+)")
    reads_tagrefined = pd.Series(reads_tagrefined_pattern.findall(text)).astype(int)

    reads_assembled_pattern = re.compile(r"Reads used in clonotypes, percent of total: (\d+)")
    reads_assembled = pd.Series(reads_assembled_pattern.findall(text)).astype(int)

    clones_pattern = re.compile(r"Final clonotype count: (\d+)")
    clones = pd.Series(clones_pattern.findall(text)).astype(int)
    
    pd.DataFrame({
        'Batch':batch,
        'Library': libraries,
        'Barcode': barcodes,
        'Isotype': isotypes,
        'Raw Reads': raw_reads,
        'Reads After Length Filter': reads_after_length_filter,
        'Reads After MiXCR Filter': reads_after_MiXCR_filter,
        'Reads After Primer Filter': reads_after_primer_filter,
        'Reads Aligned': reads_aligned,
        'Reads TagRefined': reads_tagrefined,
        'Reads Assembled': reads_assembled,
        'Clones': clones,
    }).to_csv(output, index=False)
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse MiXCR stdout logs and extract basic statistics into a CSV report')
    parser.add_argument('--batch', required=True, help='Batch name')
    parser.add_argument('--input', nargs='+', required=True, help='input MiXCR stdout log files (supports glob patterns, e.g., *.out)')
    parser.add_argument('--output', default='report.csv', help='Saved CSV file name')
    
    args = parser.parse_args()
    
    # Process glob patterns such as *.out
    files = []
    for i in args.input:
        files.extend(glob.glob(i))
    extract_data_from_slurms_out(args.batch, files, args.output)