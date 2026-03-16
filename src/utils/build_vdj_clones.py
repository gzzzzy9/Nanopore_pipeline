import numpy as np
import pandas as pd
import os
import argparse

def build_vdj_clones(result_dir:str, library:str, barcode:str, isotype:str, output:str='vdj_clone_info.csv'):
    """Build VDJ clones according to MiXCR output clones.

    Args:
        result_dir (str): Directory where the MiXCR put its results in. Batch dir is incorporated.
        library (str): Library ID
        barcode (str): Barcode ID
        isotype (str): Isotypes such as IgM, IgG, IgK, IgL, TRA, TRB
        output (str, optional): Path to output the CSV file. Defaults to 'vdj_clone_info.csv'.
    """
    
    if isotype in ['IgM','IgG']:
        assemble_df = pd.read_csv(os.path.join(result_dir, library, barcode, 'mixcr_assemble_IGH.tsv'),sep='\t')
    else:
        assemble_df = pd.read_csv(os.path.join(result_dir, library, barcode, 'mixcr_assemble_{}.tsv'.format(isotype)),sep='\t')
        

    assemble_df = assemble_df.loc[
        assemble_df['isProductiveFR1'] & assemble_df['isProductiveFR2'] & assemble_df['isProductiveFR3'] & 
        assemble_df['isProductiveCDR1'] & assemble_df['isProductiveCDR2'] & assemble_df['isProductiveCDR3']
    ]

    assemble_df = assemble_df.loc[
        (assemble_df['aaMutationsCountFR1Substitutions']!='region_not_covered') & (assemble_df['aaMutationsCountFR2Substitutions']!='region_not_covered') & 
        (assemble_df['aaMutationsCountFR3Substitutions']!='region_not_covered') & (assemble_df['aaMutationsCountCDR1Substitutions']!='region_not_covered') & 
        (assemble_df['aaMutationsCountCDR2Substitutions']!='region_not_covered')
    ]
    
    # In case assemble_df is empty after filtering 
    if assemble_df.empty: 
        return None

    assemble_df['nt_SHM'] = assemble_df['nMutationsCountFR1Substitutions'].astype(int) + \
        assemble_df['nMutationsCountFR2Substitutions'].astype(int) + \
        assemble_df['nMutationsCountFR3Substitutions'].astype(int) + \
        assemble_df['nMutationsCountCDR1Substitutions'].astype(int) + \
        assemble_df['nMutationsCountCDR2Substitutions'].astype(int)
        
    assemble_df['aa_SHM'] = assemble_df['aaMutationsCountFR1Substitutions'].astype(int) + \
        assemble_df['aaMutationsCountFR2Substitutions'].astype(int) + \
        assemble_df['aaMutationsCountFR3Substitutions'].astype(int) + \
        assemble_df['aaMutationsCountCDR1Substitutions'].astype(int) + \
        assemble_df['aaMutationsCountCDR2Substitutions'].astype(int)
        
    f1_f3_length = assemble_df['refPoints'].apply(lambda x:int(x.split(':')[9]) - int(x.split(':')[4]))
    assemble_df['nt_SHM_rate'] = assemble_df['nt_SHM'] / f1_f3_length
    assemble_df['aa_SHM_rate'] = assemble_df['aa_SHM'] / (f1_f3_length//3)

    # Aggregate 
    assemble_df['nt_SHM_weighted'] = assemble_df['nt_SHM'] * assemble_df['uniqueMoleculeCount']
    assemble_df['nt_SHM_rate_weighted'] = assemble_df['nt_SHM_rate'] * assemble_df['uniqueMoleculeCount']
    assemble_df['aa_SHM_weighted'] = assemble_df['aa_SHM'] * assemble_df['uniqueMoleculeCount']
    assemble_df['aa_SHM_rate_weighted'] = assemble_df['aa_SHM_rate'] * assemble_df['uniqueMoleculeCount']

    vdj_clone_df = assemble_df.groupby(['bestVHit','bestJHit','aaSeqCDR3']).agg(
        readCount=('readCount','sum'),
        uniqueMoleculeCount=('uniqueMoleculeCount','sum'),
        targetSequences=('targetSequences', lambda x: x.value_counts().index[0] if not x.empty else None),
        nt_SHM_weighted_sum=('nt_SHM_weighted','sum'),
        nt_SHM_rate_weighted_sum=('nt_SHM_rate_weighted','sum'),
        aa_SHM_weighted_sum=('aa_SHM_weighted','sum'),
        aa_SHM_rate_weighted_sum=('aa_SHM_rate_weighted','sum'),
    )
    
    # Calculate weighted SHM rate
    vdj_clone_df['nt_SHM'] = vdj_clone_df['nt_SHM_weighted_sum'] / vdj_clone_df['uniqueMoleculeCount']
    vdj_clone_df['nt_SHM_rate'] = vdj_clone_df['nt_SHM_rate_weighted_sum'] / vdj_clone_df['uniqueMoleculeCount']
    vdj_clone_df['aa_SHM'] = vdj_clone_df['aa_SHM_weighted_sum'] / vdj_clone_df['uniqueMoleculeCount']
    vdj_clone_df['aa_SHM_rate'] = vdj_clone_df['aa_SHM_rate_weighted_sum'] / vdj_clone_df['uniqueMoleculeCount']
    
    vdj_clone_df = vdj_clone_df[['targetSequences','readCount','uniqueMoleculeCount','nt_SHM','nt_SHM_rate','aa_SHM','aa_SHM_rate']]
    vdj_clone_df = vdj_clone_df.sort_values('uniqueMoleculeCount',ascending=False).reset_index()
    
    vdj_clone_df.to_csv(output, index=False)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build VDJ clones from MiXCR assembly results.")
    
    # 定义必需参数
    parser.add_argument("--result-dir", required=True, help="Base directory of MiXCR results.")
    parser.add_argument("--library", required=True, help="Library ID.")
    parser.add_argument("--barcode", required=True, help="Barcode ID.")
    parser.add_argument("--isotype", required=True, help="Isotype (e.g., IgM, IgG, TRB).")
    
    # 定义可选参数
    parser.add_argument("--output", default="vdj_clone_info.csv", help="Output CSV path (default: vdj_clone_info.csv).")

    args = parser.parse_args()

    # 执行核心逻辑
    build_vdj_clones(
        result_dir=args.result_dir,
        library=args.library,
        barcode=args.barcode,
        isotype=args.isotype,
        output=args.output
    )