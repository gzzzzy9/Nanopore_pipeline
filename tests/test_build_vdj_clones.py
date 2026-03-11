import pytest
import pandas as pd
import numpy as np
import os
import sys

# 假设你的代码在 src/utils/build_vdj_clones.py
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from utils.build_vdj_clones import build_vdj_clones

@pytest.fixture
def mock_data_env(tmp_path):
    """模拟服务器上的文件夹结构"""
    result_dir, batch, lib, bc, iso = "nanopore_results", "20250101", "Lib01", "NB01", "IgG"
    # IgG 逻辑会寻找 IGH.tsv
    path = tmp_path / result_dir / batch / lib / bc
    path.mkdir(parents=True)
    tsv_file = path / "mixcr_assemble_IGH.tsv"
    
    # 构造模拟数据
    df = pd.DataFrame({
        'bestVHit': ['V1', 'V1', 'V2'],
        'bestJHit': ['J1', 'J1', 'J2'],
        'aaSeqCDR3': ['CASR', 'CASR', 'CASS'],
        'readCount': [100, 200, 50],
        'uniqueMoleculeCount': [10, 20, 5],
        'targetSequences': ['ATGC', 'ATGC', 'GGCC'],
        'refPoints': ['0:0:0:0:10:0:0:0:0:100:0']*3, # length = 100-10 = 90
        'isProductiveFR1': [True]*3, 'isProductiveFR2': [True]*3, 'isProductiveFR3': [True]*3,
        'isProductiveCDR1': [True]*3, 'isProductiveCDR2': [True]*3, 'isProductiveCDR3': [True]*3,
        'nMutationsCountFR1Substitutions': [1]*3, 'nMutationsCountFR2Substitutions': [1]*3,
        'nMutationsCountFR3Substitutions': [1]*3, 'nMutationsCountCDR1Substitutions': [1]*3,
        'nMutationsCountCDR2Substitutions': [1]*3,
        'aaMutationsCountFR1Substitutions': [1]*3, 'aaMutationsCountFR2Substitutions': [1]*3,
        'aaMutationsCountFR3Substitutions': [1]*3, 'aaMutationsCountCDR1Substitutions': [1]*3,
        'aaMutationsCountCDR2Substitutions': [1]*3,
    })
    for col in ['FR1','FR2','FR3','CDR1','CDR2']:
        colname = 'nMutationsCount{}Substitutions'.format(col)
        df[colname] = df[colname].astype(object)
        colname = 'aaMutationsCount{}Substitutions'.format(col)
        df[colname] = df[colname].astype(object)
    df.to_csv(tsv_file, sep='\t', index=False)
    
    return {"tmp": tmp_path, "args": (result_dir, batch, lib, bc, iso)}

def test_build_vdj_clones_logic(mock_data_env):
    tmp_path = mock_data_env["tmp"]
    output_csv = tmp_path / "output.csv"
    
    # 切换当前路径到 tmp_path 模拟服务器环境
    os.chdir(tmp_path)
    
    build_vdj_clones(*mock_data_env["args"], output=str(output_csv))
    
    assert output_csv.exists()
    res = pd.read_csv(output_csv)
    
    # 验证分组聚合：V1-J1-CASR 应该合并成一行
    assert len(res) == 2
    v1_row = res[res['bestVHit'] == 'V1'].iloc[0]
    assert v1_row['uniqueMoleculeCount'] == 30 # 10 + 20
    # nt_SHM 应该是 5 列之和 = 5
    assert v1_row['nt_SHM'] == 5.0

def test_filter_region_not_covered(mock_data_env):
    """测试过滤 region_not_covered 是否有效"""
    tmp_path = mock_data_env["tmp"]
    # 修改文件，让一行失效
    tsv = os.path.join(tmp_path, * mock_data_env["args"][:-1], "mixcr_assemble_IGH.tsv")
    print(tsv)
    df = pd.read_csv(tsv, sep='\t')
    df.loc[2, 'aaMutationsCountFR1Substitutions'] = 'region_not_covered'
    df.to_csv(tsv, sep='\t', index=False)
    
    os.chdir(tmp_path)
    out = tmp_path / "out_filtered.csv"
    build_vdj_clones(*mock_data_env["args"], output=str(out))
    
    res = pd.read_csv(out)
    # V2 对应的那行被过滤了，所以只剩 V1
    assert len(res) == 1
    assert res.loc[0, 'bestVHit'] == 'V1'