import pytest
import os
import sys
import pandas as pd

# Add src to PATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from utils.parsers import extract_data_from_slurms_out

# Test
@pytest.fixture
def test_setup():
    return {
        "file": os.path.join(os.path.dirname(__file__), 'mixcr_stdout.txt'),
        "batch": "20250930",
        "output": "test_mixcr_report.csv"
    }

def test_extraction_logic(test_setup):
    """验证数据提取是否准确"""
    # Execute function
    extract_data_from_slurms_out(
        test_setup["batch"], 
        test_setup["file"], 
        output=test_setup["output"]
    )
    
    # reads csv
    df = pd.read_csv(test_setup["output"])
    
    # Assertions
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 19
    assert df.loc[0, 'Library'] == 'sample1'
    assert df.loc[0, 'Raw Reads'] == 1196851
    
    # 检查最后一行 (使用 iloc)
    assert df.iloc[-1]['Library'] == 'sample1'
    assert df.iloc[-1]['Reads TagRefined'] == 546030
    assert df.iloc[-1]['Clones'] == 828
    
    # 清理测试生成的垃圾文件（可选）
    if os.path.exists(test_setup["output"]):
        os.remove(test_setup["output"])

def test_file_not_found(test_setup, tmp_path):
    """
    测试文件不存在时：
    1. 是否触发警告
    2. 是否【没有】错误地生成输出文件
    """
    # 1. 定义一个临时路径，防止污染你的实际项目目录
    temp_output = tmp_path / "should_not_exist.csv"
    
    # 2. 检查警告
    with pytest.warns(UserWarning, match="not found"):
        # 调用函数（假设它现在只写文件，不返回数据）
        extract_data_from_slurms_out(
            test_setup["batch"], 
            'non_existent.out', 
            output=str(temp_output)
        )
    
    # 3. 验证文件是否【没有】被创建
    # 如果文件没找到，逻辑上不应该生成空的 CSV
    assert not os.path.exists(temp_output), "当输入文件不存在时，不应创建输出文件"

def test_extraction_saves_file(test_setup, tmp_path):
    """
    测试正常情况下：文件是否成功保存且内容正确
    """
    temp_output = tmp_path / "test_report.csv"
    
    # 执行函数
    extract_data_from_slurms_out(
        test_setup["batch"], 
        test_setup["file"], 
        output=str(temp_output)
    )
    
    # 验证文件确实存在
    assert temp_output.exists()