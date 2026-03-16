import streamlit as st
import subprocess
import pandas as pd
import plotly.express as px
import os
import re
import time

st.set_page_config(page_title="Nanopore4RepSeq UI", layout="wide")
st.title("🧬 Nanopore4RepSeq 控制中心")
st.markdown("---")

# --- 侧边栏配置 ---
st.sidebar.header("📋 运行参数配置")
batch = st.sidebar.text_input("Batch ID", value="20260209")
log_path = st.sidebar.text_input("标准输出文件", value="../slurms_out/test.txt")
input_dir = st.sidebar.text_input("输入数据路径 (Data Dir)", value="../../Nanopore_data/test_data")
output_dir = st.sidebar.text_input("输出结果路径 (Results Dir)", value="../test")
species = st.sidebar.selectbox("物种 (Species)", ["Human", "Mouse"])

st.sidebar.header("⚙️ 运行模式")
mode = st.sidebar.radio("选择运行环境", ["本地/登录节点 (Bash)", "计算集群 (Slurm/Sbatch)"])

# 匹配形如: (Wed Mar  4 09:27:13 CST 2026) Processing sample1 - NB03 - IgM - Mouse
LOG_PATTERN = re.compile(r"\(.+?\) Processing .+")

# Load styles
def load_css():
    with open("src/static/style.css", "r") as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

load_css()

def tail_slurm_log(log_path, job_id):
    """轮询读取 Slurm 日志，实时展示匹配行，直到 job 结束"""
    st.info(f"📡 监听日志: `{log_path}`  |  Job ID: `{job_id}`")
    
    log_placeholder = st.empty()
    status_placeholder = st.empty()
    matched_lines = []

    # 等待日志文件出现（最多 60s）
    for _ in range(60):
        if os.path.exists(log_path):
            break
        status_placeholder.warning("⏳ 等待日志文件生成...")
        time.sleep(1)
    else:
        st.error("❌ 日志文件未出现，请检查路径或任务状态。")
        return

    with open(log_path, "r") as f:
        while True:
            line = f.readline()
            if line:
                if LOG_PATTERN.search(line.strip()):
                    matched_lines.append(line.strip())
                    # 实时更新展示框
                    log_placeholder.code("\n".join(matched_lines), language="")
            else:
                # 检查 job 是否还在运行
                check = subprocess.run(
                    f"squeue -j {job_id} -h",
                    shell=True, capture_output=True, text=True
                )
                if check.stdout.strip() == "":
                    status_placeholder.success("✅ Job 已完成！")
                    break
                time.sleep(3)  # 没有新行时等待 3s 再轮询

st.subheader("📊 主程序")
if st.button("🚀 执行分析"):
    os.makedirs("slurms_out", exist_ok=True)

    if mode == "本地/登录节点 (Bash)":
        cmd = f"bash nanopore_pipeline.sh -b {batch} -i {input_dir} -o {output_dir}"
        with st.status("正在运行流程...", expanded=True):
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            matched_lines = []
            placeholder = st.empty()
            for line in iter(process.stdout.readline, b''):
                decoded = line.decode().strip()
                if LOG_PATTERN.search(decoded):
                    matched_lines.append(decoded)
                    placeholder.code("\n".join(matched_lines), language="")
        st.success("✅ 本地运行完成！")

    else:
        cmd = f"sbatch -o {log_path} nanopore_pipeline.sh -b {batch} -i {input_dir} -o {output_dir}"
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if res.returncode == 0:
            job_id = res.stdout.strip().split()[-1]
            st.success(f"✅ 任务已提交！Job ID: `{job_id}`")
            tail_slurm_log(log_path, job_id)
        else:
            st.error(f"❌ 提交失败：{res.stderr}")

# --- 结果快览 ---
st.markdown("---")
st.subheader("📊 结果快览")

# 初始化 session state
if "summary_generated" not in st.session_state:
    st.session_state.summary_generated = False

if st.button("🔍 生成分析摘要"):
    summary_path = os.path.join(output_dir, "summary.csv")

    if not os.path.exists(log_path):
        st.error(f"❌ 日志文件不存在：{log_path}")
    else:
        with st.spinner("正在生成摘要..."):
            res = subprocess.run(
                f"python src/utils/parsers.py "
                f"--batch {batch} "
                f"--input {log_path} "
                f"--output {summary_path}",
                shell=True, capture_output=True, text=True
            )
        if res.returncode == 0:
            st.session_state.summary_generated = True
            st.toast("摘要生成完成！", icon="✅")
        else:
            st.error(f"❌ 生成失败：{res.stderr}")

# 只有点击按钮成功后才展示
if st.session_state.summary_generated:
    summary_csv = os.path.join(output_dir, "summary.csv")
    example_csv = os.path.join(output_dir, batch, "vdj_clone_info.csv")

    if os.path.exists(summary_csv):
        st.write("**分析摘要：**")
        st.dataframe(pd.read_csv(summary_csv), use_container_width=True)

    if os.path.exists(example_csv):
        df = pd.read_csv(example_csv)
        st.write("**Top 10 Clones (by UMI Count)：**")
        st.dataframe(df.head(10), use_container_width=True)
        fig = px.histogram(df, x="aa_SHM_rate", nbins=30, title="AA SHM Rate Distribution",
                           color_discrete_sequence=['indianred'])
        st.plotly_chart(fig, use_container_width=True)