import streamlit as st
import subprocess
import pandas as pd
import plotly.express as px
import os
import re
import time
from postSeqAnalysis import read_retention, clone_summary

st.set_page_config(page_title="Nanopore4RepSeq UI", layout="wide")
st.title("🧬 Nanopore4RepSeq Console")
st.markdown("---")

# --- Sidebar configuration ---
# 初始化 session_state
if "batch" not in st.session_state:
    st.session_state.batch = ""
if "log_path" not in st.session_state:
    st.session_state.log_path = "../slurms_out/.txt"
if "input_dir" not in st.session_state:
    st.session_state.input_dir = "../../Nanopore_data/"
if "output_dir" not in st.session_state:
    st.session_state.output_dir = "../"

def on_batch_change():
    b = st.session_state.batch
    st.session_state.log_path = f"../slurms_out/{b}.txt"
    st.session_state.input_dir = f"../../Nanopore_data/{b}"
    st.session_state.output_dir = f"../{b}"

st.sidebar.header("📋 Params")
batch = st.sidebar.text_input("Batch ID", key="batch", on_change=on_batch_change)
log_path = st.sidebar.text_input("Out file", key="log_path")
input_dir = st.sidebar.text_input("Data Dir", key="input_dir")
output_dir = st.sidebar.text_input("Results Dir", key="output_dir")
species = st.sidebar.selectbox("Species", ["Mouse", "Human"])

st.sidebar.header("⚙️ Run mode")
mode = st.sidebar.radio("Select run mode", ["Local", "Slurms"])

# 匹配形如: (Wed Mar  4 09:27:13 CST 2026) Processing sample1 - NB03 - IgM - Mouse
LOG_PATTERN = re.compile(r"\(.+?\) Processing .+")

# Load styles
def load_css():
    with open("src/static/style.css", "r") as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

load_css()

LOG_PATTERN = re.compile(
    r"\((\w{3} \w{3} +\d{1,2} \d{2}:\d{2}:\d{2} \w+ \d{4})\) Processing (\S+) - (\S+) - (\S+) - (\S+)"
)

def render_log_table(df: pd.DataFrame):
    """用 st.markdown 渲染带颜色徽章的 HTML 表格"""

    def status_badge(status):
        if status == "done":
            return '<span style="background:#d4edda;color:#155724;padding:2px 10px;border-radius:99px;font-size:11px;">done</span>'
        elif status == "running":
            return '<span style="background:#cce5ff;color:#004085;padding:2px 10px;border-radius:99px;font-size:11px;">running</span>'
        else:
            return ""

    done_count = (df["status"] == "done").sum()
    running_count = (df["status"] == "running").sum()
    total = len(df)

    rows_html = ""
    for _, row in df.iterrows():
        opacity = "opacity:0.35;" if row["status"] == "pending" else ""
        rows_html += f"""
        <tr style="{opacity}">
            <td>{row['time']}</td>
            <td>{row['library']}</td>
            <td>{row['barcode']}</td>
            <td><code style="font-size:12px;">{row['isotype']}</code></td>
            <td>{row['species']}</td>
            <td>{status_badge(row['status'])}</td>
        </tr>"""

    html = f"""
    <div style="margin-bottom:8px;font-size:13px;color:gray;">
        {total} samples &nbsp;·&nbsp; {done_count} done &nbsp;·&nbsp; {running_count} running
    </div>
    <div style="border:0.5px solid #e0e0e0;border-radius:10px;overflow:hidden;">
    <table style="width:100%;border-collapse:collapse;font-size:16px;">
        <thead>
            <tr style="border-bottom:0.5px solid #e0e0e0;">
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;">Time</th>
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;">Library</th>
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;">Barcode</th>
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;">Isotype</th>
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;">Species</th>
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;">Status</th>
            </tr>
        </thead>
        <tbody>{rows_html}</tbody>
    </table>
    </div>
    """
    st.markdown(html, unsafe_allow_html=True)

def parse_log_to_df(sample_table_path: str) -> pd.DataFrame:
    """只读 sample_table 生成基准表，不读日志"""
    df = pd.read_csv(sample_table_path, header=0)
    df.columns = ["library", "barcode", "isotype", "species"]
    df["time"] = "—"
    df["status"] = "pending"
    return df


def tail_slurm_log(log_path, job_id, sample_table_path):
    st.info(f"📡 监听日志: `{log_path}`  |  Job ID: `{job_id}`")

    table_placeholder = st.empty()
    status_placeholder = st.empty()

    # 等待日志文件出现
    for _ in range(60):
        if os.path.exists(log_path):
            break
        status_placeholder.warning("⏳ 等待日志文件生成...")
        time.sleep(1)
    else:
        st.error("❌ 日志文件未出现")
        return

    # 基准表
    df = parse_log_to_df(sample_table_path)

    with open(log_path, "r") as f:
        last_key = None
        while True:
            line = f.readline()
            if line:
                m = LOG_PATTERN.search(line.strip())
                if m:
                    timestamp, library, barcode, isotype, species = m.groups()
                    time_short = re.search(r"\d{2}:\d{2}:\d{2}", timestamp).group()
                    key = (library, barcode)

                    # 把上一个 running 改成 done
                    if last_key and last_key != key:
                        df.loc[
                            (df["library"] == last_key[0]) & (df["barcode"] == last_key[1]),
                            "status"
                        ] = "done"

                    # 更新当前行为 running
                    mask = (df["library"] == library) & (df["barcode"] == barcode)
                    df.loc[mask, "time"] = time_short
                    df.loc[mask, "status"] = "running"
                    last_key = key

                    # 立刻刷新表格
                    with table_placeholder.container():
                        render_log_table(df)

            else:
                # 检查 job 是否结束
                check = subprocess.run(
                    f"squeue -j {job_id} -h",
                    shell=True, capture_output=True, text=True
                )
                if check.stdout.strip() == "":
                    # job 结束，把最后一个 running 改成 done
                    df.loc[df["status"] == "running", "status"] = "done"
                    with table_placeholder.container():
                        render_log_table(df)
                    status_placeholder.success("✅ Job 已完成！")
                    break
                time.sleep(3)

st.markdown("# ⚙️ Main")
if st.button("🚀 Run MiXCR pipeline"):
    os.makedirs("slurms_out", exist_ok=True)

    if mode == "Local":
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
            st.success(f"✅ Submitted! Job ID: `{job_id}`")
            sample_table_path = os.path.join(input_dir, "sample_table.csv")
            tail_slurm_log(log_path, job_id, sample_table_path)
        else:
            st.error(f"❌ Fail：{res.stderr}")

st.markdown("---")
st.markdown("# 📊 Results")

# 第一级 pills
category = st.pills(
    "Select analysis",
    options=[
        "🔍 Read Analysis",
        "📊 Clone Analysis",
    ],
    selection_mode="single"
)

def show_shm_distribution():
    example_csv = os.path.join(output_dir, batch, "vdj_clone_info.csv")
    if not os.path.exists(example_csv):
        st.info("暂无数据，请先完成分析。")
        return
    df = pd.read_csv(example_csv)
    fig = px.histogram(df, x="aa_SHM_rate", nbins=30, title="AA SHM Rate Distribution",
                       color_discrete_sequence=["indianred"])
    st.plotly_chart(fig, use_container_width=True)

def show_top_clones():
    st.info("🚧 待实现")

def show_v_gene_usage():
    st.info("🚧 待实现")

def show_cdr3_length():
    st.info("🚧 待实现")

# ── Router ───────────────────────────────────────────────

if category == "🔍 Read Analysis":
    read_retention.show(batch, log_path, output_dir)

elif category == "📊 Clone Analysis":
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📋 Summary",
        "🧬 SHM Distribution",
        "🏆 Top Clones",
        "📉 V Gene Usage",
        "🔬 CDR3 Length",
    ])
    with tab1:
        clone_summary.show(batch, log_path, output_dir)