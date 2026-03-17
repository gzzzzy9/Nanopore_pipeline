import streamlit as st
import subprocess
import pandas as pd
import plotly.express as px
import os
import re
import time
import re
import pandas as pd
import plotly.graph_objects as go

st.set_page_config(page_title="Nanopore4RepSeq UI", layout="wide")
st.title("🧬 Nanopore4RepSeq Console")
st.markdown("---")

# --- 侧边栏配置 ---
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

def parse_log_to_df(log_path: str, sample_table_path: str) -> pd.DataFrame:
    """
    读取 sample_table.csv 作为基准表，
    再用日志更新每行的时间戳和状态。
    """
    # 读取样本表作为基准
    df = pd.read_csv(sample_table_path, header=0)
    df.columns = ["library", "barcode", "isotype", "species"]
    df["time"] = "—"
    df["status"] = "pending"

    if not os.path.exists(log_path):
        return df

    # 解析日志中已处理的行
    processed = {}  # key: (library, barcode)
    with open(log_path, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        m = LOG_PATTERN.search(line.strip())
        if m:
            timestamp, library, barcode, isotype, species = m.groups()
            # 只取时间部分 HH:MM:SS
            time_short = re.search(r"\d{2}:\d{2}:\d{2}", timestamp).group()
            key = (library, barcode)
            # 判断是否是最后一条匹配行（即当前正在跑的）
            is_last = not any(
                LOG_PATTERN.search(l.strip()) for l in lines[i+1:]
            )
            processed[key] = {
                "time": time_short,
                "status": "running" if is_last else "done"
            }

    # 更新 df
    for idx, row in df.iterrows():
        key = (row["library"], row["barcode"])
        if key in processed:
            df.at[idx, "time"] = processed[key]["time"]
            df.at[idx, "status"] = processed[key]["status"]

    return df


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
    <table style="width:100%;border-collapse:collapse;font-size:13px;">
        <thead>
            <tr style="border-bottom:0.5px solid #e0e0e0;">
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;font-size:12px;">Time</th>
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;font-size:12px;">Library</th>
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;font-size:12px;">Barcode</th>
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;font-size:12px;">Isotype</th>
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;font-size:12px;">Species</th>
                <th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;font-size:12px;">Status</th>
            </tr>
        </thead>
        <tbody>{rows_html}</tbody>
    </table>
    </div>
    """
    st.markdown(html, unsafe_allow_html=True)

def tail_slurm_log(log_path, job_id):
    """轮询读取 Slurm 日志，实时展示匹配行，直到 job 结束"""
    st.info(f"📡 监听日志: `{log_path}`  |  Job ID: `{job_id}`")
    
    table_placeholder  = st.empty()
    status_placeholder = st.empty()

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
                    # 每匹配一行就重新解析整个日志刷新表格
                    df = parse_log_to_df(log_path, sample_table_path)
                    with table_placeholder.container():
                        render_log_table(df)
            else:
                check = subprocess.run(
                    f"squeue -j {job_id} -h",
                    shell=True, capture_output=True, text=True
                )
                if check.stdout.strip() == "":
                    # job 结束，最后刷新一次把 running 改为 done
                    df = parse_log_to_df(log_path, sample_table_path)
                    df.loc[df["status"] == "running", "status"] = "done"
                    with table_placeholder.container():
                        render_log_table(df)
                    status_placeholder.success("✅ Job 已完成！")
                    break
                time.sleep(3)

st.subheader("📊 Main")
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
            tail_slurm_log(log_path, job_id)
        else:
            st.error(f"❌ Fail：{res.stderr}")

st.markdown("---")
st.subheader("📊 Results")

action = st.pills(
    "Select analysis",
    options=[
        "🔍 Read Retention",
        "📊 Clone Summary",
        "🧬 SHM Distribution",
        "🏆 Top Clones",
        "📉 V Gene Usage",
        "🔬 CDR3 Length",
    ],
    selection_mode="single"
)

st.divider()

# ── 各分析模块 ──────────────────────────────────────────

def show_read_retention():
    summary_path = os.path.join(output_dir, "summary.csv")

    if not os.path.exists(summary_path):
        with st.spinner("正在生成摘要..."):
            res = subprocess.run(
                f"python src/utils/parsers.py "
                f"--batch {batch} "
                f"--input {log_path} "
                f"--output {summary_path}",
                shell=True, capture_output=True, text=True
            )
        if res.returncode != 0:
            st.error(f"❌ 生成失败：{res.stderr}")
            return

    df = pd.read_csv(summary_path)
    df["Retention %"] = (df["Reads Assembled"] / df["Raw Reads"] * 100).round(1)

    STEPS = [
        "Raw Reads",
        "Reads After Length Filter",
        "Reads After MiXCR Filter",
        "Reads After Primer Filter",
        "Reads Aligned",
        "Reads TagRefined",
        "Reads Assembled",
    ]

    # ── 方案一：多样本总览表 ──────────────────────────
    st.markdown("##### All samples overview")

    def retention_color(pct):
        if pct >= 50: return "#378ADD"
        if pct >= 30: return "#BA7517"
        return "#E24B4A"

    def build_overview_table(df):
        header_cells = "".join(
            f'<th style="text-align:left;padding:8px 12px;font-weight:500;color:gray;font-size:12px;border-bottom:0.5px solid #e0e0e0;">{s}</th>'
            for s in ["Library", "Barcode"] + STEPS + ["Retention %"]
        )

        rows_html = ""
        for _, row in df.iterrows():
            cells = f"""
                <td style="padding:8px 12px;font-weight:500;">{row['Library']}</td>
                <td style="padding:8px 12px;">{row['Barcode']}</td>
            """

            counts = [int(row[s]) for s in STEPS]
            for i, (step, count) in enumerate(zip(STEPS, counts)):
                loss_str = ""
                if i > 0:
                    prev = counts[i - 1]
                    loss_pct = (prev - count) / prev * 100
                    loss_str = f'<br><span style="color:#E24B4A;font-size:11px;">▼ {loss_pct:.1f}%</span>'
                cells += f'<td style="padding:8px 12px;white-space:nowrap;">{count:,}{loss_str}</td>'

            ret = row["Retention %"]
            color = retention_color(ret)
            cells += f"""
                <td style="padding:8px 12px;">
                    <div style="display:flex;align-items:center;gap:6px;">
                        <div style="width:60px;height:8px;background:#f0f0f0;border-radius:2px;">
                            <div style="width:{ret}%;height:100%;background:{color};border-radius:2px;"></div>
                        </div>
                        <span style="color:{color};font-weight:500;font-size:12px;">{ret}%</span>
                    </div>
                </td>
            """

            rows_html += f"<tr>{cells}</tr>"

        return f"""
        <div style="border:0.5px solid #e0e0e0;border-radius:10px;overflow:auto;">
          <table style="width:100%;border-collapse:collapse;font-size:13px;">
            <thead><tr>{header_cells}</tr></thead>
            <tbody>{rows_html}</tbody>
          </table>
        </div>
        """

    st.markdown(build_overview_table(df), unsafe_allow_html=True)

    st.divider()

    # ── 方案二：step-down bar chart ───────────────────
    st.markdown("##### Step-down view")

    col1, col2 = st.columns(2)
    with col1:
        selected_library = st.selectbox(
            "Library",
            options=sorted(df["Library"].unique().tolist()),
            key="retention_library"
        )
    with col2:
        filtered_barcodes = df[df["Library"] == selected_library]["Barcode"].tolist()
        selected_barcode = st.selectbox(
            "Barcode",
            options=filtered_barcodes,
            key="retention_barcode"
        )

    row = df[(df["Library"] == selected_library) & (df["Barcode"] == selected_barcode)].iloc[0]
    counts = [int(row[s]) for s in STEPS]

    colors = [
        "#378ADD", "#1D9E75", "#1D9E75",
        "#BA7517", "#BA7517", "#533AB7", "#533AB7"
    ]

    fig_step = go.Figure()
    for i, (label, count) in enumerate(zip(STEPS, counts)):
        loss_pct = f"▼ {((counts[i-1] - count) / counts[i-1] * 100):.1f}%" if i > 0 else ""
        fig_step.add_trace(go.Bar(
            name=label,
            x=[count],
            y=[label],
            orientation="h",
            marker_color=colors[i],
            marker_opacity=0.75,
            text=f"{count:,}  {loss_pct}",
            textposition="outside",
            showlegend=False,
        ))

    fig_step.update_layout(
        xaxis_title="Read count",
        yaxis=dict(autorange="reversed"),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(t=20, b=20, r=160),
        height=320,
        bargap=0.35,
    )
    st.plotly_chart(fig_step, use_container_width=True)

def show_clone_summary():
    example_csv = os.path.join(output_dir, batch, "vdj_clone_info.csv")
    if not os.path.exists(example_csv):
        st.info("暂无数据，请先完成分析。")
        return
    df = pd.read_csv(example_csv)
    st.write("Top 10 Clones (by UMI Count)")
    st.dataframe(df.head(10), use_container_width=True)

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

# ── 路由 ───────────────────────────────────────────────

if action == "🔍 Read Retention":
    show_read_retention()
elif action == "📊 Clone Summary":
    show_clone_summary()
elif action == "🧬 SHM Distribution":
    show_shm_distribution()
elif action == "🏆 Top Clones":
    show_top_clones()
elif action == "📉 V Gene Usage":
    show_v_gene_usage()
elif action == "🔬 CDR3 Length":
    show_cdr3_length()