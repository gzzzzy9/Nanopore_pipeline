import streamlit as st
import subprocess
import pandas as pd
import plotly.express as px
import os

st.set_page_config(page_title="Nanopore4RepSeq UI", layout="wide")

st.title("🧬 Nanopore4RepSeq 控制中心")
st.markdown("---")

# --- 侧边栏配置 ---
st.sidebar.header("📋 运行参数配置")
batch = st.sidebar.text_input("Batch ID", value="20260209")
input_dir = st.sidebar.text_input("输入数据路径 (Data Dir)", value="./Nanopore_data")
output_dir = st.sidebar.text_input("输出结果路径 (Results Dir)", value="./Nanopore_results")
species = st.sidebar.selectbox("物种 (Species)", ["Human", "Mouse"])

st.sidebar.header("⚙️ 运行模式")
mode = st.sidebar.radio("选择运行环境", ["本地/登录节点 (Bash)", "计算集群 (Slurm/Sbatch)"])

if st.button("🚀 执行分析"):
    if mode == "本地/登录节点 (Bash)":
        # 直接运行，实时获取输出
        cmd = f"bash nanopore_pipeline.sh -b {batch} -i {input_dir} -o {output_dir}"
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        # ... 实时打印逻辑 ...
        
    else:
        # 提交到 Slurm
        log_path = f"slurms_out/mixcr_{batch}.log"
        cmd = f"sbatch -o {log_path} nanopore_pipeline.sh -b {batch} -i {input_dir} -o {output_dir}"
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if res.returncode == 0:
            job_id = res.stdout.strip().split()[-1]
            st.success(f"✅ 任务已提交至集群！Job ID: {job_id}")
            st.info(f"请稍后在侧边栏刷新结果，或查看日志: {log_path}")

# --- 主界面：任务提交 ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("🚀 启动分析流程")
    if st.button("开始运行全流程 (Run Pipeline)"):
        with st.status("正在运行 MiXCR 流程...", expanded=True) as status:
            # 这里的命令对应你之前的 sbatch 脚本（可以改为直接运行 sh）
            cmd = f"bash nanopore_pipeline.sh -b {batch} -i {input_dir} -o {output_dir}"
            
            # 使用 subprocess 运行并捕获输出
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            for line in iter(process.stdout.readline, b''):
                st.code(line.decode().strip()) # 在界面实时打印日志
            
            status.update(label="分析完成！", state="complete", expanded=False)
        st.success("✅ 所有样本已处理完毕！")

# --- 主界面：结果可视化 ---
with col2:
    st.subheader("📊 结果快览")
    # 假设分析完成后，尝试读取其中一个克隆信息文件
    example_csv = os.path.join(output_dir, batch, "vdj_clone_info.csv") # 简化路径逻辑
    
    if os.path.exists(example_csv):
        df = pd.read_csv(example_csv)
        
        # 展示 Top 10 克隆
        st.write("Top 10 Clones (by UMI Count):")
        st.dataframe(df.head(10), use_container_width=True)
        
        # 画一个 SHM 突变率分布图
        fig = px.histogram(df, x="aa_SHM_rate", nbins=30, title="AA SHM Rate Distribution",
                           color_discrete_sequence=['indianred'])
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("等待分析完成后，此处将自动显示结果图表。")
