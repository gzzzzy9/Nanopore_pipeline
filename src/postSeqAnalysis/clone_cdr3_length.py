import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os


def show(batch, log_path, output_dir):

    # ── 样本选择器 ────────────────────────────────────
    summary_path = os.path.join(output_dir, "summary.csv")
    if not os.path.exists(summary_path):
        st.info("请先在 Read Analysis 中生成 summary.csv")
        return

    df_summary = pd.read_csv(summary_path)

    col1, col2 = st.columns(2)
    with col1:
        selected_library = st.selectbox(
            "Library",
            options=sorted(df_summary["Library"].unique().tolist()),
            key="cdr3_library"
        )
    with col2:
        filtered_barcodes = df_summary[df_summary["Library"] == selected_library]["Barcode"].tolist()
        selected_barcode = st.selectbox(
            "Barcode",
            options=filtered_barcodes,
            key="cdr3_barcode"
        )

    vdj_path = os.path.join(output_dir, selected_library, selected_barcode, "vdj_clone_info.csv")
    if not os.path.exists(vdj_path):
        st.error(f"❌ 找不到文件：{vdj_path}")
        return

    df = pd.read_csv(vdj_path)

    if "aaSeqCDR3" not in df.columns:
        st.error("❌ vdj_clone_info.csv 中没有 aaSeqCDR3 列")
        return

    # ── 计算 CDR3 长度 ────────────────────────────────
    df["CDR3_length"] = df["aaSeqCDR3"].str.len()

    st.divider()

    # ── Metric cards ──────────────────────────────────
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Total Clones", f"{len(df):,}")
    m2.metric("Median CDR3 Length", f"{df['CDR3_length'].median():.0f} aa")
    m3.metric("Min CDR3 Length", f"{df['CDR3_length'].min()} aa")
    m4.metric("Max CDR3 Length", f"{df['CDR3_length'].max()} aa")

    st.divider()

    # ── 图一：Clone 数量分布（按 CDR3 长度）─────────────
    st.markdown("##### CDR3 length distribution — clone count")

    length_counts = df.groupby("CDR3_length").size().reset_index(name="clone_count")

    fig1 = px.bar(
        length_counts,
        x="CDR3_length",
        y="clone_count",
        labels={"CDR3_length": "CDR3 aa length", "clone_count": "Clone count"},
        color_discrete_sequence=["#378ADD"],
    )
    fig1.update_layout(
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        bargap=0.2,
        xaxis=dict(dtick=1),
        height=350,
    )
    st.plotly_chart(fig1, width='stretch')

    st.divider()

    # ── 图二：UMI 加权分布 ────────────────────────────
    st.markdown("##### CDR3 length distribution — UMI weighted")

    length_umi = df.groupby("CDR3_length")["uniqueMoleculeCount"].sum().reset_index(name="umi_count")

    fig2 = px.bar(
        length_umi,
        x="CDR3_length",
        y="umi_count",
        labels={"CDR3_length": "CDR3 aa length", "umi_count": "UMI count"},
        color_discrete_sequence=["#1D9E75"],
    )
    fig2.update_layout(
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        bargap=0.2,
        xaxis=dict(dtick=1),
        height=350,
    )
    st.plotly_chart(fig2, width='stretch')

    st.divider()

    # ── 图三：多样本对比（所有 barcode）──────────────────
    st.markdown("##### CDR3 length distribution — all samples")

    all_rows = []
    for _, row in df_summary.iterrows():
        p = os.path.join(output_dir, row["Library"], row["Barcode"], "vdj_clone_info.csv")
        if not os.path.exists(p):
            continue
        tmp = pd.read_csv(p)
        tmp["CDR3_length"] = tmp["aaSeqCDR3"].str.len()
        tmp["sample"] = f"{row['Library']} - {row['Barcode']}"
        all_rows.append(tmp[["sample", "CDR3_length", "uniqueMoleculeCount"]])

    if all_rows:
        all_df = pd.concat(all_rows, ignore_index=True)
        fig3 = px.box(
            all_df,
            x="sample",
            y="CDR3_length",
            labels={"CDR3_length": "CDR3 aa length", "sample": ""},
            color="sample",
            points=False,
        )
        fig3.update_layout(
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            showlegend=False,
            xaxis_tickangle=-30,
            height=400,
        )
        st.plotly_chart(fig3, width='stretch')