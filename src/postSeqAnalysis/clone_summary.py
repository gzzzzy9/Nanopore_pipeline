import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import os, subprocess

def show(batch, log_path, output_dir):
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
    
    # ── Clone chart per sample ───────────────────
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
    
    vdj_info_file = os.path.join(output_dir,selected_library,selected_barcode,'vdj_clone_info.csv')
    if not os.path.exists(summary_path):
        st.error(f"Clone file not exist!：{res.stderr}")
    vdj_clone_df = pd.read_csv(vdj_info_file)
    vdj_clone_df['fullAASeq'] = vdj_clone_df['aaFR1Seq'] + vdj_clone_df['aaCDR1Seq'] + \
                                vdj_clone_df['aaFR2Seq'] + vdj_clone_df['aaCDR2Seq'] + \
                                vdj_clone_df['aaFR3Seq'] + vdj_clone_df['aaCDR3Seq'] + \
                                vdj_clone_df['aaFR4Seq']

    # ── Metric cards ─────────────────────────────
    m1, m2, m3= st.columns(3)
    m1.metric("Total Clones", f"{len(vdj_clone_df):,}")
    m2.metric("Total UMI", f"{vdj_clone_df['uniqueMoleculeCount'].sum():,}")
    m3.metric("Read depth", f"{(vdj_clone_df['readCount'].sum() / vdj_clone_df['uniqueMoleculeCount'].sum()):.1f}")

    st.divider()

    # ── 表格 ──────────────────────────────────────
    st.markdown("## Clone table")

    # 只展示关键列，太多列会很挤
    display_cols = [
        "cloneId", "uniqueMoleculeCount",
        "aaSeqCDR3", "bestVHit", "bestJHit",'fullAASeq'
    ]
    # 只保留存在的列，避免 KeyError
    display_cols = [c for c in display_cols if c in vdj_clone_df.columns]

    st.dataframe(
        vdj_clone_df[display_cols]
            .sort_values("uniqueMoleculeCount", ascending=False)
            .reset_index(drop=True),
        use_container_width=True,
        hide_index=True,
        column_config={
            "uniqueMoleculeCount": st.column_config.NumberColumn("UMI Count"),
            "aaSeqCDR3": st.column_config.NumberColumn("CDR3 aa sequence"),
            "bestVHit": st.column_config.NumberColumn("V usage"),
            "bestJHit": st.column_config.NumberColumn("J usage"),
            "fullAASeq": st.column_config.NumberColumn("Full aa sequence"),
        }
    )
    