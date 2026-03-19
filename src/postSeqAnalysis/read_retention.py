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
    st.markdown("## All samples overview")

    def retention_color(pct):
        if pct >= 50: return "#378ADD"
        if pct >= 30: return "#BA7517"
        return "#E24B4A"

    def build_overview_table(df):
        header_cells = "".join(
            f'<th style="padding:8px 12px;font-weight:500;color:gray;border-bottom:0.5px solid #e0e0e0;">{s}</th>'
            for s in ["Library", "Barcode"] + STEPS + ["Reads Left (%)"]
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
                    loss_str = f'<br><span style="color:#E24B4A;">▼ {loss_pct:.1f}%</span>'
                cells += f'<td style="padding:8px 12px;white-space:nowrap;">{count:,}{loss_str}</td>'

            ret = row["Retention %"]
            color = retention_color(ret)
            cells += f"""
                <td style="padding:8px 12px;">
                    <div style="display:flex;align-items:center;gap:6px;">
                        <div style="width:60px;height:8px;background:#f0f0f0;border-radius:2px;">
                            <div style="width:{ret}%;height:100%;background:{color};border-radius:2px;"></div>
                        </div>
                        <span style="color:{color};font-weight:500;">{ret}%</span>
                    </div>
                </td>
            """

            rows_html += f"<tr>{cells}</tr>"

        return f"""
        <div style="border:0.5px solid #e0e0e0;border-radius:10px;overflow:auto;">
          <table style="width:100%;border-collapse:collapse;font-size:16px;text-align:center;">
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
        "#378ADD", "#BA7517", "#BA7517", "#BA7517",
        "#533AB7", "#1D9E75","#1D9E75"
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