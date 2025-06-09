import pandas as pd
import numpy as np
import plotly.graph_objects as go
import xml.etree.ElementTree as ET
import streamlit as st

from shiny_app_prot import (
    parse_path_data,
    svg_files,
    expression_data_filtered,
    expression_data_clavata_filtered,
)


def parse_and_plot(svg_files, values):
    points = []
    for svg_file, val in zip(svg_files, values):
        tree = ET.parse(svg_file)
        root = tree.getroot()
        namespace = {"svg": "http://www.w3.org/2000/svg"}
        path_data = root.find(".//svg:path", namespace).attrib["d"]
        coords = parse_path_data(path_data)
        if not coords:
            continue
        xs, ys = zip(*coords)
        centroid = (sum(xs) / len(xs), sum(ys) / len(ys))
        points.append((centroid[0], centroid[1], val))

    if not points:
        return ""

    x, y, z = zip(*points)
    fig = go.Figure(
        data=[
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="markers",
                marker=dict(size=12, color=z, colorscale="Viridis", opacity=0.8),
            )
        ]
    )
    fig.update_layout(
        scene=dict(xaxis_visible=False, yaxis_visible=False, zaxis_title="Expression"),
        margin=dict(l=0, r=0, b=0, t=0),
        template="plotly_dark",
    )
    return fig.to_html(include_plotlyjs="cdn", full_html=False)

st.set_page_config(page_title="AI Spider Explorer", layout="wide")
st.title("AI Spider Explorer")

st.markdown("Explore spider gene expression in an interactive 3D interface.")

orthogroup = st.text_input("Orthogroup", "OG0000050")
spider_choice = st.selectbox("Choose Model Spider", ["Bridge Spider", "Clavata", "Both"]) 

if st.button("Visualize"):
    if spider_choice == "Bridge Spider":
        df = expression_data_filtered[expression_data_filtered["Orthogroup"] == orthogroup]
        if df.empty:
            st.warning("Orthogroup not found.")
        else:
            values = df.iloc[0,1:].values
            tissues = df.columns[1:]
            html = parse_and_plot(svg_files, values)
            st.components.v1.html(html, height=600)
            idx = np.argmax(values)
            st.success(f"Highest expression in {tissues[idx]} with level {values[idx]:.2f}.")
    elif spider_choice == "Clavata":
        df = expression_data_clavata_filtered[expression_data_clavata_filtered["Orthogroup"] == orthogroup]
        if df.empty:
            st.warning("Orthogroup not found.")
        else:
            values = df.iloc[0,1:].values
            tissues = df.columns[1:]
            html = parse_and_plot(svg_files, values)
            st.components.v1.html(html, height=600)
            idx = np.argmax(values)
            st.success(f"Highest expression in {tissues[idx]} with level {values[idx]:.2f}.")
    else:
        df1 = expression_data_filtered[expression_data_filtered["Orthogroup"] == orthogroup]
        df2 = expression_data_clavata_filtered[expression_data_clavata_filtered["Orthogroup"] == orthogroup]
        if df1.empty or df2.empty:
            st.warning("Orthogroup not found in both datasets.")
        else:
            values1 = df1.iloc[0,1:].values
            values2 = df2.iloc[0,1:].values
            tissues = df1.columns[1:]
            html1 = parse_and_plot(svg_files, values1)
            html2 = parse_and_plot(svg_files, values2)
            st.components.v1.html(html1, height=600)
            st.components.v1.html(html2, height=600)
