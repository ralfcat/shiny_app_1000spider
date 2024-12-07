from shiny import App, ui, render, reactive
import pandas as pd
import subprocess
import os
import tempfile
from matplotlib import pyplot as plt, patches, cm
from PIL import Image
import xml.etree.ElementTree as ET
import plotly.express as px
import plotly.io as pio

# Paths to data files
blast_db_path = "C:\\Users\\victo\\Documents\\1000spider\\shiny_app\\combined_genes_blastdb\\combined_genes"
orthogroup_data_path = "C:/Users/victo/Documents/1000spider/shiny_app/filtered_with_ids_no_dupes.tsv"
expression_data_path = "C:/Users/victo/Documents/1000spider/shiny_app/summed_lasc_expression.txt"
spider_image_path = "C:/Users/victo/Documents/1000spider/shiny_app/spider_image.png"  # Replace with your image path
svg_files = [
    "Flagelliform",  # Replace with your SVG file paths
    "head",
    "Major_ampullate",
    "minor_ampullate",
    "Aggregate",
    "aciniform_pyriform",
    "Tubuliform",
]

# Load datasets
orthogroup_data = pd.read_csv(orthogroup_data_path, sep="\t", index_col=0)
expression_data = pd.read_csv(expression_data_path, sep="\t")

# Helper functions
def find_orthogroup(best_match_id, orthogroup_data):
    for orthogroup, row in orthogroup_data.iterrows():
        for cell in row:
            if pd.notna(cell) and best_match_id in cell.split(", "):
                return orthogroup
    return None

def create_spider_visualization(expression_data, orthogroup, svg_files):
    # Normalize mock expression values
    norm = plt.Normalize(min(expression_data), max(expression_data))
    colormap = cm.viridis

    # Parse SVG files for coordinates
    def parse_path_data(path_data):
        path_commands = path_data.replace("M", "").replace("C", "").split()
        coordinates = []
        for cmd in path_commands:
            try:
                coords = cmd.split(",")
                if len(coords) == 2:
                    x, y = float(coords[0]), float(coords[1])
                    coordinates.append((x, y))
            except ValueError:
                continue
        return coordinates

    areas = []
    for svg_file in svg_files:
        tree = ET.parse(svg_file)
        root = tree.getroot()
        path_data = root.find(".//{http://www.w3.org/2000/svg}path").attrib["d"]
        coordinates = parse_path_data(path_data)
        areas.append(coordinates)

    # Plot the spider visualization
    fig, ax = plt.subplots(figsize=(10, 10))
    spider_image = Image.open(spider_image_path)
    ax.imshow(spider_image, extent=[0, spider_image.width, spider_image.height, 0])
    
    for coords, value in zip(areas, expression_data):
        color = colormap(norm(value))
        polygon = patches.Polygon(coords, closed=True, facecolor=color, edgecolor="black", alpha=0.7)
        ax.add_patch(polygon)

    ax.axis("off")
    temp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    fig.savefig(temp_file.name, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    return temp_file.name

# Define UI
app_ui = ui.page_fluid(
    ui.h2("Spider Gene Visualization App"),
    ui.input_text_area("fasta_input", "Enter FASTA Sequence:", rows=5, placeholder="Paste your FASTA sequence here..."),
    ui.input_action_button("submit", "Submit"),
    ui.input_select("visualization_type", "Select Visualization Type:", choices=["Bar Chart", "Spider Visualization"]),
    ui.output_text_verbatim("blast_result"),
    ui.output_text_verbatim("orthogroup_result"),
    ui.output_table("expression_table"),
    ui.output_ui("dynamic_visualization"),
)

# Server logic
def server(input, output, session):
    session_data = reactive.Value({})  # Reactive storage for session data

    @reactive.Effect
    @reactive.event(input.submit)
    def run_blast():
        fasta_sequence = input.fasta_input()
        if not fasta_sequence:
            session_data.set({"error": "Please enter a FASTA sequence."})
            return

        query_file = "query.fasta"
        with open(query_file, "w") as f:
            f.write(fasta_sequence)

        blast_output = "blast_results.txt"
        try:
            subprocess.run(
                [
                    "blastp",
                    "-query", query_file,
                    "-db", blast_db_path,
                    "-out", blast_output,
                    "-outfmt", "6",
                    "-max_target_seqs", "1",
                ],
                check=True,
            )
        except subprocess.CalledProcessError:
            session_data.set({"error": "BLAST failed. Please check your input."})
            return

        with open(blast_output) as f:
            lines = f.readlines()

        if not lines:
            session_data.set({"error": "No match found."})
            return

        best_match = lines[0].split("\t")
        session_data.set({"best_match_id": best_match[1], "e_value": best_match[10], "orthogroup": None})

    @reactive.Effect
    def find_and_visualize_orthogroup():
        data = session_data.get()
        if "best_match_id" not in data or data["orthogroup"] is not None:
            return

        orthogroup = find_orthogroup(data["best_match_id"], orthogroup_data)
        if orthogroup is None:
            session_data.set({"error": "Orthogroup not found."})
            return

        data["orthogroup"] = orthogroup
        session_data.set(data)

    @output
    @render.text
    def blast_result():
        data = session_data.get()
        if not data or "error" in data:
            return data.get("error", "No BLAST result available.")
        if "best_match_id" not in data:
            return "No BLAST result available."
        return f"Best match: {data['best_match_id']} (E-value: {data['e_value']})"

    @output
    @render.text
    def orthogroup_result():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return "No orthogroup available."
        return f"Gene {data['best_match_id']} belongs to orthogroup {data['orthogroup']}"

    @output
    @render.table
    def expression_table():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return None
        return expression_data[expression_data["Orthogroup"] == data["orthogroup"]]

    @output
    @render.ui
    def dynamic_visualization():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return ui.tags.div("No visualization available. Please ensure the BLAST process completes successfully.")

        orthogroup = data["orthogroup"]
        expr_data = expression_data[expression_data["Orthogroup"] == orthogroup]

        if expr_data.empty:
            return ui.tags.div("No expression data found for the selected orthogroup.")

        if input.visualization_type() == "Bar Chart":
            # Prepare bar chart visualization using Plotly
            melted_data = expr_data.melt(
                id_vars=["Orthogroup"],
                var_name="Tissue",
                value_name="Expression Level",
            )

            # Create a Plotly bar chart
            fig = px.bar(
                melted_data,
                x="Tissue",
                y="Expression Level",
                title=f"Tissue Expression Levels for {orthogroup}",
            )

            # Save to a temporary file using kaleido
            temp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
            pio.write_image(fig, temp_file.name, format="png")

            return ui.tags.img(src=temp_file.name, style="width: 800px; height: 600px;")

        elif input.visualization_type() == "Spider Visualization":
            # Prepare spider visualization
            expr_values = expr_data.iloc[0, 1:].values
            img_path = create_spider_visualization(expr_values, orthogroup, svg_files)
            return ui.tags.img(src=img_path, style="width: 800px; height: 600px;")



# Create Shiny App
app = App(app_ui, server)



