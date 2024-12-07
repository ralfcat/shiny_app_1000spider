from shiny import App, ui, render, reactive
import pandas as pd
import subprocess
import os
import tempfile
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from PIL import Image
import xml.etree.ElementTree as ET
import numpy as np

# Paths to data files
blast_db_path = "C:\\Users\\victo\\Documents\\1000spider\\shiny_app\\combined_genes_blastdb\\combined_genes"
orthogroup_data_path = "C:/Users/victo/Documents/1000spider/shiny_app/filtered_with_ids_no_dupes.tsv"
expression_data_path = "C:/Users/victo/Documents/1000spider/shiny_app/summed_lasc_expression.txt"
spider_image_path = "C:/Users/victo/Documents/1000spider/shiny_app/spider_image.png"
# Define the mapping between expression columns and SVG files
tissue_to_svg = {
    "Flagelliform-gland": "C:/Users/victo/Documents/1000spider/shiny_app/Flagelliform",
    "head": "C:/Users/victo/Documents/1000spider/shiny_app/head",
    "Major-ampullate-gland": "C:/Users/victo/Documents/1000spider/shiny_app/Major_ampullate",
    "Minor-ampullate-gland": "C:/Users/victo/Documents/1000spider/shiny_app/minor_ampullate",
    "Aggregate-gland": "C:/Users/victo/Documents/1000spider/shiny_app/Aggregate",
    "Aciniform-Pyriform-gland": "C:/Users/victo/Documents/1000spider/shiny_app/aciniform_pyriform",
    "Tubuliform-gland": "C:/Users/victo/Documents/1000spider/shiny_app/Tubuliform",
}

# Ensure the order of SVG files matches the columns in expression_data
svg_files = [tissue_to_svg[tissue] for tissue in tissue_to_svg.keys()]


# Load datasets
orthogroup_data = pd.read_csv(orthogroup_data_path, sep="\t", index_col=0)
expression_data = pd.read_csv(expression_data_path, sep="\t")

# Exclude specific columns
columns_to_exclude = ["wholebody", "Duct", "Sac", "Tail"]
filtered_columns = [col for col in expression_data.columns if col not in columns_to_exclude]

# Ensure "Orthogroup" is included
if "Orthogroup" not in filtered_columns:
    filtered_columns.insert(0, "Orthogroup")  # Add it at the start

expression_data_filtered = expression_data[filtered_columns]



# Helper functions
def find_orthogroup(best_match_id, orthogroup_data):
    for orthogroup, row in orthogroup_data.iterrows():
        for cell in row:
            if pd.notna(cell) and best_match_id in cell.split(", "):
                return orthogroup
    return None

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

def create_spider_visualization_with_animation(expression_values, svg_files, spider_image_path):
    spider_image = Image.open(spider_image_path)

    # Parse all SVG files
    areas = []
    for svg_file in svg_files:
        tree = ET.parse(svg_file)
        root = tree.getroot()
        namespace = {"svg": "http://www.w3.org/2000/svg"}
        path_data = root.find(".//svg:path", namespace).attrib["d"]
        coordinates = parse_path_data(path_data)
        areas.append(coordinates)

    norm = plt.Normalize(min(expression_values), max(expression_values))
    colormap = cm.viridis
    max_value_index = np.argmax(expression_values)

    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Display the image without borders
    ax.imshow(spider_image, aspect='auto')  # Removed `extent` to let it auto-scale
    ax.axis("off")  # Turn off axes to remove the border

    polygons = []
    for coords, value in zip(areas, expression_values):
        color = colormap(norm(value))
        polygon = patches.Polygon(coords, closed=True, edgecolor="black", facecolor=color, alpha=0.7, linewidth=1.5)
        ax.add_patch(polygon)
        polygons.append(polygon)

    # Animation function to make the max expression area blink and add a red border
    def animate(frame):
        alpha = 0.1 + 0.7 * np.abs(np.sin(frame * 0.1))  # Oscillates between 0.3 and 1.0
        for i, polygon in enumerate(polygons):
            if i == max_value_index:
                polygon.set_alpha(alpha)
                polygon.set_edgecolor("red")  # Set red border
                polygon.set_linewidth(3)  # Make the border thicker
            else:
                polygon.set_alpha(0.7)
                polygon.set_edgecolor("black")  # Reset other polygons' borders to black
                polygon.set_linewidth(1)  # Reset the border thickness

    # Create the animation
    ani = FuncAnimation(fig, animate, frames=100, interval=50, repeat=True)

    # Add a colorbar for reference
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", fraction=0.04, pad=0.05)
    cbar.set_label("Expression Level")

    # Save the animation as a GIF
    temp_file = tempfile.NamedTemporaryFile(suffix=".gif", delete=False)
    ani.save(temp_file.name, writer="pillow")  # Save using Pillow
    plt.close(fig)
    return temp_file.name



# Define UI
app_ui = ui.page_fluid(
    ui.h2("Spider Gene Visualization App"),
    ui.input_text_area("fasta_input", "Enter FASTA Sequence:", rows=5, placeholder="Paste your FASTA sequence here..."),
    ui.input_action_button("submit", "Submit"),
    ui.output_text_verbatim("blast_result"),
    ui.output_text_verbatim("orthogroup_result"),
    ui.output_table("expression_table"),
    ui.output_image("expression_plot"),
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
        if "error" in data:
            return data["error"]
        if "best_match_id" not in data:
            return "No BLAST result available."
        return f"Best match: {data['best_match_id']} (E-value: {data['e_value']})"

    @output
    @render.text
    def orthogroup_result():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return "No orthogroup available."
        return f"Orthogroup: {data['orthogroup']}"

    @output
    @render.table
    def expression_table():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return None
        return expression_data_filtered[expression_data_filtered["Orthogroup"] == data["orthogroup"]]

    @output
    @render.image
    def expression_plot():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return None

        orthogroup = data["orthogroup"]
        expr_data = expression_data_filtered[expression_data_filtered["Orthogroup"] == orthogroup]
        if expr_data.empty:
            return None

        expr_values = expr_data.iloc[0, 1:].values  # Exclude "Orthogroup" column
        img_path = create_spider_visualization_with_animation(expr_values, svg_files, spider_image_path)
        return {"src": img_path, "width": 800, "height": 800}


# Create Shiny App
app = App(app_ui, server)