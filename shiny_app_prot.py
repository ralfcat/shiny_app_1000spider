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
from pathlib import Path


# Paths to data files
blast_db_path = "combined_genes_blastdb\\combined_genes"
orthogroup_data_path = "datasets/filtered_with_ids_no_dupes.tsv"
expression_data_path = "datasets/summed_lasc_expression.txt"
spider_image_path = "visualizations/spider_image.png"
# Define the mapping between expression columns and SVG files
tissue_to_svg = {
    "Flagelliform-gland": "svg/Flagelliform",
    "head": "svg/head",
    "Major-ampullate-gland": "svg/Major_ampullate",
    "Minor-ampullate-gland": "svg/minor_ampullate",
    "Aggregate-gland": "svg/Aggregate",
    "Aciniform-Pyriform-gland": "svg/aciniform_pyriform",
    "Tubuliform-gland": "svg/Tubuliform",
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
def find_orthogroup_and_species(best_match_id, orthogroup_data):
    for orthogroup, row in orthogroup_data.iterrows():
        for column in orthogroup_data.columns[1:]:  # Skip the first column (Orthogroup)
            cell = row[column]
            if pd.notna(cell) and best_match_id in str(cell).split(", "):
                species = column.split(".")[0]  # Extract species name from column (before '.fasta')
                return orthogroup, species
    return None, "Unknown"


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

    fig, ax = plt.subplots(figsize=(8, 8))
    
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



app_ui = ui.page_fluid(
    ui.include_css("css/bootstrap.css"),
    ui.busy_indicators.options(
        spinner_type=Path("www/spinner.svg"),
        spinner_color="black",
        spinner_size="100px",
        spinner_delay="0s",
    ),
    # Flexbox wrapper for sticky footer and main content
    ui.Tag(
        "div",
        {"class": "d-flex flex-column min-vh-100"},
        # Header Section
        ui.div(
            ui.h2("Spider Gene Visualization App", class_="text-center mt-4"),
            ui.p(
                "Analyze and visualize spider gene expression based on FASTA input.",
                class_="text-center text-muted",
            ),
        ),
        # Main Content Wrapper
        ui.div(
            {
                "class": "d-flex justify-content-center align-items-start flex-grow-1",
                "style": "position: relative; margin-top: 30px;",
            },
            # About Box Positioned to the Left
            ui.Tag(
                "div",
                {
                    "class": "card p-3 shadow-lg",
                    "style": (
                        "width: 20%; max-width: 300px; border: 3px solid #86d3c9; "
                        "position: absolute; left: 10%; top: 0;"
                    ),
                },
                ui.h4("About", class_="text-center mb-3"),
                ui.p(
                    "This application allows you to input a FASTA sequence to identify the corresponding orthogroup "
                    "based on a dataset containing 9000+ orthogroups from 292 spider species.",
                    class_="mb-2",
                ),
                ui.p(
                    "It visualizes the orthogroups' expression levels in tissues of the model spider Larinioides sclopetarius, "
                    "with the most expressed tissue highlighted.",
                    class_="mb-2",
                ),
            ),
            # Input Box Centered Horizontally
            ui.Tag(
                "div",
                {
                    "class": "card p-4 shadow-lg",
                    "style": (
                        "width: 40%; max-width: 600px; border: 3px solid #86d3c9; margin: 0 auto;"
                    ),
                },
                ui.h4("Input Your Data", class_="text-center mb-4"),
                ui.div(
                    ui.input_text_area(
                        "fasta_input",
                        "Enter FASTA Sequence:",
                        rows=6,
                        placeholder="Paste your FASTA sequence here...",
                    ),
                    style="display: flex; justify-content: center; width: 100%;",
                ),
                ui.div(
                    ui.input_action_button(
                        "submit", "Submit", class_="btn btn-secondary btn-lg btn-block mt-3"
                    ),
                    style="display: flex; justify-content: center; width: 100%;",
                ),
                ui.div(
                    ui.input_select(
                        "view_selector",
                        "",
                        choices=["Spider Visualization", "Expression Boxplot"],
                    ),
                    style="margin-top: 20px; display: flex; justify-content: center; width: 100%;",
                ),
            ),
        ),
                # Results Section: BLAST Results, Expression Data, and Visualization
        ui.div(
            ui.h3("BLAST Results", class_="text-center mt-4", style="display: none;"),
            ui.div(ui.output_ui("blast_result"), class_="container"),
        ),
        ui.div(
            ui.h3(
                "Expression Data",
                class_="text-center mt-4",
                style="display: none;",
            ),  # Initially hidden
            ui.div(ui.output_ui("expression_table"), class_="container"),
        ),
        ui.div(
            ui.h3(
                "Visualization",
                class_="text-center mt-4",
                style="display: none;",
            ),  # Initially hidden
            ui.div(ui.output_image("expression_plot"), class_="text-center", style="margin-bottom: 250px;",),
        ),
        # Footer Section
        ui.Tag(
            "footer",
            {
                "class": "text-center mt-auto py-2",
                "style": (
                    "position: fixed; bottom: 0; left: 0; width: 100%; "
                    "height: 40px; color: white; background-color: #78c2ad; "
                    "padding-top: 10px; z-index: 1000; border: none; margin: 0; box-shadow: none;"
                ),
            },
            ui.div(
                ui.span("1000 Spider Project 2024", class_="text-white"),
                ui.span(" | Developed by Victor Englöf", class_="text-white"),
                class_="small",
            ),
        ),
    ),
)








def create_expression_boxplot(expression_values, tissue_names):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(tissue_names, expression_values, color="#f3969a")
    ax.set_ylabel("Expression Levels")
    ax.set_title("Expression Levels by Tissue")
    ax.tick_params(axis="x", rotation=45)

    temp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    fig.savefig(temp_file.name, bbox_inches="tight")
    plt.close(fig)
    return temp_file.name


def server(input, output, session):
    session_data = reactive.Value({})  # Reactive storage for session data

    @reactive.Effect
    @reactive.event(input.submit)
    def run_blast():
        fasta_sequence = input.fasta_input()
        if not fasta_sequence:
            session_data.set({"error": "Please enter a FASTA sequence."})
            return

        query_file = "datasets/blast_results/query.fasta"
        with open(query_file, "w") as f:
            f.write(fasta_sequence)

        blast_output = "datasets/blast_results/blast_results.txt"
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
        session_data.set(
            {
                "best_match_id": best_match[1],
                "e_value": best_match[10],
                "species": "Unknown",
                "orthogroup": None,
            }
        )

    @reactive.Effect
    def find_and_visualize_orthogroup():
        data = session_data.get()
        if "best_match_id" not in data or data.get("orthogroup") is not None:
            return

        orthogroup, species = find_orthogroup_and_species(data["best_match_id"], orthogroup_data)
        if orthogroup is None:
            session_data.set({"error": "Orthogroup not found."})
            return

        data["orthogroup"] = orthogroup
        data["species"] = species
        session_data.set(data)

    @output
    @render.ui
    def blast_result():
        data = session_data.get()

        # Check for errors and display them
        if "error" in data:
            return ui.div(
                ui.strong("Error: "),
                ui.span(data["error"], style="color: red;"),
            )

        # If no BLAST result is available yet
        if "best_match_id" not in data:
            return None  # No UI rendered when no results are available

        # Display the best match results in a styled table
        species = data.get("species", "Unknown")
        df = pd.DataFrame(
            [
                {
                    "Best Match ID": data["best_match_id"],
                    "E-Value": data["e_value"],
                    "Species": species,
                }
            ]
        )

        table_html = df.to_html(
            classes="table table-primary table-striped",
            index=False,  # Hide index column
            border=0,
        )

        # Dynamically display the header with the table
        return ui.TagList(
            ui.h3("BLAST Search Results", class_="mt-4"),  # Add header dynamically
            ui.HTML(table_html),
        )

    @output
    @render.ui
    def expression_table():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return None

        df = expression_data_filtered[expression_data_filtered["Orthogroup"] == data["orthogroup"]]
        if df.empty:
            return None

        # Render the table with a custom "table-primary" class and dynamically display the header
        table_html = df.to_html(
            classes="table table-primary table-striped",
            index=False,  # Hide the index column
            border=0,
        )
        return ui.TagList(
            ui.h3("Expression Data", class_="mt-4"),  # Add header dynamically
            ui.HTML(table_html),
        )

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
        tissue_names = expr_data.columns[1:]  # Exclude "Orthogroup"

        if input.view_selector() == "Spider Visualization":
            img_path = create_spider_visualization_with_animation(expr_values, svg_files, spider_image_path)
        else:
            img_path = create_expression_boxplot(expr_values, tissue_names)

        # Return the expected dictionary with src, width, and height
        return {
            "src": img_path,
            "width": "600px",
            "height": "600px",
        }


# Create Shiny App
app = App(app_ui, server)
