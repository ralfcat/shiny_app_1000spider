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


#paths
blast_db_path = "combined_genes_blastdb\\combined_genes"
orthogroup_data_path = "datasets/filtered_with_ids_no_dupes.tsv"
expression_data_path = "datasets/summed_lasc_expression.txt"
clavata_expression_path = "datasets/summed_clavata_expression.txt"
spider_image_path = "visualizations/spider_image.png"
#define the mapping between expression columns and SVG files

tissue_to_svg = {
    "Flagelliform-gland": "svg/Flagelliform",
    "Major-ampullate-gland": "svg/Major_ampullate",
    "Minor-ampullate-gland": "svg/minor_ampullate",
    "Aggregate-gland": "svg/Aggregate",
    "Aciniform-Pyriform-gland": "svg/aciniform_pyriform",
    "Tubuliform-gland": "svg/Tubuliform",
}

#load the data
orthogroup_data = pd.read_csv(orthogroup_data_path, sep="\t", index_col=0)
expression_data = pd.read_csv(expression_data_path, sep="\t")
expression_data_clavata = pd.read_csv(clavata_expression_path, sep="\t")

columns_to_exclude = [
    "wholebody", "Duct", "Sac", "Tail", "Female-whole-body", "Male-whole-body", 
    "Hemolymph", "Pedipalp", "Leg", "Epidermis", "Venom-gland", "Fat-body", "head", "Ovary", "Orthogroup"
]
#match the svg files with the columns in the expression data
filtered_columns = sorted([col for col in expression_data.columns if col not in columns_to_exclude])
filtered_columns_clavata = sorted([col for col in expression_data_clavata.columns if col not in columns_to_exclude])

svg_files = [tissue_to_svg[tissue] for tissue in filtered_columns if tissue in tissue_to_svg]
svg_files_clavata = [tissue_to_svg[tissue] for tissue in filtered_columns_clavata if tissue in tissue_to_svg]

if "Orthogroup" not in filtered_columns_clavata:
    filtered_columns_clavata.insert(0, "Orthogroup")
if "Orthogroup" not in filtered_columns:
    filtered_columns.insert(0, "Orthogroup")  

expression_data_filtered = expression_data[filtered_columns]
expression_data_clavata_filtered = expression_data_clavata[filtered_columns_clavata]


#helper funcitons
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

def create_spider_visualization_with_animation(svg_files, spider_image_path, expression_values=None, expression_values_clavata=None):
    def generate_polygons(ax, expression_values, title):
        """Helper function to plot polygons with correct heatmap scaling on a given axis."""
        spider_image = Image.open(spider_image_path)

        # Parse the SVG coordinates
        areas = []
        for svg_file in svg_files:
            tree = ET.parse(svg_file)
            root = tree.getroot()
            namespace = {"svg": "http://www.w3.org/2000/svg"}
            path_data = root.find(".//svg:path", namespace).attrib["d"]
            coordinates = parse_path_data(path_data)
            areas.append(coordinates)

        # Normalize values independently for this dataset
        norm = plt.Normalize(min(expression_values), max(expression_values))
        colormap = cm.viridis
        max_value_index = np.argmax(expression_values)

        # Plot spider image
        ax.imshow(spider_image, aspect='auto')
        ax.axis("off")
        ax.set_title(title, fontsize=12, weight="bold")

        # Add polygons
        polygons = []
        for coords, value in zip(areas, expression_values):
            color = colormap(norm(value))
            polygon = patches.Polygon(coords, closed=True, edgecolor="black", facecolor=color, alpha=0.7, linewidth=1.5)
            ax.add_patch(polygon)
            polygons.append(polygon)

        # Animation function for blinking effect
        def animate(frame):
            alpha = 0.3 + 0.4 * np.abs(np.sin(frame * 0.1))
            for i, polygon in enumerate(polygons):
                if i == max_value_index:
                    polygon.set_alpha(alpha)
                    polygon.set_edgecolor("red")
                    polygon.set_linewidth(3)
                else:
                    polygon.set_alpha(0.7)
                    polygon.set_edgecolor("black")
                    polygon.set_linewidth(1)

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        return sm, animate

    # Create a figure with two subplots if both datasets are provided
    if expression_values is not None and expression_values_clavata is not None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={"wspace": 0.4})
        
        # Generate both visualizations independently
        sm1, animate1 = generate_polygons(axes[0], expression_values, "Bridge Spider")
        sm2, animate2 = generate_polygons(axes[1], expression_values_clavata, "Clavata")

        # Add colorbars for each subplot
        fig.colorbar(sm1, ax=axes[0], orientation="vertical", fraction=0.05, pad=0.04, label="Expression Level")
        fig.colorbar(sm2, ax=axes[1], orientation="vertical", fraction=0.05, pad=0.04, label="Expression Level")

        # Combine the two animations
        ani = FuncAnimation(fig, lambda frame: (animate1(frame), animate2(frame)), frames=100, interval=50, repeat=True)

    elif expression_values is not None:
        fig, ax = plt.subplots(figsize=(6, 6))
        sm1, animate1 = generate_polygons(ax, expression_values, "Bridge Spider")
        fig.colorbar(sm1, ax=ax, orientation="vertical", fraction=0.05, pad=0.04, label="Expression Level")
        ani = FuncAnimation(fig, animate1, frames=100, interval=50, repeat=True)

    elif expression_values_clavata is not None:
        fig, ax = plt.subplots(figsize=(6, 6))
        sm2, animate2 = generate_polygons(ax, expression_values_clavata, "Clavata")
        fig.colorbar(sm2, ax=ax, orientation="vertical", fraction=0.05, pad=0.04, label="Expression Level")
        ani = FuncAnimation(fig, animate2, frames=100, interval=50, repeat=True)

    else:
        raise ValueError("At least one of 'expression_values' or 'expression_values_clavata' must be provided.")

    # Save the animation as a single GIF
    temp_file = tempfile.NamedTemporaryFile(suffix=".gif", delete=False)
    ani.save(temp_file.name, writer="pillow")
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
                "Analyze and visualize spider gene expression based on FASTA input or Orthogroup number.",
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
                    "This application allows you to input either a FASTA sequence or an Orthogroup number "
                    "to analyze and visualize gene expression in spider tissues.",
                    class_="mb-2",
                ),
                ui.p(
                    "It visualizes expression levels in tissues of the model spider Larinioides sclopetarius, "
                    "with highlights for the most expressed tissues.",
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
                # Toggle Input Type
                ui.div(
                    ui.input_radio_buttons(
                        "input_type",
                        "Select Input Type:",
                        choices=["FASTA Sequence", "Orthogroup Number"],
                        selected="FASTA Sequence",
                        inline=True,
                    ),
                    style="margin-bottom: 20px; display: flex; justify-content: center; width: 100%;",
                ),
                # Dynamic Input Area
                ui.div(ui.output_ui("input_area"), style="display: flex; justify-content: center; width: 100%;",),
                ui.div(
                    ui.input_action_button(
                        "submit", "Submit", class_="btn btn-secondary btn-lg btn-block mt-3"
                    ),
                    style="display: flex; justify-content: center; width: 100%;",
                ),
                # Visualization Options
                ui.div(
                    ui.input_select(
                        "view_selector",
                        "Choose Visualization:",
                        choices=["Spider Visualization", "Expression Boxplot"],
                    ),
                    style="margin-top: 20px; display: flex; justify-content: center; width: 100%;",
                ),
                ui.div(
                    ui.input_select(
                        "spider_selector",
                        "Choose Model Spider:",
                        choices=["Bridge Spider", "Clavata", "Both"],
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
            ui.h3("Expression Data", class_="text-center mt-4", style="display: none;"),
            ui.div(ui.output_ui("expression_table"), class_="container"),
        ),
        
        ui.div(
            ui.h3("Visualization", class_="text-center mt-4", style="display: none;"),
            ui.div(
                ui.output_image("expression_plot"),
                class_="text-center",
                style="margin-bottom: 250px;",
            ),
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









def create_expression_boxplot(expression_values=None, tissue_names=None, expression_values_clavata=None):
    # Set up figure with two axes if both datasets are provided
    if expression_values_clavata is not None and expression_values is not None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # 1 row, 2 columns for side-by-side plots
        
        # Bar plot for Bridge Spider
        axes[0].bar(tissue_names, expression_values, color="#f3969a")
        axes[0].set_title("Bridge Spider")
        axes[0].set_ylabel("Expression Levels")
        axes[0].tick_params(axis="x", rotation=45)
        
        # Bar plot for Clavata
        axes[1].bar(tissue_names, expression_values_clavata, color="#7cb9e8")
        axes[1].set_title("Clavata")
        axes[1].set_ylabel("Expression Levels")
        axes[1].tick_params(axis="x", rotation=45)
    
    else:
        # If only one dataset is provided, create a single bar plot
        fig, ax = plt.subplots(figsize=(8, 6))
        if expression_values is not None:
            ax.bar(tissue_names, expression_values, color="#f3969a")
            ax.set_title("Bridge Spider")
        elif expression_values_clavata is not None:
            ax.bar(tissue_names, expression_values_clavata, color="#7cb9e8")
            ax.set_title("Clavata")
        
        ax.set_ylabel("Expression Levels")
        ax.tick_params(axis="x", rotation=45)
    
    # Save to a temporary file and return its path
    temp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    fig.savefig(temp_file.name, bbox_inches="tight")
    plt.close(fig)
    return temp_file.name


def server(input, output, session):
    session_data = reactive.Value({})  #reactive storage for session data
    @output
    @render.ui
    def input_area():
        """Conditionally render the input box based on the selected input type."""
        if input.input_type() == "FASTA Sequence":
            return ui.input_text_area(
                "fasta_input",
                "Enter FASTA Sequence:",
                rows=6,
                placeholder="Paste your FASTA sequence here...",
            )
        else:
            return ui.input_text_area(
                "orthogroup_input",
                "Enter Orthogroup Number:",
                placeholder="e.g., OG0000050",
            )

    @reactive.Effect
    @reactive.event(input.submit)
    def process_input():
        session_data.set({})  # Reset session data on new submission

        # Check input type selected by the user
        if input.input_type() == "FASTA Sequence":
            # Handle FASTA Sequence Input
            fasta_sequence = input.fasta_input()
            if not fasta_sequence:
                session_data.set({"error": "Please enter a FASTA sequence."})
                return

            # Write FASTA sequence to a temporary file
            query_file = "datasets/blast_results/query.fasta"
            with open(query_file, "w") as f:
                f.write(fasta_sequence)

            # Run BLASTP
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

            # Parse BLAST output
            with open(blast_output) as f:
                lines = f.readlines()

            if not lines:
                session_data.set({"error": "No match found from BLAST."})
                return

            # Extract best match ID and e-value
            best_match = lines[0].split("\t")
            best_match_id = best_match[1]
            e_value = best_match[10]

            # Find the corresponding Orthogroup and Species
            orthogroup, species = find_orthogroup_and_species(best_match_id, orthogroup_data)
            if orthogroup is None:
                session_data.set({"error": "Orthogroup not found."})
                return

            # Save data to session
            session_data.set({
                "best_match_id": best_match_id,
                "e_value": e_value,
                "orthogroup": orthogroup,
                "species": species,
            })

        else:
            # Handle Orthogroup Input
            orthogroup = input.orthogroup_input()
            if not orthogroup:
                session_data.set({"error": "Please enter a valid Orthogroup number."})
                return

            # Check if Orthogroup exists in the data
            if orthogroup not in expression_data_filtered["Orthogroup"].values and \
            orthogroup not in expression_data_clavata_filtered["Orthogroup"].values:
                session_data.set({"error": f"No data found for Orthogroup '{orthogroup}'."})
                return

            # Save data to session
            session_data.set({
                "orthogroup": orthogroup,
                "species": "N/A (Direct Orthogroup Input)"
            })

 
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

        if "error" in data:
            return ui.div(
                ui.strong("Error: "),
                ui.span(data["error"], style="color: red;"),
            )

        
        if "best_match_id" not in data:
            return None  

        
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
            index=False,  
            border=0,
        )

        
        return ui.TagList(
            ui.h3("BLAST Search Results", class_="mt-4"),  
            ui.HTML(table_html),
        )

    @output
    @render.ui
    def expression_table():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return None

        spider_choice = input.spider_selector()
        
        if spider_choice == "Bridge Spider":
            # Filter for Bridge Spider
            df = expression_data_filtered[expression_data_filtered["Orthogroup"] == data["orthogroup"]]
            if df.empty:
                return None

            table_html = df.to_html(
                classes="table table-primary table-striped",
                index=False,
                border=0,
            )
            return ui.TagList(
                ui.h3("Expression Data", class_="mt-4"),
                ui.h4("Bridge Spider", class_="mt-2", style="text-align: center;"),
                ui.HTML(table_html),
            )
        
        elif spider_choice == "Clavata":
            # Filter for Clavata
            df = expression_data_clavata_filtered[expression_data_clavata_filtered["Orthogroup"] == data["orthogroup"]]
            if df.empty:
                return None

            table_html = df.to_html(
                classes="table table-primary table-striped",
                index=False,
                border=0,
            )
            return ui.TagList(
                ui.h3("Expression Data", class_="mt-4"),
                ui.h4("Clavata", class_="mt-2", style="text-align: center;"),
                ui.HTML(table_html),
            )

        else:
            # Filter for both datasets
            df = expression_data_filtered[expression_data_filtered["Orthogroup"] == data["orthogroup"]]
            df_clavata = expression_data_clavata_filtered[expression_data_clavata_filtered["Orthogroup"] == data["orthogroup"]]
            
            if df.empty and df_clavata.empty:
                return None
            
            # Generate tables
            table_html = ""
            if not df.empty:
                table_html = df.to_html(
                    classes="table table-primary table-striped",
                    index=False,
                    border=0,
                )
            table_html_clavata = ""
            if not df_clavata.empty:
                table_html_clavata = df_clavata.to_html(
                    classes="table table-primary table-striped",
                    index=False,
                    border=0,
                )

            # Return combined results
            tag_list = [ui.h3("Expression Data", class_="mt-4")]

            if table_html:
                tag_list.append(ui.h4("Bridge Spider", class_="mt-2", style="text-align: center;"))
                tag_list.append(ui.HTML(table_html))
            if table_html_clavata:
                tag_list.append(ui.h4("Clavata", class_="mt-2", style="text-align: center;"))
                tag_list.append(ui.HTML(table_html_clavata))

            return ui.TagList(*tag_list)
        
    @output
    @render.ui
    def error_message():
        data = session_data.get()
        if "error" in data:
            return ui.div(
                ui.strong("Error: "),
                ui.span(data["error"], style="color: red; font-weight: bold;"),
                style="""
                    color: red; 
                    text-align: center; 
                    margin-top: 20px; 
                    display: flex; 
                    justify-content: center; 
                    align-items: center; 
                    height: 100px;  /* Optional height for vertical centering */
                """,
            )
        return None


    @output
    @render.image
    def expression_plot():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return None

        orthogroup = data["orthogroup"]
        spider_choice = input.spider_selector()

        # Check for Bridge Spider
        if spider_choice == "Bridge Spider":
            expr_data = expression_data_filtered[expression_data_filtered["Orthogroup"] == orthogroup]
            if expr_data.empty:
                return None  # Return cleanly without setting session_data
            expr_values = expr_data.iloc[0, 1:].values
            tissue_names = expr_data.columns[1:]

        # Check for Clavata
        elif spider_choice == "Clavata":
            expr_data_clavata = expression_data_clavata_filtered[
                expression_data_clavata_filtered["Orthogroup"] == orthogroup
            ]
            if expr_data_clavata.empty:
                return None  # Return cleanly without setting session_data
            expr_values_clavata = expr_data_clavata.iloc[0, 1:].values
            tissue_names_clavata = expr_data_clavata.columns[1:]

        # Check for Both
        else:
            expr_data = expression_data_filtered[expression_data_filtered["Orthogroup"] == orthogroup]
            expr_data_clavata = expression_data_clavata_filtered[
                expression_data_clavata_filtered["Orthogroup"] == orthogroup
            ]

            if expr_data.empty and expr_data_clavata.empty:
                return None
            if expr_data.empty or expr_data_clavata.empty:
                return None  # Return cleanly if one dataset is empty

            expr_values = expr_data.iloc[0, 1:].values
            expr_values_clavata = expr_data_clavata.iloc[0, 1:].values
            tissue_names = expr_data.columns[1:]

        # Generate plots
        if input.view_selector() == "Spider Visualization":
            if spider_choice == "Bridge Spider":
                img_path = create_spider_visualization_with_animation(
                    svg_files, spider_image_path, expression_values=expr_values
                )
                return {"src": img_path, "width": "600px", "height": "600px"}
            elif spider_choice == "Clavata":
                img_path = create_spider_visualization_with_animation(
                    svg_files, spider_image_path, expression_values_clavata=expr_values_clavata
                )
                return {"src": img_path, "width": "600px", "height": "600px"}
            else:
                img_path = create_spider_visualization_with_animation(
                    svg_files, spider_image_path, expression_values=expr_values,
                    expression_values_clavata=expr_values_clavata
                )
                return {"src": img_path, "width": "1400px", "height": "600px"}
        else:
            if spider_choice == "Bridge Spider":
                img_path = create_expression_boxplot(expression_values=expr_values, tissue_names=tissue_names)
                return {"src": img_path, "width": "600px", "height": "600px"}
            elif spider_choice == "Clavata":
                img_path = create_expression_boxplot(expression_values_clavata=expr_values_clavata, tissue_names=tissue_names_clavata)
                return {"src": img_path, "width": "600px", "height": "600px"}
            else:
                img_path = create_expression_boxplot(
                    expression_values=expr_values, tissue_names=tissue_names,
                    expression_values_clavata=expr_values_clavata
                )
                return {"src": img_path, "width": "1100px", "height": "600px"}


       


#run the app
app = App(app_ui, server)
