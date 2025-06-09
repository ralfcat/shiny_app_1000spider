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
import plotly.graph_objects as go
import plotly.io as pio


# paths
# Use os.path.join for portability across operating systems
blast_db_path = os.path.join("combined_genes_blastdb", "combined_genes")
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


def create_spider_visualization_3d(svg_files, expression_values):
    points = []
    for svg_file, val in zip(svg_files, expression_values):
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
    return pio.to_html(fig, full_html=False, include_plotlyjs="cdn")




app_ui = ui.page_fluid(
    ui.include_css("css/bootstrap.css"),
    ui.include_css("css/cyberpunk.css"),
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
                        choices=["Spider Visualization", "Expression Boxplot", "3D Expression"],
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
                ui.output_ui("expression_plot"),
                class_="text-center",
                style="margin-bottom: 30px;",
            ),
            ui.div(ui.output_ui("summary_text"), class_="mt-3"),
            ui.div(ui.output_ui("similar_orthogroups"), class_="mt-3"),
            ui.div(ui.output_ui("similarity_plot"), class_="mt-3"),
            ui.div(ui.output_ui("ai_report_text"), class_="mt-3"),
            ui.div(
                ui.input_text_area(
                    "chat_prompt",
                    "Ask the AI:",
                    rows=2,
                    placeholder="Type your question here...",
                ),
                class_="mt-4",
            ),
            ui.div(
                ui.input_action_button(
                    "chat_submit",
                    "Ask AI",
                    class_="btn btn-secondary btn-sm mt-2",
                ),
                class_="mb-3",
            ),
            ui.div(ui.output_ui("chat_response"), class_="mt-3"),
            ui.div(
                ui.input_action_button(
                    "report_button",
                    "Download Report",
                    class_="btn btn-primary btn-sm mt-2",
                ),
                class_="mb-3",
            ),
            ui.div(ui.output_ui("report_download"), class_="mt-3"),
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
                ui.span(" | Developed by the Spider Project Team", class_="text-white"),
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


def create_expression_boxplot_interactive(
    expression_values=None,
    tissue_names=None,
    expression_values_clavata=None,
    predicted_values=None,
):
    """Return an interactive bar plot with optional AI predictions."""
    fig = go.Figure()
    if expression_values is not None:
        fig.add_trace(
            go.Bar(
                x=tissue_names,
                y=expression_values,
                name="Bridge Spider",
                marker_color="#f3969a",
            )
        )
    if expression_values_clavata is not None:
        fig.add_trace(
            go.Bar(
                x=tissue_names,
                y=expression_values_clavata,
                name="Clavata",
                marker_color="#7cb9e8",
            )
        )
    if predicted_values is not None:
        fig.add_trace(
            go.Scatter(
                x=tissue_names,
                y=predicted_values,
                mode="lines+markers",
                name="AI Prediction",
                line=dict(color="cyan"),
            )
        )
    fig.update_layout(
        barmode="group",
        xaxis_tickangle=-45,
        yaxis_title="Expression Levels",
        template="plotly_dark",
    )
    return pio.to_html(fig, full_html=False, include_plotlyjs="cdn")

# AI helper functions

def find_similar_orthogroups(target_og, df, n=5):
    """Return the top n orthogroups most correlated with the target."""
    row = df[df["Orthogroup"] == target_og]
    if row.empty:
        return []
    target_vals = row.iloc[0, 1:].astype(float)
    numeric = df.set_index("Orthogroup").astype(float)
    corr = numeric.apply(lambda r: np.corrcoef(r.values, target_vals)[0,1], axis=1)
    corr = corr.drop(target_og, errors="ignore")
    return corr.sort_values(ascending=False).head(n).index.tolist()


def create_similarity_plot(target_vals, tissues, similar_dict):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=tissues, y=target_vals, mode="lines+markers", name="Query"))
    for og, vals in similar_dict.items():
        fig.add_trace(go.Scatter(x=tissues, y=vals, mode="lines", name=og))
    fig.update_layout(template="plotly_dark", xaxis_tickangle=-45, yaxis_title="Expression")
    return pio.to_html(fig, full_html=False, include_plotlyjs="cdn")


def generate_ai_report(expression_values, tissues):
    """Return a natural language summary using a text generation model if available."""
    try:
        from transformers import pipeline
    except Exception:
        return "Install the 'transformers' library to enable advanced summaries."

    text = ", ".join(f"{t}:{v:.2f}" for t, v in zip(tissues, expression_values))
    prompt = (
        "Summarize the following spider tissue expression levels in an engaging way: "
        + text
        + ""
    )
    try:
        gen = pipeline("text-generation", model="gpt2", framework="tf")
        out = gen(prompt, max_length=60, num_return_sequences=1)
        return out[0]["generated_text"]
    except Exception:
        return "AI summary generation failed."


def predict_expression(sequence, tissues):
    """Return pseudo AI-predicted expression levels for the sequence."""
    if not sequence:
        return None
    rng = np.random.default_rng(abs(hash(sequence)) % (2**32))
    return rng.random(len(tissues))


def generate_html_report(orthogroup, tissues, expr_values, predicted_values=None):
    """Create a basic HTML report with plots and AI summary."""
    plot_html = create_expression_boxplot_interactive(
        expression_values=expr_values,
        tissue_names=tissues,
        predicted_values=predicted_values,
    )
    summary = generate_ai_report(expr_values, tissues)
    html = f"""
    <html><head><title>Report {orthogroup}</title></head>
    <body>
    <h2>Expression Report for {orthogroup}</h2>
    {plot_html}
    <p>{summary}</p>
    </body></html>
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".html")
    with open(tmp.name, "w") as f:
        f.write(html)
    return tmp.name


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
            session_data.set(
                {
                    "best_match_id": best_match_id,
                    "e_value": e_value,
                    "orthogroup": orthogroup,
                    "species": species,
                    "sequence": fasta_sequence,
                }
            )

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
    @render.ui
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
        view = input.view_selector()
        pred_vals = None
        if "sequence" in data:
            if spider_choice == "Bridge Spider":
                pred_vals = predict_expression(data["sequence"], tissue_names)
            elif spider_choice == "Clavata":
                pred_vals = predict_expression(data["sequence"], tissue_names_clavata)
            else:
                pred_vals = predict_expression(data["sequence"], tissue_names)

        if view == "Spider Visualization":
            if spider_choice == "Bridge Spider":
                img_path = create_spider_visualization_with_animation(
                    svg_files, spider_image_path, expression_values=expr_values
                )
                return ui.img(src=img_path, width="600px", height="600px")
            elif spider_choice == "Clavata":
                img_path = create_spider_visualization_with_animation(
                    svg_files, spider_image_path, expression_values_clavata=expr_values_clavata
                )
                return ui.img(src=img_path, width="600px", height="600px")
            else:
                img_path = create_spider_visualization_with_animation(
                    svg_files, spider_image_path, expression_values=expr_values,
                    expression_values_clavata=expr_values_clavata
                )
                return ui.img(src=img_path, width="1400px", height="600px")
        elif view == "3D Expression":
            if spider_choice == "Bridge Spider":
                html = create_spider_visualization_3d(svg_files, expr_values)
                return ui.HTML(html)
            elif spider_choice == "Clavata":
                html = create_spider_visualization_3d(svg_files, expr_values_clavata)
                return ui.HTML(html)
            else:
                html1 = create_spider_visualization_3d(svg_files, expr_values)
                html2 = create_spider_visualization_3d(svg_files, expr_values_clavata)
                return ui.TagList(ui.HTML(html1), ui.HTML(html2))
        else:
            if spider_choice == "Bridge Spider":
                html = create_expression_boxplot_interactive(
                    expression_values=expr_values,
                    tissue_names=tissue_names,
                    predicted_values=pred_vals,
                )
                return ui.HTML(html)
            elif spider_choice == "Clavata":
                html = create_expression_boxplot_interactive(
                    expression_values_clavata=expr_values_clavata,
                    tissue_names=tissue_names_clavata,
                    predicted_values=pred_vals,
                )
                return ui.HTML(html)
            else:
                html = create_expression_boxplot_interactive(
                    expression_values=expr_values,
                    tissue_names=tissue_names,
                    expression_values_clavata=expr_values_clavata,
                    predicted_values=pred_vals,
                )
                return ui.HTML(html)

    @output
    @render.ui
    def summary_text():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return None

        orthogroup = data["orthogroup"]
        spider_choice = input.spider_selector()

        if spider_choice == "Bridge Spider":
            expr_df = expression_data_filtered[expression_data_filtered["Orthogroup"] == orthogroup]
            if expr_df.empty:
                return None
            values = expr_df.iloc[0, 1:].values
            tissues = expr_df.columns[1:]
            idx = np.argmax(values)
            summary = f"Highest expression in {tissues[idx]} with level {values[idx]:.2f}."
        elif spider_choice == "Clavata":
            expr_df = expression_data_clavata_filtered[expression_data_clavata_filtered["Orthogroup"] == orthogroup]
            if expr_df.empty:
                return None
            values = expr_df.iloc[0, 1:].values
            tissues = expr_df.columns[1:]
            idx = np.argmax(values)
            summary = f"Highest expression in {tissues[idx]} with level {values[idx]:.2f}."
        else:
            expr_df1 = expression_data_filtered[expression_data_filtered["Orthogroup"] == orthogroup]
            expr_df2 = expression_data_clavata_filtered[expression_data_clavata_filtered["Orthogroup"] == orthogroup]
            if expr_df1.empty or expr_df2.empty:
                return None
            val1 = expr_df1.iloc[0, 1:].values
            val2 = expr_df2.iloc[0, 1:].values
            tissues1 = expr_df1.columns[1:]
            tissues2 = expr_df2.columns[1:]
            idx1 = np.argmax(val1)
            idx2 = np.argmax(val2)
            summary = (
                f"Bridge Spider highest in {tissues1[idx1]} ({val1[idx1]:.2f}); "
                f"Clavata highest in {tissues2[idx2]} ({val2[idx2]:.2f})."
            )

        return ui.TagList(ui.h4("AI Summary", class_="mt-4"), ui.p(summary))


    @output
    @render.ui
    def similar_orthogroups():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return None
        spider_choice = input.spider_selector()
        df = expression_data_filtered if spider_choice != "Clavata" else expression_data_clavata_filtered
        sims = find_similar_orthogroups(data["orthogroup"], df)
        if not sims:
            return None
        html = "<ul>" + "".join(f"<li>{og}</li>" for og in sims) + "</ul>"
        return ui.TagList(ui.h4("Similar Orthogroups", class_="mt-4"), ui.HTML(html))

    @output
    @render.ui
    def similarity_plot():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return None
        spider_choice = input.spider_selector()
        df = expression_data_filtered if spider_choice != "Clavata" else expression_data_clavata_filtered
        row = df[df["Orthogroup"] == data["orthogroup"]]
        if row.empty:
            return None
        target_vals = row.iloc[0, 1:].astype(float).values
        tissues = row.columns[1:]
        sims = find_similar_orthogroups(data["orthogroup"], df)
        if not sims:
            return None
        sim_dict = {og: df[df["Orthogroup"] == og].iloc[0, 1:].astype(float).values for og in sims}
        html = create_similarity_plot(target_vals, tissues, sim_dict)
        return ui.HTML(html)

    @output
    @render.ui
    def ai_report_text():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return None
        spider_choice = input.spider_selector()
        df = expression_data_filtered if spider_choice != "Clavata" else expression_data_clavata_filtered
        row = df[df["Orthogroup"] == data["orthogroup"]]
        if row.empty:
            return None
        values = row.iloc[0, 1:].astype(float).values
        tissues = row.columns[1:]
        summary = generate_ai_report(values, tissues)
        return ui.TagList(ui.h4("AI Generated Report", class_="mt-4"), ui.p(summary))

    @output
    @render.ui
    @reactive.event(input.chat_submit)
    def chat_response():
        prompt = input.chat_prompt()
        if not prompt:
            return ui.p("Enter a question for the AI chat.")
        try:
            from transformers import pipeline
        except Exception:
            return ui.p("Install the 'transformers' library to enable the chat feature.")

        try:
            gen = pipeline("text-generation", model="gpt2", framework="tf")
            out = gen(prompt, max_length=50, num_return_sequences=1)
            text = out[0]["generated_text"]
        except Exception:
            text = "AI chat generation failed."

        return ui.TagList(ui.h4("AI Response", class_="mt-4"), ui.p(text))

    @output
    @render.ui
    @reactive.event(input.report_button)
    def report_download():
        data = session_data.get()
        if "orthogroup" not in data or not data["orthogroup"]:
            return ui.p("No results to summarize.")
        spider_choice = input.spider_selector()
        df = (
            expression_data_filtered
            if spider_choice != "Clavata"
            else expression_data_clavata_filtered
        )
        row = df[df["Orthogroup"] == data["orthogroup"]]
        if row.empty:
            return ui.p("Expression data not found.")
        expr_values = row.iloc[0, 1:].astype(float).values
        tissues = row.columns[1:]
        pred_vals = None
        if "sequence" in data:
            pred_vals = predict_expression(data["sequence"], tissues)
        path = generate_html_report(
            data["orthogroup"], tissues, expr_values, predicted_values=pred_vals
        )
        return ui.a("Download Report", href=path, download="report.html")

       


#run the app
app = App(app_ui, server)
