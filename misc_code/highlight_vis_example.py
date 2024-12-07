import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import xml.etree.ElementTree as ET
from PIL import Image
import numpy as np
from matplotlib.animation import FuncAnimation

# Load the spider image
image_path = "shiny_app/spider_image.png"  # Replace with your image path
spider_image = Image.open(image_path)

# Parse the SVG file to extract coordinates
# Parse the SVG files to extract coordinates
svg_files = [
    "shiny_app/minor_ampullate",  # Replace with your SVG file paths
    "shiny_app/Major_ampullate",
    "shiny_app/head",
]

# Mock expression data for each area
expression_values = [8.5, 4.2, 12.0]  # Example values for each area

# Normalize expression values for the colormap
norm = plt.Normalize(min(expression_values), max(expression_values))
colormap = plt.get_cmap("viridis")  # Choose a colormap, e.g., 'viridis'

# Function to parse SVG path data into x, y coordinates
def parse_path_data(path_data):
    path_commands = path_data.replace("M", "").replace("C", "").split()  # Remove 'M' and 'C' commands
    coordinates = []
    for cmd in path_commands:
        try:
            coords = cmd.split(",")
            if len(coords) == 2:  # Ensure it has both x and y values
                x, y = float(coords[0]), float(coords[1])
                coordinates.append((x, y))
        except ValueError:
            # Skip invalid entries
            continue
    return coordinates

# Parse all SVG files and store their coordinates
areas = []
for svg_file in svg_files:
    tree = ET.parse(svg_file)
    root = tree.getroot()
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    path_data = root.find(".//svg:path", namespace).attrib["d"]
    coordinates = parse_path_data(path_data)
    areas.append(coordinates)

# Identify the area with the highest expression value
max_value_index = expression_values.index(max(expression_values))

# Plot the image with heatmap-like visualization
fig, ax = plt.subplots(figsize=(10, 10))
ax.imshow(spider_image, extent=[0, spider_image.width, spider_image.height, 0])

# Add polygons with colors based on mock data
polygons = []
for coords, value in zip(areas, expression_values):
    color = colormap(norm(value))  # Map the value to a color
    polygon = patches.Polygon(
        coords, closed=True, edgecolor="black", facecolor=color, alpha=0.7, linewidth=1.5
    )
    ax.add_patch(polygon)
    polygons.append(polygon)

# Adjust plot limits
ax.set_xlim(0, spider_image.width)
ax.set_ylim(spider_image.height, 0)

# Animation function to make the max expression area blink
def animate(frame):
    alpha = 0.3 + 0.7 * np.abs(np.sin(frame * 0.1))  # Oscillates between 0.3 and 1.0
    for i, polygon in enumerate(polygons):
        if i == max_value_index:
            polygon.set_alpha(alpha)
        else:
            polygon.set_alpha(0.7)

# Create the animation
ani = FuncAnimation(fig, animate, frames=100, interval=50, repeat=True)

# Add a colorbar for reference
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation="vertical", fraction=0.02, pad=0.02)
cbar.set_label("Expression Level")

# Display the animation
plt.axis("off")
plt.show()