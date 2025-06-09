import subprocess
import sys
import os

def install_packages():
    """
    Installs the required dependencies for the Shiny app.
    """
    required_packages = [
        "shiny",
        "pandas",
        "matplotlib",
        "Pillow",
        "numpy",
        "plotly",
        "streamlit"
    ]
    
    print("Installing required packages...")
    for package in required_packages:
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
            print(f"Successfully installed: {package}")
        except subprocess.CalledProcessError:
            print(f"Error installing: {package}. Please check your internet connection or package name.")
    
    print("\nAll necessary packages are installed!")

def run_app(app_name):
    """
    Runs the Shiny app.
    """
    print(f"\nStarting the Spider Gene Visualization App ({app_name})...")
    try:
        subprocess.check_call([sys.executable, "-m", "shiny", "run", app_name])
    except subprocess.CalledProcessError as e:
        print(f"Error running the Shiny app: {e}. Make sure '{app_name}' is in the same directory.")

if __name__ == "__main__":
    print("Welcome to the Spider Gene Visualization App Setup!")
    print("This script will install all necessary dependencies and run the app.\n")
    
    install_packages()

    # Dynamically find the Shiny app file in the current directory
    app_filename = "shiny_app_prot.py"
    if os.path.isfile(app_filename):
        input("\nPress Enter to start the application...")
        run_app(app_filename)
    else:
        print(f"Error: '{app_filename}' not found in the current directory.")
