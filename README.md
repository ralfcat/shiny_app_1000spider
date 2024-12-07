Spider Gene Visualization App - Installation Instructions
Step 1: Clone or Download the Project

To get the project files, you can either clone the repository using Git or download the files manually.

    Using Git:
    Open your terminal or command prompt and run the following commands:

    git clone https://github.com/yourusername/spider-gene-visualization.git  
    cd spider-gene-visualization  

    Manual Download:
    Go to the repository on GitHub, click on "Download ZIP", and extract the contents to your desired location.

Step 2: Install Dependencies

The project requires Python libraries to run. To install these, follow these steps:

    Ensure you have Python 3.10 or higher installed.
        Download Python
        During installation, check "Add Python to PATH".

    Install the required Python libraries:

    Open your terminal or command prompt in the project folder and run:

    pip install -r requirements.txt  

    If pip is not installed, follow the official guide here:
    Install pip for Python.

Step 3: Install NCBI BLAST+

The app uses BLAST+ to perform sequence alignment searches.
1. Download BLAST+

    Windows:
        Download the BLAST+ installer:
        Download BLAST+ for Windows
        Run the installer and install BLAST+ on your system.

    Mac/Linux:
        Download the BLAST+ archive:
        Download BLAST+
        Extract the downloaded archive to your desired location.

2. Add BLAST+ to System PATH

To ensure the app can find BLAST+, you need to add its bin folder to your system PATH.

    Windows:
        Open the Start Menu and search for "Environment Variables".
        Click "Edit the system environment variables".
        In the System Properties window, click "Environment Variables".
        Under "System Variables", find the Path variable and click Edit.
        Click New and paste the path to the BLAST bin folder (e.g., C:\Program Files\NCBI\blast-2.13.0+\bin).
        Click OK and restart your computer.

    Mac/Linux:
    Add the following line to your shell profile (~/.bashrc or ~/.zshrc):

export PATH=/path/to/blast/bin:$PATH  

Replace /path/to/blast/bin with the actual path where BLAST+ was extracted.

Apply changes by running:

    source ~/.bashrc  

3. Verify BLAST Installation

Open your terminal/command prompt and run the following command:

blastp -version  

If the installation is successful, it will show the installed BLAST version.
Step 4: Run the Application

    Open a terminal or command prompt in the project folder.

    Run the application:

python app.py  

If you encounter version issues, use:

python3 app.py  

Once the app starts successfully, it will provide a local URL:

    http://127.0.0.1:8000  

    Open this URL in your browser to access the app.

How to Use the App

    Enter FASTA Sequence:
        Paste a valid nucleotide or protein FASTA sequence into the input box.

    Select BLAST Type:
        Choose between:
            BLASTP: For protein sequences.
            BLASTX: For nucleotide sequences translated to protein.

    Submit:
        Click the "Submit" button to start the BLAST search.

    View Results:
        BLAST Results: Displays the closest match and its associated E-value.
        Expression Data: Shows expression levels for the identified orthogroup.
        Spider Visualization: Highlights expression in specific tissues or glands.

Troubleshooting
"BLAST Not Found" Error

    Ensure BLAST+ is installed and added to your system PATH.
    Verify installation by running:

    blastp -version  

Dependencies Not Installed

    If the required libraries fail to install, manually run:

    pip install shiny pandas matplotlib pillow numpy  

Incorrect Python Version

    Install Python 3.10 or later and ensure it is added to your system PATH.

Credits

    Developer: Victor Englöf
    Project: Spider Gene Visualization App (2024)
    Dependencies: BLAST+, Shiny for Python, Matplotlib, Pandas, Pillow, Numpy
