# Spider Gene Visualization App - Installation Instructions

## Step 1: Download the blast database
The blast database consists of fasta files of all of the filtered orthogroups. This database is REQUIRED for the Shiny app to run. The zip-file needs to be downloaded from this UPPMAX-directory:
*/proj/uppstore2019013/nobackup/private/1000spider_master_project/busco/busco_results_filtered/combined_genes_blastdb.zip*

Extract this zip folder in the directory where you have cloned this git-repository.
The folder strucutre should look like this:
![bild](https://github.com/user-attachments/assets/b6c4f132-a02a-422b-8bd3-0387f10214b6)


## Step 2: Install NCBI BLAST+

The app uses BLAST+ to perform sequence alignment searches.
### Download BLAST+
Windows:
    Download the BLAST+ installer:
    Download BLAST+ for Windows
    Run the installer and install BLAST+ on your system.

Mac/Linux:
    Use homebrew to download it. Run the command: *brew install blast*. 
    If you dont have homebrew installed, run this: 
    */bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)*


### Add BLAST+ to System PATH

To ensure the app can find BLAST+, you need to add its bin folder to your system PATH. *On macOS it did automatically when using homebrew.*

Windows:
    Open the Start Menu and search for "Environment Variables".
    Click "Edit the system environment variables".
    In the System Properties window, click "Environment Variables".
    Under "System Variables", find the Path variable and click Edit.
    Click New and paste the path to the BLAST bin folder (e.g., C:\Program Files\NCBI\blast-2.13.0+\bin).
    Click OK and restart your computer.

    
### Verify BLAST Installation

Open your terminal/command prompt and run the following command:

blastp -version  

If the installation is successful, it will show the installed BLAST version.
## Step 3: Install Dependencies

The project requires Python libraries to run. 
Run the file *setup_and_run.py*. It should install all the required packages and you will be prompted to press enter to start the application.


## Step 4: Run the Application

Run *setup_and_run.py* and press enter when prompted. If you already have the dependencies installed, you can stand in the directory of *shiny_app_prot.py* and run this command: *shiny run shiny_app_prot.py*. 
This should start the Shiny application. If it doesn't, make sure that you are standing in the correct directory which contains the app (*shiny_app_prot.pt*)

Open a browser and go to the adress prompted by the terminal. It is usually http://127.0.0.1:8000 


## How to Use the App

Enter FASTA Sequence:
    Paste a protein FASTA sequence in the box. There are random FASTA sequences located in the folder *test_fasta_sequences*

Submit:
    Click the "Submit" button to start the BLAST search.

View Results:
    BLAST Results: Displays the closest match and its associated E-value.
    Expression Data: Shows expression levels for the identified orthogroup.
    Spider Visualization: Highlights expression in specific tissues or glands.

## Troubleshooting

If the terminal outputs *INFO connection lost* after running the blast search and seeing the Spider visualization, run this command in the terminal in the directory where the app is located:
*uvicorn shiny_app_prot:app --timeout-keep-alive 240*


"BLAST Not Found" Error

    Ensure BLAST+ is installed and added to your system PATH.
    Verify installation by running:

    blastp -version  

Dependencies Not Installed

    If the required libraries fail to install, manually run:

    pip install shiny pandas matplotlib pillow numpy  

Incorrect Python Version

    Install Python 3.10 or later and ensure it is added to your system PATH.

## Credits

    Developer: Victor Englöf
    Project: Spider Gene Visualization App (2024)
    Dependencies: BLAST+, Shiny for Python, Matplotlib, Pandas, Pillow, Numpy
