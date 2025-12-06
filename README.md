# ComTarget
ComTarget is a free, cross-platform, scientific software for  Chemical Target Prediction with Combinatorial Modeling. It is a command line oriented program. It can be executed in a shell (UNIX/Linux systems) or at a Windows or DOS command prompt. 
# Software Requirement
ComTarget requires the assistance of Obabel (Open Babel version 3.0.0 or higher) for format conversion and conformation search operations. Before using ComTarget, obabel must be installed in advance with environment variables properly configured, and it should be ensured that the confab module of obabel can be called normally. ComTarget also requires the assistance of AutoDock Vina to perform reverse docking operations. Please ensure that Vina has been installed and is callable.
# Database Download
ComTarget requires a database of protein structures. Please download the protein database from Figshare:

Download URL: https://figshare.com/articles/dataset/Reverse_Target_Fishing_Protein_PDBQT_Database/30812579

DOI: 10.6084/m9.figshare.30812579

After downloading, extract the files to create the following directory structure (/info_data/receptor_data/):
```markdown
- ComTarget/
  - ComTarget.c
  - LICENSE
  - README.md
  - info_data/
    - box.data
    - discriptor.txt
    - ele_data.txt
    - score_data.txt
    - targets.tsv
    - receptor_data/          # ← Extract downloaded files here
      - 1ABC.pdbqt
      - 2DEF.pdbqt
      - ... (and other PDBQT files)
```
## **Installation**
The ComTarget program is written in C and has been tested on Linux and Windows platforms. It is provided as source codes. To install the command line version of ComTarget please follow the instructions below.
1. Install GNU C compiler if you haven't already. On Ubuntu you may use sudo apt-get install gcc, for other Linux users you may need to use a different method.
2. Downloading the source code and info_data
3. compiling source code

```
gcc ComTarget.c -o ComTarget -lm -O3         
```
  or

```
gcc ComTarget.c -o ComTarget -lm 
```
If you are a Windows user, you can compile using gcc on the command prompt. Alternatively, install the free, portable, fast and simple C/C++ IDE Dev-C++, then Using Dev-C++, open the source code, click the compile option to compile and obtain the binary executable file.


## **How to Use**
Add ComTarget as an environment variable.
Open Terminal (Mac & Linux) or PowerShell (Windows).
The general synopsis for using ComTarget is:
```
ComTarget inputFile 
```
## **Supported Input Formats:**
MOL2 (.mol2) - Recommended format
SDF (.sdf) - Structure Data File
PDB (.pdb) - Protein Data Bank format
XYZ (.xyz) - Coordinate format
## **Workflow:**
1. Input your molecule file
2. Program will prompt for conformation generation method (Monte Carlo or Confab)
3. Select Boltzmann distribution threshold (recommended: 0.5-5.0%)
4. Choose whether to perform molecular docking (for secondary screening)
5. Select number of targets for docking (based on similarity ranking)
6. Results will be saved as:
    similarity_results_[input_name].txt
    docking_results_[input_name].csv
    docking_details_[input_name].txt
## **Output Interpretation:**
The CSV file contains prediction results including:

  Target PDB ID
  Docking scores (best, average, worst)
  Delta scores (difference from original ligand)
  Protein and gene names
  Uniprot ID
  Prediction category (HIGH/GOOD/MEDIUM/LOW/VERY LOW)
## **Citation**
If you use ComTarget in your research, please cite:

Database:
Li, Y. (2024). Reverse Target Fishing Protein PDBQT Database. figshare. https://doi.org/10.6084/m9.figshare.30812579

Software:
Li, Y. (2024). ComTarget: Chemical Target Prediction with Combinatorial Modeling. GitHub repository.

## **License**
ComTarget is distributed under the MIT License. See LICENSE file for details.

## **Contact**
For more information or support, please contact: drugyuzhu@163.com
