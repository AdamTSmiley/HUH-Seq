**HUH-Seq Analysis**

This repository contains the code used for analyzing Next-generationg sequencing data from "Sequence-Directed Covalent Protein-RNA Linkages in a Single Step Using Engineered HUH-Tags". 
The primary focus of this analysis is to identify unique k-mers and calculate their percent reduction value with respect to reference to characterize sequence specificity.

**Overview**

This script performs the following tasks:

Parses sequences from a FASTA file.
Identifies and counts unique k-mers within reads.
Calculates percent reduction values.
Exports as .CSV

**Getting Started**

This code was developed and executed in a local IDE environment. It is not specifically designed as a user-friendly tool but rather as a repository for the exact script used in the manuscript. 
However, it is made available under the MIT License, so others can take and modify the code if they wish.

**Prerequisites**

Python 3.10.13
The following Python packages:
pandas 2.1.4
biopython 1.78
numpy 1.26.3

**Usage**

To use this script:

Ensure all required packages are installed.
Run the script in your preferred Python environment.
Since this script was used for a specific analysis, it will require modification to fit other datasets or use cases.

**License**

This project is licensed under the MIT License - see the LICENSE file for details.

**Acknowledgments**

The code makes use of the Biopython library for sequence parsing.

**Data Availability**

The data files that this script analyzes are quite large (100s of GB) and as such are not readily included here.
However, these files will be made available upon request (adamtreysmiley@gmail.com).
