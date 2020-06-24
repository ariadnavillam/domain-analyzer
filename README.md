# Domain analyzer 🧬💻

Python package to analyze homologous proteins. 

## Description

The aim of this program is to perform a study of a certain type of protein, introduced as the query, searching for other similar proteins, and then showing their characteristic domains.  The script follows the following steps in order to generate the results: 

1. Create a **MULTIFASTA**. 
2. Create data base.
3. Run a **BLASTP**. 
4. Do a multiple sequence alignment (**MSA**) using Muscle.
5. Draw a **phylogenetic tree** using Muscle and Phylo (Biopython module)
6. Search for **PROSITE domains** using Expasy Prosite module.
7. Create a **graphical representation** of the results using matplotlib.

One or more proteins can be introduced by the user (in the same file), so the script will create the result files for each of them.

## Requirements

Platform requirements: **Linux**.

The following programs must be installed in your computer.

- Python3.
- Biopython.
- Muscle.
- Blast. 

Files

-`prosite.dat`: must be in the same path as `main.py`. Download: ftp://ftp.expasy.org/databases/prosite/prosite.dat

Python libraries used: os, system, shutil, subprocess, re, matplotlib, numpy.

## Usage

````bash
python main.py input_query(s) subjects_folder/ [coverage] [identity]
````

| Parameter           | Description                                                  |
| ------------------- | ------------------------------------------------------------ |
| input_query(s)      | Fasta file containing one or more query sequences.           |
| subject_folder/     | Folder containing protein genbank files. These are the subject sequences. |
| Coverage (optional) | Query coverage: how much of the query sequence is covered by the subject sequence. The coverage cut-off is the minimum percentage of the query sequence that aligns with the subject (default: 50%) |
| Identity (optional) | Identity: how many characters in each sequence are identical.  The identity cut-off is the minimum percentage of identity between the query and the subject sequences (default: 25%). |

Note: the query file, the folder with the genbank files and the `prosite.dat` file must be in the same path as the main script, as well as all the python files of the package

## Output

### Folder structure

After running the script these folders and files will be created:

```bash
path/
    ├── Results/
    │   ├── [query].fa
    │   ├── [query]_m_results.fa
    │   ├── [query]_tree.nw
    │   ├── [query]_domains.txt
    │   └── Figures/
    		 ├── [query]_tree.png
    │		 ├── Domains.png
    │		 └── Blastp.png
    ├── Data/
    	├── [input_query].fa
    	├──	[subject_folder].fa
        └── DBs/
        	└── [subject_folder DB] (6 files)
    
```

#### Files in `Data` folder

- `[input_query].fa` - Input file introduced by the user. This file must contain one or more query sequences in fasta format. In order to obtain results easier to understand, its better to asign a explanatory name to the query sequences (for example: avoid using "query", "query sequence 1"..., use the name of the protein instead). The name of this file will also be the name for the blastp results file. 

- `[subject_folder].fa` - Multifasta file that contains all the genbank proteins introduced by the user. The name you assign to the folder will also be the name of the file. 

- `DBs/` - Folder that contains the files of the DB created with the script. Again, the name of the folder will also be the name of the DB too. 

#### Files in `Results` folder

- `[input_query]_results.tsv` - Blastp results with the cut-offs introduced. The e-value is by default 1.0e-5. The file contains (for each hit): [query id] [subject id] [identity] [query coverage] [e-value] [subject sequence]
- `[query].fa` - File containing the query sequence and the ones obtained in the local alignment (homologous proteins).
- `[query]_m_results.fa` - File containing the multiple sequence alignment. 
- `[query]_tree.nw` - Phylogenetic tree file obtained with muscle.
- `[query]_domains.txt` - File with the domains found for each of the sequences of the alignment. The file contains: [protein ID] [domain] [accession] [description] [pattern found].
- `Figures/` - Folder containing phylogenetic trees, a blastp results plot and a bar graphic with the most common domains.

Note: [query] is the name of the query sequence so for each of them all the previous files will be created.

### Graphical representation

All the figures are created with the matplotlib library.

- <u>Phylogenetic trees</u>. For each query introduced a tree is generated. The tree is automatically saved in `Results/Figures/`. The figure will not automatically pop-up, to visualize the tree open the png image.

- <u>Blastp plot</u>. Plot of the blastp results, representing coverage and identity. Results for homologous proteins (same query) share the same color and the e-value is represented by the size of the dot. The matches are more significant if the dots are on the top right of the graph (more coverage and identity) and the size of the dot is small (lower e-value). 

- <u>Domain analysis</u>. For each set of homologous proteins a bar graphic is created with the 5 most common domains found for that set. Each bar represents a domain and it shows how many proteins have that same domain. 

