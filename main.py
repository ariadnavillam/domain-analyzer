#!/usr/bin/env python
# -*- coding: utf-8 *-*
import os
import sys
import shutil
from subprocess import DEVNULL, STDOUT, run

from Bio import SeqIO
from Bio.ExPASy import Prosite,Prodoc
from Bio.Align.Applications import MuscleCommandline

import blast as bs 
import muscle as mc 
import prosite as ps 
import graphics as gf

def help_msg():
    print("\nUSAGE: main.py [query(s)_file] [genbank_folder/] [cov(>=50)] [id(>=25)]")
    print("Use -h or --help to print descriptions of the arguments.\n")

def help_ext():
    print("\nUSAGE \n\t main.py [query(s)_file] [genbank_folder] [cov(>=50)] [id(>=25)]\n")
    print("DESCRIPTION")
    print("\tThis script runs a blastp search against genbank proteins to obtain")
    print("\thomologous proteins of a query sequence. Then creates a MSA, ")
    print("\ta phylogetic tree and extracts domains of the set of proteins. The files")
    print("\tgenerated are saved in Results/.\n")
    print("ARGUMENTS")
    print("\t-query(s)_file \n\t   fasta file with one or more proteins to analise.")
    print("\t-genbank_folder/ \n\t   folder with genbank files to create a data base to search for ")
    print("\t   homologous proteins.\n")
    print("· Optional arguments")
    print("\t-cov \n\t   coverage cutoff. Minimum percentage of query coverage for a hit.")
    print("\t-id \n\t   identity cutoff. Minimun percentage of identical characters in the")
    print("\t   two sequences.")
    print("\t*When not provided the default values are: cov = 50, id = 25\n")
    print("See ReadMe.md for further information.\n")
    sys.exit()

def check_file(check_file, form_file):
    """
    Check if a file is genbank or fasta.
    If the file format is not correct the script ends.
    """
    try: 
        f = open(check_file, "r")
        a = list(SeqIO.parse(f, form_file))
        b = a[0]

    except:
        print("File "+ check_file + " is not in " + form_file + " format.")
        help_msg()
        sys.exit()

    f.close()

def arguments():
    """
    Number of arguments must be 2 or 4. 
    Value of cov and id is default if not introduced by the user.
    """
    if len(sys.argv) == 3:
        cov = 50
        iden = 25

    elif len(sys.argv) == 5:
        if sys.argv[3].isnumeric() and sys.argv[4].isnumeric():
            cov = sys.argv[3]
            iden = sys.argv[4]
        else:
            print("Error. Coverage cut-off and identity cut-off must be numbers.")
            help_msg()
            sys.exit()

    else:
        print("Error. Incorrect number of arguments.")
        help_msg()
        sys.exit()  

    return (cov, iden) 

if sys.argv[1] in ('-h', '--help'):
    help_ext()

cov, iden = arguments() #default values or the ones set by the user

query_file = sys.argv[1]
check_file(query_file, "fasta") 

subject_folder = sys.argv[2]
if os.path.isdir(subject_folder) == False: #check if the argument introduced is a folder
    print("Error. No folder found.")
    help_msg()

#name of the folder must end with "/". If not, a "/" is added 
if subject_folder[-1] != "/":
    subject_folder = subject_folder + "/"

subject_fasta = subject_folder[:-1] + ".fa" #name for the multifasta is the same as the folder

#dict to save all locus names and the organism
locus_org=dict()

if os.path.exists(subject_fasta):
    os.remove(subject_fasta)

#create directories. if folder exists a new one is not created.
os.makedirs("Data/", exist_ok=True)
os.makedirs("Data/DBs/", exist_ok=True)
os.makedirs("Results/", exist_ok=True)

for file_p in os.listdir(subject_folder): #check if files in the folder are genbank files
    check_file(subject_folder + file_p, "genbank")
    locus_dict=bs.convert_seq(subject_folder + file_p, subject_fasta, locus_org)

print("---------------OUTPUT--------------")
print("The arguments introduced are correct. \n")
print("Files generated: \nMultifasta file: Data/" + subject_fasta)

dbname = bs.create_db(subject_fasta) #create database and return name

shutil.copy(query_file, "Data/"+ query_file)
shutil.move(subject_fasta, "Data/" + subject_fasta)

print("Database: Data/DBs/" + dbname + ".\n")

blast_output, range_cov, range_iden, hits = bs.blast_alignment("Data/" + query_file, "Data/DBs/" + dbname, query_file[:-3] + "_result.out", cov, iden)

if os.stat(blast_output).st_size == 0:
    print("No hits found for these parameters.")
    help_msg()
    sys.exit()

print("···BLASTP RESULTS···")
print("Blastp output file: " + blast_output)
print("Identity cut-off: " + str(iden))
print("Coverage cut-off: " + str(cov))
print("Query file: Data/" + query_file)
print("Subject multifasta file: Data/" + subject_fasta)
print("Identity range: [ " + range_iden + " ]")
print("Coverage range: [ " + range_cov + " ]")
print("Hits (number of sequences aligned): " + str(hits) + "\n")

query_sequences = bs.query_dict(query_file) #dictionary containing name and sequence for each query in the file

pattern_dict = ps.c_patterns() #dictionary containing patterns from prosite

print("···MSA RESULTS···")
print("Query sequences introduced: " + str(len(query_sequences)) )
pro_domains = list() #list containing lists of prosite domains found for each query
for key in query_sequences:
    print(key)
    m_f = mc.muscle_fasta(key, query_sequences, blast_output, cov, iden, locus_dict)
    align_file, tree = mc.MSA(m_f)

    pro_domains.append(ps.find_patterns(m_f, pattern_dict)) 

print("\nFiles created for each query: ")
print("\tMultifasta file: Results/[query].fa")
print("\tAlignment file: Results/[query]_m_results.fa")
print("\tTree file: Results/[query]_tree.nw\n")

print("···DOMAINS RESULTS···")
print("Prosite domains found for the sequences align with query sequence.")
print("File: Results/[query]_domains.txt")

query_names = [key for key in query_sequences]   

#graficos
os.makedirs("Results/Figures/", exist_ok=True)
print("\n···FIGURES···")
print("*Close the windows with graphics in order to end the script.*")
print("All the images are automatically saved.")
tr_im  = gf.graph_tree(query_names)
bl_im  = gf.graph_blastp(blast_output, query_names)
dom_im = gf.graph_domains(pro_domains, query_names)
print("\nPhylogenetic trees: " + tr_im)
print("Figure showing blastp results: " + bl_im)
print("Figure showing most common domains per query sequence : " + dom_im)





