import os
import shutil
from subprocess import run, DEVNULL, STDOUT

from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio import Phylo


def muscle_fasta(query_n, seq, blast_output, cov, iden, locus_dict):
    """
    From the local alignment file this function creates a new fasta file to do an alignment with muscle. Filtered by coverage and identity.
    """
    m_input_file = query_n + ".fa"
    f = open(m_input_file, "w") 
    f.write(">" + query_n + "\n" + seq[query_n] + "\n") #seq is the dictionary created with the query sequences

    input_handle = open(blast_output, "r")
    for line in input_handle:
        columns = line.split("\t")
        
        if columns[0] == query_n:
            name = columns[1].split("|")
            f.write(">" + name[0] + "|" + locus_dict[name[1]] + "\n" + columns[5])

    input_handle.close()
    f.close()
    shutil.move(m_input_file, "Results/" + m_input_file)
    return "Results/" + m_input_file


def MSA(input_f):
    """
    Create a multi sequence alignment using muscle.
    And a phylogenetic tree.
    input_f: file containing homologous sequences
    """
    output_f = input_f[:input_f.index(".")] + "_m_result.fa"
    cmd="muscle -in " + input_f + " -out " + output_f
    arg_cmd = cmd.split(" ") 
    run(arg_cmd, stdout=DEVNULL, stderr=STDOUT)

    file_tree = output_f[:output_f.index("_")] + "_tree.nw"
    cmd="muscle -maketree -in " + output_f + " -out " + file_tree + " -cluster neighborjoining"
    arg_cmd = cmd.split(" ") 
    run(arg_cmd, stdout=DEVNULL, stderr=STDOUT)

    return output_f, file_tree
    
