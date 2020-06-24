import os
import sys
import shutil
from subprocess import DEVNULL, STDOUT, run
from Bio import Seq, SeqIO
from Bio.ExPASy import Prosite,Prodoc
from Bio.Blast.Applications import NcbiblastpCommandline as blastp

def convert_seq(input_file, output_file, locus_dict, input_form="genbank", output_form="fasta"):
    """
    Generate a fasta file from a genbank file. 
    """ 
    output_handle=open(output_file, "a") 

    with open(input_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            org_name = record.annotations["organism"]
            org_name = org_name[:org_name.index(" PCC")]
            org_name = org_name.replace(" ", "_")
            locus_dict[record.id] = org_name #dictionary with the locus id and the corresponding organism

            for feature in record.features:
                if feature.type == "CDS":

                    try:
                        output_handle.write(">" + feature.qualifiers["locus_tag"][0] + "|" + record.id + "\n" + feature.qualifiers["translation"][0] + "\n")
                        
                    except:
                        pass
                        
    output_handle.close()
    input_handle.close()

    return locus_dict

def create_db(multifasta):
    """
    Create a data base for the multifasta file. The DB is saved in /Data/DBs/
    """
    if os.path.exists("Data/DBs/"): 
        try: 
            dbname = multifasta[:-3]
            cmd = "makeblastdb -in " + multifasta + " -dbtype prot -parse_seqids -out Data/DBs/" + dbname
            arg_cmd = cmd.split(" ") #pasamos el comando en forma de lista para que no de error
            run(arg_cmd, stdout = DEVNULL, stderr = STDOUT)

        except:
            print("Error. Could not create data base.")
            sys.exit()

    else:
        print("Data/DBs/ not found.") #this error shouldn occur if main.py is executed
        sys.exit()

    return dbname

def blast_alignment(query, db, output, cov, iden):
    """
    Blastp alignment with query sequence and data base
    """
    try:
        cmd = 'blastp -query ' + query + ' -db ' + db + ' -evalue 1.0e-5 -outfmt "6 qseqid sseqid pident qcovs evalue" -out ' + output 
        os.system(cmd)

    except:
        print("Error in blastp alignment.")
        sys.exit()

    #filter file by coverage and iden cut-offs
    input_handle = open(output, "r")
    output_handle = open("filtered", "w")
    l_cov = list()
    l_iden = list()
    for line in input_handle:
        columns = line.split("\t")
        if float(columns[2]) >= float(iden) and float(columns[3]) >= float(cov):
            l_cov.append(columns[3])
            l_iden.append(columns[2])
            output_handle.write(columns[0] + "\t" + columns[1] + "\t" + columns[2] + "\t" + columns[3] + "\t" + columns[4])
    
    input_handle.close()
    output_handle.close()
    
    #range of identity and coverage
    range_cov = str(min(l_cov)) + " - " + str(max(l_cov))
    range_id = str(min(l_iden)) + " - " + str(max(l_iden))
    len_b = len(l_cov) #number of hits
    
    #save the subjects ids in a file in order to extract the sequence from the data base created
    blast_out = open("filtered", "r")
    s_file = open("subjects_ids", "w")

    for line in blast_out:
        columns = line.split("\t")
        s_file.write(columns[1] + '\n')
    
    blast_out.close()
    s_file.close()
    
    try:
        cmd = 'blastdbcmd -entry_batch subjects_ids -db ' + db + ' -outfmt "%s" -out sequences'
        os.system(cmd)
    except:
        print("Could not access data base.")
        sys.exit()

    #join sequences file with blast output file. 
    #Subject sequences are saved in a list and then print in a file with the blast output
    file_seq = open("sequences", "r")
    seq = list() 
    for line in file_seq:
        seq.append(line)
    file_seq.close()

    file_blast = open("filtered", "r")
    output_handle = open("blast_out", "w")
    c = 0
    
    for line in file_blast:
        line = line.strip()
        output_handle.write(line + "\t" + seq[c])
        c += 1

    file_blast.close()
    output_handle.close()

    os.remove(output)
    os.remove("filtered")
    os.remove("sequences")
    os.remove("subjects_ids")

    output_path = "Results/" + output[:-4] + ".tsv"
    shutil.move("blast_out", output_path)

    return output_path, range_cov, range_id, len_b

def query_dict(query_file):
    """
    Create a dictionary with query names as keys and sequence as values
    {query_name: seq}
    """
    query_seq = dict()
    q_seq = str()
    file_q = open(query_file, "r")
    c = 0
    for line in file_q:

        if line.startswith(">"): 
            if c != 0:
                query_seq[q_id] = q_seq
                q_seq = str()

            q_id = line.strip()
            q_id = q_id[1:]
            c = 1

        elif line.startswith("\n") == False: 
            q_seq = q_seq + line.strip()

        else:
            continue
    
    query_seq[q_id] = q_seq
    
    file_q.close()
    return query_seq






