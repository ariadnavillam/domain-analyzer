# -*- coding: utf-8 *-*
import re
from Bio.ExPASy import Prosite, Prodoc

def c_patterns():
	"""
	Create a dictionary with names and pattern of the each record.
	And convert prosite patterns into REGEX pattern
	"""
	handle = open("prosite.dat","r")
	records = Prosite.parse(handle)
	pattern_dict = dict()
	for record in records:
		if record.pattern != "":
			pattern_dict[record.accession] = [record.pattern, record.name, record.description]
	handle.close()

	#change prosite patters to REGEX patterns
	for key in pattern_dict:
		p = pattern_dict[key][0]
		p = p.replace("-", "")
		p = p.replace("x", ".")
		p = p.replace("(", "{")
		p = p.replace(")", "}")
		pattern_dict[key][0] = p
	
	return pattern_dict
	
def find_patterns(input_file, p_dict):
	"""
	For each sequence in a fasta file this function finds domains from prosite 
	present in the protein sequence.
	Output: File with the protein and domain information. 
	Each protein may have more than one domain, but each domain is printed in one line.
	Returns a list of all domains found (later used to make a graphic) 
	"""
	input_handle = open(input_file, "r")
	output_handle = open(input_file[:-3] + "_domains.txt", "w")
	output_handle.write("Protein ID(|Organism)\tDomain name\tDomain accession\tDomain description\tPattern found\n")
	f_patterns = list()
	for line in input_handle:
		if line.startswith(">"):
			prot_name = line.strip()

		else:
			sequence = line.strip()
			for key in p_dict:
				p_regex = re.compile(p_dict[key][0])
				p_found = p_regex.search(sequence)
				if  p_found != None:
					f_patterns.append(p_dict[key][1])
					output_handle.write(prot_name[1:] + "\t" + p_dict[key][1] + "\t" + key + "\t" + p_dict[key][2][:-1] + "\t" + str(p_found.group()) +"\n")
	
	return f_patterns
				


