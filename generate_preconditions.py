# Usage: generate_preconditions.py -i /path/to/folder -o /path/for/output

import sys, string, os, glob, subprocess, pprint
from SCHEMA_RASPP import pdb
from shutil import copyfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import ClustalOmegaCommandline

class cd:
	# Context manager for changing the current working directory
	def __init__(self, newPath):
		self.newPath = os.path.expanduser(newPath)

	def __enter__(self):
		self.savedPath = os.getcwd()
		os.chdir(self.newPath)

	def __exit__(self, etype, value, traceback):
		os.chdir(self.savedPath)

def parse_arguments(args):
	# Turn linear arguments into a dictionary of (option, [values,...]) pairs
	arg_dict = {}
	key = None
	for arg in args[1:]:
		if arg[0] == '-':
			key = arg[1:]
			arg_dict[key] = None
		else:
			if arg_dict.has_key(key):
				if arg_dict[key]:
					if type(arg_dict[key]) is list:
						arg_dict[key] = arg_dict[key]+[arg]
					else:
						arg_dict[key] = [arg_dict[key],arg]
				else:
					arg_dict[key] = arg
			else:
				arg_dict[key] = arg
	return arg_dict

def main(args):
	# Obtain arguments parsed into dictionary format.
	parsed_arguments = parse_arguments(args)

	print("Step 1: Building potential parents based on input directory.")

	subdirs = next(os.walk(parsed_arguments["i"]))[1]

	# Iterate through subdirectories in input directory and build all possible parents.
	parents = []
	for subdir in subdirs:

		# Ignore hidden directories
		if subdir.startswith("."):
			break

		# Define a dictionary for the parent structure with the format:
		# {"id" : string,
		#  "sequence" : string,
		#  "fasta_path" : string,
		#  "pdb_path" : string,
		#  "pdb_seq" : string,}

		# Use a context to temporarily enter the subdirectory
		with cd(parsed_arguments["i"] + "/" + subdir):

			parent = {}
			# Get the FASTA file containing the sequence.
			for file in glob.glob("*.fasta"):
				# Set the parent sequence to the first sequence in the fasta file.
				for fasta in SeqIO.parse(open(file),'fasta'):
					parent["seq"] = str(fasta.seq)
					parent["id"] = fasta.id
					parent["fasta_path"] = os.getcwd() + "/" + file

			# Get the PDB file containing the structure.
			for file in glob.glob("*.pdb"):
				parent["pdb_seq"] = pdb.get(file)
				parent["pdb_id"] = pdb.File().getIDCode(open(file,'r'))
				parent["pdb_path"] = os.getcwd() + "/" + file
			parents.append(parent)

	print("Finished building " + str(len(parents)) + " parents.")

	print("Please select a parent sequence and all additional sequences to compute.")

	parent = parents[0]
	children = parents[1:]

	print("Step 2: Building FASTA file using all sequences.")
	fasta_files = []
	allsequences = children
	allsequences.append(parent)
	allsequences_path = ""

	for seq in allsequences:
		fasta_files.append(seq["fasta_path"])
	with cd(parsed_arguments["o"]):
		with open('all_unaligned.fasta', 'w') as w_file:
			for filen in fasta_files:
				with open(filen, 'rU') as o_file:
					seq_records = SeqIO.parse(o_file, 'fasta')
					SeqIO.write(seq_records, w_file, 'fasta')
					allsequences_path = os.path.abspath("all_unaligned.fasta")
	
	print("Step 3: Obtaining sequence alignments using ClustalOmega.")

	in_file = allsequences_path
	out_file = parsed_arguments["o"] + "/" + "allsequences_aligned.fasta"
	clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True, force=True, outfmt="clu")
	clustalomega_cline()
	print("Completed sequence alignment.")

	print("Step 4: Building FASTA file using raw and PDB parent sequences.")

	parents_unaligned_path = ""
	with cd(parsed_arguments["o"]):
		with open('parent_unaligned.fasta', 'w') as w_file:
			with open(parent["fasta_path"], 'rU') as o_file:
				seq_records = SeqIO.parse(o_file, 'fasta')
				SeqIO.write(seq_records, w_file, 'fasta')
				parents_unaligned_path = os.path.abspath("parent_unaligned.fasta")

			pdb_sequence = SeqRecord(Seq(parent["pdb_seq"], IUPAC.protein), id=parent["pdb_id"], description="")
			SeqIO.write(pdb_sequence, w_file, 'fasta')

	print("Step 5: Obtaining sequence alignment for master parent using ClustalOmega.")

	in_file = parents_unaligned_path
	out_file = parsed_arguments["o"] + "/" + "parent_aligned.fasta"
	clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True, force=True, outfmt="clu")
	clustalomega_cline()
	print("Completed sequence alignment.")

	print("Step 6: Copying the selected parent structure file to output directory.")
	copyfile(parent["pdb_path"], parsed_arguments["o"]+"/parent.pdb")

	print("Preconditions successfully generated! Run the following command to generate contacts:")
	print("python ../SCHEMA_RASPP/schemacontacts.py -pdb parent.pdb -msa allsequences_aligned.fasta -pdbal parent_aligned.fasta -o contacts.txt")

def main_wrapper():
	main(sys.argv)

main_wrapper()