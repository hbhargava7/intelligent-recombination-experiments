# Usage: generate_preconditions.py -i /path/to/folder -o /path/for/output

import sys, string, os, glob, subprocess, pprint
from picker import *
from SCHEMA_RASPP import pdb
from shutil import copyfile
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import ClustalOmegaCommandline
import compute_chimeras

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

def computeSequenceIdentityForAlignment(alignmentFile):
	inputFile = open(alignmentFile, "rU")
	alignment = AlignIO.read(inputFile, "clustal")

	j=0 # counts positions in first sequence
	i=0 # counts identity hits 
	for record in alignment:
	    for amino_acid in record.seq:
	        if amino_acid == '-':
	            pass
	        else:
	            if amino_acid == alignment[0].seq[j]:
	                i += 1
	        j += 1
	    j = 0
	    seq = str(record.seq)
	    gap_strip = seq.replace('-', '')
	    percent = 100*i/len(gap_strip)
	    print record.id+' '+str(percent)
	    i=0

def main(args):

	# Obtain arguments parsed into dictionary format.
	parsed_arguments = parse_arguments(args)
        
        # Create the output directory if it does not exist.
        output_path = parsed_arguments["o"]
        if not os.path.exists(output_path):
            os.makedirs(output_path)

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
			if len(glob.glob("*.pdb")) > 0:
				for file in glob.glob("*.pdb"):
					parent["pdb_seq"] = pdb.get(file)
					parent["pdb_id"] = pdb.File().getIDCode(open(file,'r'))
					parent["pdb_path"] = os.getcwd() + "/" + file

			parents.append(parent)

	print("Finished building " + str(len(parents)) + " parents.")

	print("Please select a parent sequence and all additional sequences to compute.")

	pickerOptions = []
	potentialParents = []

	for parentCandidate in parents:
		if "pdb_path" in parentCandidate:
			pickerOptions.append(str(parentCandidate["id"] + " | " + parentCandidate["pdb_id"]))
			potentialParents.append(parentCandidate)

	parentPicker = Picker(
		title = "Please select a single sequence to use as a parent (showing only those with pdb structure file).",
		options = pickerOptions
		).getSelected()

	if parentPicker == False:
		print("Aborting, Parent selection was cancelled.")
		return

	parent = potentialParents[pickerOptions.index(parentPicker[0])]

	pickerOptions = []

	potentialChildren = []
	for childCandidate in parents:
		if "fasta_path" in childCandidate:
			if childCandidate["fasta_path"] is parent["fasta_path"]:
				continue
			pickerOptions.append(childCandidate["id"])
			potentialChildren.append(childCandidate)

	childPicker = Picker(
		title= "Please select the homologous sequences to be used for recombination (structures not required)",
		options = pickerOptions
		).getSelected()
	
	if childPicker == False:
		print("Aborting, sequence selection was cancelled.")
		return

	children = []
	for selectedOption in childPicker:
		children.append(potentialChildren[pickerOptions.index(childPicker[childPicker.index(selectedOption)])])

	print("Successfullly selected Parent: " + parent["id"])
	print("Successfully selected Children: " + str(len(children)))

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

	print("Step 6: Computing sequence identity for aligned sequences.")
	computeSequenceIdentityForAlignment(parsed_arguments["o"] + "/" + "parent_aligned.fasta")

	print("Step 7: Copying the selected parent structure file to output directory.")
	copyfile(parent["pdb_path"], parsed_arguments["o"]+"/parent.pdb")

	print("Preconditions successfully generated! Run the following command to generate contacts:")
	print("python ../SCHEMA_RASPP/schemacontacts.py -pdb parent.pdb -msa allsequences_aligned.fasta -pdbal parent_aligned.fasta -o contacts.txt")

	print("Step 8: Computing contacts from PDB.")
	compute_chimeras.generateContacts(parsed_arguments["o"]+"/")

	print("Step 9: Generating RASPP Curve for specified crossovers.")
	nCrossovers = raw_input("Please specify number of crossover sites to compute: ")
	print("Now computing RASPP curve!")
	compute_chimeras.generateRASPPCurve(parsed_arguments["o"] + "/", nCrossovers, 5)

	print("Step 9: Computing RASPP curve")
	crossoverSites = raw_input("Please review the opt.txt file and input a set of crossover sites: ")
	compute_chimeras.computeEnergies(parsed_arguments["o"] + "/", crossoverSites)

def main_wrapper():
	main(sys.argv)

main_wrapper()
