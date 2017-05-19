import os
from SCHEMA_RASPP import schemacontacts, rasppcurve, schemaenergy


class cd:
	# Context manager for changing the current working directory
	def __init__(self, newPath):
		self.newPath = os.path.expanduser(newPath)

	def __enter__(self):
		self.savedPath = os.getcwd()
		os.chdir(self.newPath)

	def __exit__(self, etype, value, traceback):
		os.chdir(self.savedPath)

def generateContacts(inputDirectory):
	with cd(inputDirectory):
		parameters = ['placeholder', '-pdb', 'parent.pdb', '-msa', 'allsequences_aligned.fasta', '-pdbal', 'parent_aligned.fasta', '-o', 'contacts.txt']
		schemacontacts.main(parameters)

def generateRASPPCurve(inputDirectory, nCrossovers, minUnique):
	with cd(inputDirectory):
		parameters = ['placeholder', '-msa', 'allsequences_aligned.fasta', '-con', 'contacts.txt', '-xo', str(nCrossovers), '-o', 'opt.txt', '-min', str(minUnique)]
		rasppcurve.main(parameters)

def computeEnergies(inputDirectory, crossoverSites):
	with cd(inputDirectory):
		parameters = ['placeholder', '-msa', 'allsequences_aligned.fasta', '-con', 'contacts.txt', '-xo', str(crossoverSites), '-E', '-m', '-o', 'energies.txt']
		schemaenergy.main(parameters)
