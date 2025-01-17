import os
import loadFile, tree

import logging, sys
from Bio import Entrez
from statistics import median
import json

"""
This file pools functions related to the creation and conversion of fasta format data.
"""

def dict2fasta(dico):
	"""
	Function converting a dictionary associating gene names to their sequences into a writable Fasta format.

	@param dico: Dictionary associating gene names (keys) to their CCDS fasta sequence (values)
	@return txtoutput: Fasta formatted string of the dictionary
	"""

	txtoutput = ""
	for key, value in dico.items():
		txtoutput += ">{:s}\n{:s}\n".format(str(key),str(value))

	return(txtoutput)


def remoteDl(lBlastRes, queryName, apiKey):
	"""
	Function dowloading species to generate new datas.

	@param1 lBlastRes: List of accessions
	@param2 geneName: Gene name
	@param3 sequence: gene sequence
	@param4 hitsFasta: Path
	@param5 sptree: species tree
	@param6 o: Output directory
	@param7 apiKey: Key for the API of NCBI
	@param8 treerecs: Boolean
	@param9 logger: Object logging
	@return1 outCat: Path to the file containing the sequences and the new IDs
	@return2 corSG: Path
	"""
	#logger = logging.getLogger("main.accessions")
	dSpecies = {}
	dId2Seq = {}
	lTax = []

	Entrez.email="example@example.com"
	Entrez.api_key = apiKey
	
	handle = Entrez.efetch(db="nuccore", id=lBlastRes , idtype="acc", retmode="xml")
	records = list(Entrez.read(handle))
	
	for record in records:
		acc = record['GBSeq_primary-accession']
		tax = record['GBSeq_organism']
		if " x " in tax:
			tax = tax.split(" x ")[0].split(" ")
		elif " X " in tax:
			tax = tax.split(" X ")[0].split(" ")
		else:
			tax = tax.split(" ")
			
		tax = tax[0][:3].lower()+"".join([ i[:3].title() for i in tax[1:]])

		features = [record['GBSeq_feature-table'][i]['GBFeature_quals'] for i, d in enumerate(record['GBSeq_feature-table']) if 'GBFeature_quals' in d]
		
		for feat in features:
			for d in feat:
				if ('GBQualifier_name', 'gene') in d.items():
					name = d['GBQualifier_value']
					break
				else:
					name = ""
					
		if "." in name or "-" in name or name == "":
			try:
				name = "pot"+queryName.split("_")[1]
			except IndexError:
				name = "pot"
		if tax == "synCon" or 'GBSeq_sequence' not in record.keys():
			continue
		else:
			lTax.append(tax)
			dId2Seq[tax+"_"+name+"_"+acc.split(".")[0]] = record['GBSeq_sequence'].upper()
			
	handle.close()
	nbSp = len(set(lTax))
	print(f"Remote option on, downloaded gene IDs and sequences from NCBI databases ({nbSp} different species represented in the retrieved sequences).\n")
	print("Remote option on, downloaded gene IDs and sequences from NCBI databases ({} different species represented in the retrieved sequences).".format(nbSp))
	
	return(dId2Seq)


def sizeCheck(dId2Seq):
	#logger = logging.getLogger("main.accessions")

	dId2Len = {Id:len(seq) for Id, seq in dId2Seq.items()}
	m = median(dId2Len.values())
	n = 0
	
	for k, v in dId2Len.items():
		if m > 10000:
			if v > 2*m:
				try:
					del dId2Seq[k]
					#logger.debug("Deleted sequence {:s} (length {:d})".format(k, v))
					n += 1
				except KeyError:
					pass
		else:
			if v > 3*m or v > 20000:
				try:
					del dId2Seq[k]
					#logger.debug("Deleted sequence {:s} (length {:d})".format(k, v))
					n += 1
				except KeyError:
					pass
	print(f"Deleted {n} sequences due to excessive length\n")
	print("Deleted {} sequences due to excessive length.".format(n))

	return(dId2Seq)
	

def catFile(queryFile, dId2Seq, firstFasta):
	#logger = logging.getLogger("main")

	with open(queryFile, "r") as query:
		lquery = query.readlines()
		query.close()
		dId2Seq[lquery[0].strip().replace(">", "")] = lquery[1]
	
		with open(firstFasta, "w") as fasta:
			fasta.write(dict2fasta(dId2Seq))
			fasta.close()

	return(firstFasta)
	

def fastaCreation(data_dict, remote, apiKey, step, treerecs, outputfile):
	"""
	Function handling the creation of fasta files in the pipeline.

	@param1 data: basicdata object
	@param2 remote: Boolean (online database or not)
	@param3 apiKey: Key for the API of NCBI
	@param4 treerecs: Booleans
	"""

	if remote:
		dId2Seq = remoteDl(data_dict["lBlastRes"],data_dict["queryName"], apiKey)
	else: ### need to code this!!!!
		#logger = logging.getLogger("main.fasta")
		print("Local retrieval of information not yet implemented, exiting DGINN.")
		sys.exit()
	
	dId2Seq = sizeCheck(dId2Seq)
	
	firstFasta = outputfile

	if step == "blast":
		firstFasta = catFile(data_dict["queryFile"], dId2Seq, firstFasta)
	else:
		with open(firstFasta, "w") as out:
			out.write(dict2fasta(dId2Seq))
			out.close()

	data_dict["seqFile"] = firstFasta

	if treerecs:
		
		print(f"spTree : {data_dict}")
		sptree_tmp, treerecs = tree.treeCheck(data_dict["sptree"], firstFasta, treerecs)
		data_dict["sptree"] = sptree_tmp
	if treerecs:
		outCat, corSG = loadFile.filterData(data_dict["sptree"], firstFasta, "results/")
		data_dict["seqFile"] = outCat
		data_dict["cor"] = corSG
	
	return data_dict

if __name__ == "__main__" :

	with open(sys.argv[1], 'r') as config_in:
		config_dict = json.load(config_in)

	data_dict = config_dict["data"]
	params_dict = config_dict["parameters"]

	if params_dict["step"] == "fasta":
		data_dict = loadFile.getSeqEntry(data_dict)
	data_dict["sptree"] = params_dict["sptree"]
	data_dict = fastaCreation(data_dict, params_dict["remote"], params_dict["APIKey"], params_dict["step"], params_dict["duplication"],sys.argv[3])

	config_dict["data"] = data_dict
	config_dict["parameters"] = params_dict
	
	with open(sys.argv[1],'w') as config_out:
		json.dump(config_dict, config_out, indent="")
		