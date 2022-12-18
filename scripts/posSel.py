import loadFile, PSPFunc, GeneAnalysis, SiteAnalysis, BranchAnalysis, os
from time import localtime, strftime
import logging
import json
import sys,shutil

def pspAnalysis(data, parms, aln, tree):
	"""
	procedure which execute functions for psp step

	@param1 data: basicData object

        @return Output directory name
	"""
	#logger=logging.getLogger("main.positiveSelection")
	dCtrls, lModels = PSPFunc.getParams(parms["models"], 
					    parms["paml"], 
					    parms["bppml"], 
					    parms["mixedlikelihood"], 
					    parms["busted"], 
					    parms["meme"], 
					    parms["opb"], 
					    parms["gnh"])
	timeStamp = strftime("%Y%m%d%H%M", localtime())
	
	outDir = data["o"]+"positive_selection_results_"+timeStamp+"/"
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	
	# cladoFile =  PSPFunc.supBoot(outDir, data["baseName"], tree, logger)<-- Old line, to uncomment if logger is back
	cladoFile =  PSPFunc.supBoot(outDir, data["baseName"], tree)
					
	### Terminal output for user
	"""	
	logger.info("Output directory: {:s}".format(outDir))
	logger.info("Alignement: {:s}".format(aln))
	logger.info("Alignement is in {:s} format.".format(data["alnFormat"]))
	logger.info("Tree: {:s}".format(tree))"""

	### Run the different analysis as determined by control file
	"""
	logger.info("Starting positive selection analyses.")
	logger.info("POSITIVE SELECTION ANALYSIS: ")
	logger.info("Analysis to be run:")"""
	"""
	dAnalysis = {"paml": "Site (codeml)", 
		     "BUSTED":"Whole-Gene", 
		     "bppml":"Site (Bio++ - Optimization)", 
		     "bppmixedlikelihood":"Site (Bio++ - Results)", 
		     "OPB":"Branch", 
		     "GNH":"Branch-site on positively selected branches", 
		     "MEME":"Branch-site"}
	for key in dCtrls.keys():
		logger.info(dAnalysis[key])"""
	
	if "BUSTED" in dCtrls:
		#GeneAnalysis.hyphyBusted(aln, cladoFile, outDir, data["baseName"], logger) <-- Old line, to uncomment if logger is back
		GeneAnalysis.hyphyBusted(aln, cladoFile, outDir, data["baseName"])			
		"""try:		
			GeneAnalysis.hyphyBusted(aln, cladoFile, outDir, data.baseName, logger)
		except Exception:
			logger.info("BUSTED encountered an unexpected error, skipping.")"""		

	if "MEME" in dCtrls:
		try:
			# BranchAnalysis.memeBranchSite(aln, cladoFile, outDir, data["baseName"], logger)<-- Old line, to uncomment if logger is back
			BranchAnalysis.memeBranchSite(aln, 
						      cladoFile, 
						      outDir, 
						      data["baseName"])
		except Exception:
			pass
			#logger.error("MEME encountered an unexpected error, skipping.")
			
	if "bppml" in dCtrls:
#	  try:
            if not dCtrls["bppmixedlikelihood"]:
              dCtrls["bppmixedlikelihood"]=dCtrls["bppml"]
            SiteAnalysis.bppSite(dCtrls["bppml"], 
			         dCtrls["bppmixedlikelihood"], 
			         aln, 
			         data["alnFormat"], 
			         tree, 
			         lModels, 
			         outDir, 
			         data["baseName"])
#	  except Exception:
#	    logger.error("Bio++ Site encountered an unexpected error, skipping.")
	
	lPSNodes = []
	if "OPB" in dCtrls:
#		try:
			params = BranchAnalysis.bppBranch(dCtrls["OPB"], 
							  outDir, 
							  data["baseName"], 
							  aln, 
							  data["alnFormat"], 
							  tree)	
		# except Exception:
		# 	logger.error("Bio++ Branch Analysis encountered an unexpected error, skipping.")
	
	if "OPB" and "GNH" in dCtrls and len(lPSNodes) > 1:
#		try:
			BranchAnalysis.bppBranchSite(dCtrls["GNH"], lPSNodes, outDir, data["baseName"], aln, data["alnFormat"], tree)
		# except Exception:
		# 	logger.error("Bio++ Pseudo Branch-Site Analysis encountered an unexpected error, skipping.")
	
	if "paml" in dCtrls:
		SiteAnalysis.pamlSite(aln, 
				      tree, 
				      lModels, 
				      dCtrls["paml"], 
				      outDir, 
				      data["baseName"])
		"""try:
			SiteAnalysis.pamlSite(aln, tree, lModels, dCtrls["paml"], outDir, data.baseName, logger)
		except Exception:
			logger.info("PAML (codeml) Site encountered an unexpected error, skipping.")"""

	print("Finished positive selection analyses.")
	return(outDir)

if __name__ == "__main__" :	
	with open(sys.argv[1], 'r') as config_in:
		config_dict = json.load(config_in)

	parameters = config_dict["parameters"]
	data = config_dict["data"]
    
	if parameters["positiveSelection"] :  
		if parameters["step"] == "positiveSelection":
			data, data["dAlTree"] = loadFile.pspEntry(data, parameters)
    
	listArgsPosSel =  []
	filename = data["o"]+sys.argv[2]+"_files_list.txt"
	fAT = open(filename, "w")
	for aln in data["dAlTree"]:
		listArgs = [data, parameters, data["aln"], data["dAlTree"][data["aln"]]]
		listArgsPosSel.append(listArgs)
		
		if len(data["dAlTree"][data["aln"]]): # only if it exists
			fAT.write(data["aln"]+"\t"+data["dAlTree"][data["aln"]])
			outDir = pspAnalysis(data, parameters, data["aln"], data["dAlTree"][data["aln"]])
			fAT.write("\t"+outDir+"\n")

	res = shutil.copy(filename, sys.argv[3])
	
	config_dict["parameters"] = parameters
	config_dict["data"] = data

	with open(sys.argv[1],'w') as config_out:
	    json.dump(config_dict, config_out, indent="")

		