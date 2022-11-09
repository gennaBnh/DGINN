import argparse, logging, os, sys
import json
import FormatFunc, PSPFunc
from time import localtime, strftime
import shutil

version="1.0"

"""This file pools the necessary functions to initialize the pipeline."""

def paramDef(params):
	"""
	Check the parameters in the file.

	@param inf: path's file
	@return defaultParam: dico of parameters
	"""
			   
	#If infile(s) given through command line, takes priority
	params["infile"] = list(map(str.strip,params["infile"].split(",")))
	
	#If list of file given, split and check what each file is
	if len(params["infile"]) > 1:
		for entryfile in params["infile"]:
			if FormatFunc.isCCDSFasta(entryfile):
				params["queryFile"] = os.path.abspath(entryfile)
			if FormatFunc.isAln(entryfile):
				params["alnfile"] = os.path.abspath(entryfile)
			if FormatFunc.isTree(entryfile):
				params["treefile"] = os.path.abspath(entryfile)
	else:
		params["queryFile"] = os.path.abspath(params["infile"][0])
	
	if params["queryFile"] != "":
		params["infile"] = params["queryFile"]	
	elif params["alnfile"] != "":
		params["infile"] = params["alnfile"]	
	
	answers = ["Y", "YES", "T", "TRUE"]
	negAnswers = ["N", "NO", "F", "FALSE"]
	
	for param in params.keys():
		if type(params[param]) is not list:
			if str(params[param]).upper() in answers:
				params[param] = True
			elif str(params[param]).upper() in negAnswers:
				params[param] = False
		
	#Check if parameters are correct
	lSteps = ["blast", "accessions", "fasta", "orf", "alignment", "tree", "duplication", "recombination", "positiveSelection"]
			  
	if params["step"] not in lSteps:
		print("Step \""+params["step"]+"\" not available, set to blast by default.")
		params["step"] = "blast"
	
	if params["remote"] == "":
		print("Remote option needs to be a boolean, set to True by default.")
		params["remote"] = True
	
	if "positiveSelection" not in params:
		print("Positive selection analyses will not be executed, set to False by default.")
		params["positiveSelection"] = False
		
	elif params["positiveSelection"]:
		if params["step"] == "positiveSelection":
			if params["treefile"] == "":
				print("The pipeline requires a phylogenetic tree. Please provide one.")
				sys.exit()
			elif params["alnfile"] == "":
				print("The pipeline requires a codon alignment. Please provide one.")
				sys.exit()
		for opt in ["meme", "busted", "models", "paml", "bppml", "mixedlikelihood", "opb", "gnh"]:
			if opt not in ["meme", "busted", "models", "paml"]:
				if type(params[opt]) is not bool and os.path.exists(params[opt].strip("\n")):
					params[opt] = params[opt].strip("\n")
				elif params[opt]:
					path = "/".join(params["infile"].split("/")[:-1])+"/"+opt+"_params.bpp"
					PSPFunc.pspFileCreation(path, opt)
					params[opt] = path
					
			elif opt == "models":
				ltemp = []
				lmodelok = ["M0", "M1", "M2", "M7", "M8", "M8a", "M10", "DFP07_0", "DFP07"]
				for M in map(str.strip,params[opt].split(",")):
					if M == "":
						next
					elif M not in lmodelok and M[:-2] not in lmodelok:
						print(M + " isn't a valid model.")
					else:
						ltemp.append(M)
				params["models"] = ",".join(ltemp)

	elif params["step"] == "positiveSelection":
		print("Error: positiveSelection option set to false and step set to positiveSelection.")
		sys.exit()
		
	if params["step"] in ["blast","accessions","fasta"]:
		if params["infile"] == "" or params["blastdb"] == "":
			print("Infile and Blastdb are necessary.")
			sys.exit()

	return params


def initLogger(args, debug, version):
	"""
	Function initializing pipeline logger for optimal monitoring.

	@param1 args: Object containing pipeline parameters
	@param2 debug: if the option debug is set
	@param3 version: Pipeline version
	@return1 mainData: Filled basicData object
	@return2 logger: Logging object
	"""
	
	## Log
	### Set up the log directory
	timeStamp = strftime("%Y%m%d%H%M", localtime())
	if args["infile"] != "":
		if args["logfile"] == "":
			args["logfile"] = args["infile"].split(".")[0]+"_DGINN_"+timeStamp+".log"
		else:
			args["logfile"] = args["infile"].split(".")[0]+args["logfile"]
	else:
		if args["logfile"] == "":
			args["logfile"] = args["alnfile"].split(".")[0]+"_DGINN_"+timeStamp+".log"
		else:
			args["logfile"] = args["logfile"]
			
	# create logger
	#logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger("main")
	logger.setLevel(logging.INFO)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(args["logfile"])
	# create console handler with a higher log level
	ch = logging.StreamHandler()
	if debug:
		logger.setLevel(logging.DEBUG)
		ch.setLevel(logging.DEBUG)
		fh.setLevel(logging.DEBUG)
	else:
		ch.setLevel(logging.INFO)
		fh.setLevel(logging.INFO)
	#create formatter and add it to the handlers
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	ch.setFormatter(formatter)
	# add the handlers to the logger
	logger.addHandler(fh)
	logger.addHandler(ch)
	# Welcome message
	logger.info("Starting {:s} (v{:s})".format(__file__, version))

	mainData = {"queryFile": args["infile"],
	"db": args["blastdb"],
	"lBlastRes":[],
	"sptree":args["sptree"],
	"aln":args["alnfile"],
	"tree":args["treefile"],
	"queryName":args["queryName"]
	}

	if args["step"] == "duplication" and args["duplication"] == False:
		args["duplication"] = True
	elif args["step"] == "recombination" and args["recombination"] == False:
		args["recombination"] = True

	logger.info("Reading input file {:s}".format(mainData["queryFile"]))
	logger.info("Analysis will begin at the {:s} step".format(args["step"]))

	return mainData


if __name__ == "__main__" :
	with open(sys.argv[1], 'r') as json_in :
		json_dict = json.loads(json_in.read())
	parameters = json_dict["parameters"]
	print(parameters)
	data = json_dict["data"]

	parameters_complete = paramDef(parameters)
	data_filled = initLogger(parameters_complete, parameters_complete["debug"], version)

	json_dict["parameters"] = parameters_complete
	json_dict["data"] = data_filled
	json_dict_updated = json.dumps(json_dict)
	print(parameters_complete["infile"])

	#filePath = shutil.copy(parameters_complete["infile"][0], '/home/user/doc/')

	with open(sys.argv[2], 'w') as json_out :
   		json_out.write(json_dict_updated)