from snakemake.utils import min_version, validate

### CONFIG FILE DEFINITION ###
config_path = "examples/config_exemple.json" #chemin vers le fichier config de l'utilisateur 
configfile: "examples/config_exemple.json" #chemin du fichier config par défaut
##############################

container: "docker://laugueguen/dginn"
min_version("6.15.1")

#verification du config file, s'il ne correspond pas au schéma, le pipeline ne sera pas lancé.
#permet d'éviter l'ajout de paramètres érronés (ex: infile:200) 
validate(config, "config/config.schema.yaml")

#definition du dossier de stockage des résultats 
dir = config["data"]["o"]

#récupération du nom du infile (ex: "data/filename.fasta")
filename = [config["parameters"]["infile"]][0].split("data/")[1].split(".")[0]
    
#définition du dernier file à produire pour la rule detection car plusieurs fichiers de sortie possible en fonction de la dernière étape.
if config["parameters"]["positiveSelection"] :
    endfile = dir+filename+"_files_list2.txt" 
elif config["parameters"]["recombination"]:
    endfile = dir+filename+"2.gard" 
elif config["parameters"]["duplication"]:
    endfile = dir+filename+"_recs2.nhx"
else: 
    endfile = dir+filename+"_filtered2.phylip_phyml_tree.txt"
    
#règle globale permettant la définition du dernier file à produire par le pipeline. 
#Permet de ne pas avoir à indiquer l'output final:
#snakemake -c1 -p au lieu de snakemake -c1 -p "output_file_name"
rule all:
    input: endfile
        
        
#règle qui permet de lire le fichier config et détermine à partir de ce dernier le premier input du pipeline et la règle qui lui est associée. 
rule init:
    input: 
        config = "examples/config_exemple.json"
    output:
        infile = dir+config["parameters"]["step"]+"_input_"+[config["parameters"]["infile"]][0].split("/")[1]
    message: "\nInitializing pipeline, starting at step "+config["parameters"]["step"]+"."
    shell:
        "python3 scripts/Init.py {input.config} {output.infile}"

#règle qui lance le script Blast.py 
rule blast:
    input:
        infile = dir+"blast_input_{samples}.fasta"
    output:
        tsvfile = dir+"accessions_input_{samples}.tsv"
    params:
        config = config_path
    message: "\nStarting BLAST alignment on file {input.infile}, writing output in file {output.tsvfile}"
    shell:
        "python3 scripts/Blast.py {params.config} {output.tsvfile}"

#règle qui lance le script Extract.py
rule extract:
    input:
        blastres = dir+"accessions_input_{samples}.tsv" 
    output:
        accns = dir+"fasta_input_{samples}.txt" 
    params:
        config = config_path
    message: "\nStarting Accessions search from file {input.blastres}, writing accessions in file {output.accns}"
    shell:
        "python3 scripts/Extract.py {params.config} {input.blastres} {output.accns}"

#règle qui lance le script fastaRes.py 
rule fasta_res:
    input: 
        accns = dir+"fasta_input_{samples}.txt" 
    output:
        seq = dir+"orf_input_{samples}.fasta"
    params:
        config = config_path
    message: "\nStarting the retreival of FASTA sequences from accessions file {input.accns}, writing sequences in file {output.seq}"
    shell:
        "python3 scripts/fastaRes.py {params.config} {input.accns} {output.seq}"

#règle qui lance le script ORF.py 
rule get_orf:
    input:
        seqfile = dir+"orf_input_{samples}.fasta" 
    output:
        res = dir+"{samples}_allORFs.fasta",
        res_longest = dir+"align_input_{samples}.fasta"
    params:
        config = config_path
    message: "\nStarting ORFs search in FASTA sequences from file {input.seqfile}.\nWriting longest ORFs found in file {output.res_longest}"
    shell:
        "python3 scripts/ORF.py {params.config} {output.res} {output.res_longest}"

#règle qui lance l'alignement avec mafft  
rule align_mafft:
    input:
        lgst_ORFs = dir+"align_input_{samples}.fasta"
    output:
        out_covAln = dir+"{samples}_mafft.fasta"
    params:
        config = config_path,
        out_mafft = dir+"{samples}_ORFs_al_mafft.fasta"
    message: "\nStarting mafft nucleotide alignment. Aligning sequences from file {input.lgst_ORFs}, writing alignment in file {output.out_covAln}"
    run:
        shell("mafft --auto --quiet {input.lgst_ORFs} > {params.out_mafft}"),
        shell("python3 scripts/covAln.py {params.out_mafft} {output.out_covAln} {params.config}")
        
#règle qui lance le premier codon aligner choisi par l'utilisateur  
rule first_align_codon:
    input:
        in_alcodon = dir+"{samples}_mafft.fasta" if config["parameters"]["align_nt"] else dir+"align_input_{samples}.fasta"
    output:
        out_clusterIso = dir+"{samples}_clustiso_first.fasta" if len(config["parameters"]["codon_aligner"]) == 2 else dir+"{samples}_clustiso.fasta"
    params:
        config = config_path,
        extension = config["parameters"]["codon_aligner"][0],
        out_codon_param = dir+"{samples}_alcodon",
        macse_param = "-prog refineAlignment -align" if config["parameters"]["align_nt"] else "-prog alignSequences -seq"
    message: "\nStarting first codon alignment with "+config["parameters"]["codon_aligner"][0]+" from file {input.in_alcodon}, writing results in file {output.out_clusterIso}" if len(config["parameters"]["codon_aligner"]) == 2 \
    else "\nStarting codon alignment with "+config["parameters"]["codon_aligner"][0]+" from file {input.in_alcodon}, writing results in file {output.out_clusterIso}"
    run:
        if config["parameters"]["codon_aligner"][0] == "prank" :
            shell("prank -d={input.in_alcodon} -o={params.out_codon_param}_{params.extension}.fas -codon -F"),
            shell("python3 scripts/clusterIso.py {params.config} {output.out_clusterIso} {params.out_codon_param}_{params.extension}.fas")

        elif config["parameters"]["codon_aligner"][0] == "macse" :
            shell("java -jar macse_v2.06.jar {params.macse_param} {input.in_alcodon} -out_NT {params.out_codon_param}_{params.extension}.fas"),
            shell("python3 scripts/clusterIso.py {params.config} {output.out_clusterIso} {params.out_codon_param}_{params.extension}.fas")

        else :
            print(f"Le paramètre codon_aligner du fichier config n'a pas été correctement rempli.\nIl doit contenir une liste d'au moins 1 élément et d'au plus 2.\nLes deux seules valeurs qu'il peut contenir sont \"macse\" ou \"prank\".\nIl contient actuellement la valeur {config['parameters']['codon_aligner']}")

# Règle qui lance le second aligneur de codon, si il est spécifié dans le fichier config. Si un seul aligneur est spécifié, la règle renomme
# simplement le fichier de sortie de la règle précédente pour qu'il soit pris en input par la règle suivante.
rule second_align_codon:
    input:
        in_alcodon_2 = dir+"{samples}_clustiso_first.fasta" if len(config["parameters"]["codon_aligner"]) == 2 else dir+"{samples}_clustiso.fasta"
    output:
        out_clusterIso_2 = dir+"tree_input_{samples}.fasta"
    params:
        config = config_path,
        extension = config["parameters"]["codon_aligner"][1] if len(config["parameters"]["codon_aligner"]) == 2 else config["parameters"]["codon_aligner"][0],
        out_codon_param = dir+"{samples}_alcodon",
        macse_param = "-prog refineAlignment -align" if config["parameters"]["align_nt"] else "-prog alignSequences -seq"
    message: "\nStarting second codon alignment with "+config["parameters"]["codon_aligner"][1]+" from file {input.in_alcodon_2}, writing results in file {output.out_clusterIso_2}"\
    if len(config["parameters"]["codon_aligner"]) == 2 else "\nAlignment {input.in_alcodon_2} renamed : {output.out_clusterIso_2}"
    run:
        if len(config["parameters"]["codon_aligner"]) == 1 :
            shell("mv results/{wildcards.samples}_clustiso.fasta results/tree_input_{wildcards.samples}.fasta") 
            shell("python3 scripts/clusterIso.py {params.config} {output.out_clusterIso_2} {params.out_codon_param}_{params.extension}.fas")
        elif len(config["parameters"]["codon_aligner"]) == 2 :

            if config["parameters"]["codon_aligner"][1] == "prank" :
                shell("prank -d={input.in_alcodon_2} -o={params.out_codon_param}_{params.extension} -codon -F"),
                shell("python3 scripts/clusterIso.py {params.config} {output.out_clusterIso_2} {params.out_codon_param}_{params.extension}.fas")

            elif config["parameters"]["codon_aligner"][1] == "macse" :
                shell("java -jar macse_v2.06.jar {params.macse_param} {input.in_alcodon_2} -out_NT {params.out_codon_param}_{params.extension}.fas"),
                shell("python3 scripts/clusterIso.py {params.config} {output.out_clusterIso_2} {params.out_codon_param}_{params.extension}.fas")
            else :
                print(f"Le paramètre codon_aligner du fichier config n'a pas été correctement rempli.\nIl doit contenir une liste d'au moins 1 élément et d'au plus 2.\nLes deux seules valeurs qu'il peut contenir sont \"macse\" ou \"prank\".\nIl contient actuellement la valeur {config['parameters']['codon_aligner']}")

        else :
            print(f"Le paramètre codon_aligner du fichier config n'a pas été correctement rempli.\nIl doit contenir une liste d'au moins 1 élément et d'au plus 2.\nLes deux seules valeurs qu'il peut contenir sont \"macse\" ou \"prank\".\nIl contient actuellement la valeur {config['parameters']['codon_aligner']}")
  
#règle qui lance la fonction phyMLTree du script Analysis.py 
rule tree:
    input:
        infile = dir+"tree_input_"+filename+".fasta" 
    params:
        fonction = "phyMLTree",
        config = config_path,
    output:
        phymlfile = temp(dir+filename+"_filtered2.phylip_phyml_tree.txt")
    message: "\nStarting Tree building, writing output in file {output.phymlfile}" 
    shell:
        "python3 scripts/Analysis.py {params.config} {params.fonction} {output.phymlfile}"
        
#règle qui lance les étapes "duplication", "recombinaison" et "positiveSelection". Les 3 peuvent être lancés séparement ou à la suite.
#Il était plus simple de les combiner dans une seule règle afin de ne pas avoir à déterminer les différents fichier d'entrée possible.
#De plus, "duplication" et "positiveSelection" prennent en entrée des fichiers similaires. 

rule detection:
    input:
        infile =  dir+"detection_input_"+filename+".fasta" if config["parameters"]["step"] == "detection" else dir+"tree_input_"+filename+".fasta",
        tree =  dir+filename+"_filtered2.phylip_phyml_tree.txt"
    params:
        duplication_function = "checkPhyMLTree",
        recombination_function = "gardRecomb",
        config = config_path,
        filename = filename
    output:
        final_ouput = temp(endfile)
    message: "\nStarting Detection of positive selection, writing output in file {output.final_ouput}" 
    run:
        if config["parameters"]["duplication"]:
            if config["parameters"]["positiveSelection"] or config["parameters"]["recombination"]:
                shell("python3 scripts/Analysis.py {params.config} {params.duplication_function} '' ")
            else: 
                shell("python3 scripts/Analysis.py {params.config} {params.duplication_function} {output.final_ouput}")
        if config["parameters"]["recombination"]:
            if config["parameters"]["positiveSelection"]:
                shell("python3 scripts/Analysis.py {params.config} {params.recombination_function} '' ")
            else :
                shell("python3 scripts/Analysis.py {params.config} {params.recombination_function} {output.final_ouput}")
        if config["parameters"]["positiveSelection"]:
            shell("python3 scripts/posSel.py {params.config} {params.filename} {output.final_ouput}")
        


