#!/usr/bin/python3

# coding: utf-8

# Filename: launch_igor.py (based on old script 'launch_IGoR_IgG_cohortWide.sh')
# Author: Cosimo Lupo, <cosimo.lupo89@gmail.com>
# Last updated: 11 February 2023

import os as os
import yaml

with open('config_igor.yaml') as config_file:
  config = yaml.full_load(config_file)

wkdir = os.path.abspath(config['wkdir'])
if wkdir[-1] != '/':
  wkdir += '/'

if config['SHMmodel'] not in ["none","default","indip","add5mer","full5mer"]:
  print("Wrong SHM model! Please choose among ['none','default','indip','add5mer','full5mer']")
  quit()

if config['chainType'] not in ["HC","KC","LC"]:
  print("Wrong chain type! Please choose among ['HC','KC','LC']")
  quit()

prefix = config['cohort'] + "_" + config['cellType'] + "_" + config['chainType'] + "_uniqueSeqs_" + config['produc']
if config['filter_by_UMI']:
  prefix += "_UMIthr_" + str(config['UMI_thr'])
if config['edges_trimmed']:
  prefix += "_trimmed_" + str(config['n_l']) + "_" + str(config['n_r'])
if config['filter_by_indels']:
  prefix += "_noHyperIndels"

seq_file = wkdir + config['cohort'] + "_" + config['cellType'] + "/cohortWide_analysis/" + prefix + ".fasta"
out_dir = wkdir + config['cohort'] + "_" + config['cellType'] + "/cohortWide_analysis/IGoR_inferred_" + prefix + "_SHM" + config['SHMmodel'] + "/"
if not os.path.exists(out_dir):
  os.makedirs(out_dir)

template_path = config['mainIgorDir'] + "models/" + config['species'] + "/"
model_path = config['mainIgorDir']  + "models/" + config['species'] + "/"
if config['chainType']=="HC":
  template_path += "bcr_heavy/ref_genome/"
  model_path += "bcr_heavy/"
  genomicD = template_path + "genomicDs.fasta"
elif config['chainType']=="KC":
  template_path += "IGK/ref_genome/"
  model_path += "IGK/"
elif config['chainType']=="LC":
  template_path += "IGL/ref_genome/"
  model_path += "IGL/"

genomicV = template_path + "genomicVs.fasta"
genomicJ = template_path + "genomicJs.fasta"
anchorsV = template_path + "V_gene_CDR3_anchors.csv"
anchorsJ = template_path + "J_gene_CDR3_anchors.csv"
  
################################ MODEL CONSTRUCTION ########################################

if config['modelConstruction']:
  
  # Stores the particular Bayesian network

  print("# Model Construction")
  
  runcmd = config['mainIgorDir'] + "igor_src/igor"
  runcmd += " -set_wd " + out_dir
  runcmd += " -threads " + str(config['N_cores'])
  #runcmd += " -batch " + batchname
  #runcmd += " -set_genomic --V " + config['genomicV']
  #runcmd += " --D " + config['genomicD']
  #runcmd += " --J " + config['genomicJ']
  runcmd += " -run_custom"

  os.system(runcmd)

################################ SEQ READ ########################################
  
if config['seqRead']:
  
  print("# Seq Read")
  
  # build line command
  runcmd = config['mainIgorDir'] + "igor_src/igor"
  runcmd += " -set_wd " + out_dir
  runcmd += " -threads " + str(config['N_cores'])
  #runcmd += " -batch " + batchname
  runcmd += " -read_seqs " + seq_file
  if config['N_subsampling'] != None and config['N_subsampling'] != 'None':
    runcmd += " -subsample " + str(config['N_subsampling'])
  
  os.system(runcmd)

################################ ALIGNMENT ########################################
  
if config['alignment']:
  
  print("# Alignment")
  
  runcmd = config['mainIgorDir'] + "igor_src/igor"
  runcmd += " -set_wd " + out_dir
  runcmd += " -threads " + str(config['N_cores'])
  #runcmd += " -batch " + batchname
  if config['SHMmodel']=="none":
    runcmd += " -species " + config['species'] + " -chain heavy_naive"
  elif config['SHMmodel']=="default":
    runcmd += " -species " + config['species'] + " -chain heavy_memory"
  else:
    runcmd += " -set_genomic --V " + genomicV
    if config['chainType']=="HC":
      runcmd += " --D " + genomicD
    runcmd += " --J " + genomicJ
    runcmd += " -set_CDR3_anchors --V " + anchorsV + " --J " + anchorsJ
    runcmd += " -set_custom_model "
    if config['SHMmodel']=="indip":
      runcmd += model_path + "models/model_parms.txt "
      runcmd += model_path + "models/model_marginals.txt"
    elif config['SHMmodel']=="add5mer":
      runcmd += model_path + "supplementary_models/add5mer/V_model/model_parms.txt "
      runcmd += model_path + "supplementary_models/add5mer/V_model/model_marginals.txt"
    elif config['SHMmodel']=="full5mer":
      runcmd += model_path + "supplementary_models/full5mer/V_model/model_parms.txt "
      runcmd += model_path + "supplementary_models/full5mer/V_model/model_marginals.txt"
  runcmd += " -align --all"
  if config['bestGeneOnly']:
    runcmd += " ---best_gene_only true"
  #runcmd += " ---thresh " + str(Sc_thresh)
  runcmd += " ---gap_penalty " + str(config['gap_penalty'])
  #runcmd += " ---matrix " + nuc_matrix
  #runcmd += " --V ---offset_bounds " + str(config['v_low']) + " " + str(config['v_high'])
  #runcmd += " --D ---offset_bounds " + str(config['d_low']) + " " + str(config['d_high'])
  #runcmd += " --J ---offset_bounds " + str(config['j_low']) + " " + str(config['j_high'])
  runcmd += " --V ---thresh " + str(config['V_thresh'])
  if config['chainType']=="HC":
    runcmd += " --D ---thresh " + str(config['D_thresh'])
  runcmd += " --J ---thresh " + str(config['J_thresh'])
  
  os.system(runcmd)

################################ FILTERING OF ALIGNMENT ########################################

if config['filter']:
  
  print("# Filtering")
  
  os.system("mv " + out_dir + "aligns/indexed_sequences.csv " + out_dir + "aligns/indexed_sequences_old.csv");
  os.system("mv " + out_dir + "aligns/V_alignments.csv " + out_dir + "aligns/V_alignments_old.csv");
  os.system("mv " + out_dir + "aligns/D_alignments.csv " + out_dir + "aligns/D_alignments_old.csv");
  os.system("mv " + out_dir + "aligns/J_alignments.csv " + out_dir + "aligns/J_alignments_old.csv");
  
  os.system("python filter_IGoR_align.py " + out_dir + " " + str(config['align_thr']));

################################ INFER GENERATIVE MODEL ########################################

if config['inferGenModel']:
  
  print("# Infer Generative Model")
  
  runcmd = config['mainIgorDir'] + "igor_src/igor"
  runcmd += " -set_wd " + out_dir
  runcmd += " -threads " + str(config['N_cores'])
  #runcmd += " -batch " + batchname
  if config['SHMmodel']=="none":
    runcmd += " -species " + config['species'] + " -chain heavy_naive"
  elif config['SHMmodel']=="default":
    runcmd += " -species " + config['species'] + " -chain heavy_memory"
  else:
    runcmd += " -set_genomic --V " + genomicV
    if config['chainType']=="HC":
      runcmd += " --D " + genomicD
    runcmd += " --J " + genomicJ
    runcmd += " -set_CDR3_anchors --V " + anchorsV + " --J " + anchorsJ
    runcmd += " -set_custom_model "
    if config['SHMmodel']=="indip":
      runcmd += model_path + "models/model_parms.txt "
      runcmd += model_path + "models/model_marginals.txt"
    elif config['SHMmodel']=="add5mer":
      if config['initInferBySHMindip']:
        runcmd += out_dir + "../IGoR_inferred_" + prefix + "_SHMindip/inference/final_parms_for_SHMadd5mer.txt "
        runcmd += out_dir + "../IGoR_inferred_" + prefix + "_SHMindip/inference/final_marginals_for_SHMadd5mer.txt"
      else:
        runcmd += model_path + "supplementary_models/add5mer/V_model/model_parms.txt "
        runcmd += model_path + "supplementary_models/add5mer/V_model/model_marginals.txt"
    elif config['SHMmodel']=="full5mer":
      runcmd += model_path + "supplementary_models/full5mer/V_model/model_parms.txt "
      runcmd += model_path + "supplementary_models/full5mer/V_model/model_marginals.txt"
  runcmd += " -infer"
  #runcmd += " --not_infer ErrorRate"
  #runcmd += " --infer_only SingleErrorRate"
  runcmd += " --N_iter " + str(config['N_iters'])
  runcmd += " --L_thresh " + config['L_thresh']
    
  os.system(runcmd)

################################ EVALUATE GENERATIVE MODEL ########################################
  
if config['evalGenModel']:
  
  print("# Evaluate Generative Model")
  
  runcmd = config['mainIgorDir'] + "igor_src/igor"
  runcmd += " -set_wd " + out_dir
  runcmd += " -threads " + str(config['N_cores'])
  #runcmd += " -batch " + batchname
  if not config['useLastInferred']:
    if config['SHMmodel']=="none":
      runcmd += " -species " + config['species'] + " -chain heavy_naive"
    elif config['SHMmodel']=="default":
      runcmd += " -species " + config['species'] + " -chain heavy_memory"
    else:
      runcmd += " -set_genomic --V " + genomicV
      if config['chainType']=="HC":
        runcmd += " --D " + genomicD
      runcmd += " --J " + genomicJ
      runcmd += " -set_CDR3_anchors --V " + anchorsV + " --J " + anchorsJ
      runcmd += " -set_custom_model "
      if config['SHMmodel']=="indip":
        runcmd += model_path + "models/model_parms.txt "
        runcmd += model_path + "models/model_marginals.txt"
      elif config['SHMmodel']=="add5mer":
        runcmd += model_path + "supplementary_models/add5mer/V_model/model_parms.txt "
        runcmd += model_path + "supplementary_models/add5mer/V_model/model_marginals.txt"
      elif config['SHMmodel']=="full5mer":
        runcmd += model_path + "supplementary_models/full5mer/V_model/model_parms.txt "
        runcmd += model_path + "supplementary_models/full5mer/V_model/model_marginals.txt"
  elif config['useLastInferred']:
    runcmd += " -set_genomic --V " + genomicV
    if config['chainType']=="HC":
      runcmd += " --D " + genomicD
    runcmd += " --J " + genomicJ
    runcmd += " -set_CDR3_anchors --V " + anchorsV + " --J " + anchorsJ
    runcmd += " -set_custom_model "
    runcmd += out_dir + "inference/final_parms.txt "
    runcmd += out_dir + "inference/final_marginals.txt"
  runcmd += " -evaluate"
  runcmd += " --L_thresh " + config['L_thresh']
  runcmd += " -output --Pgen"
  runcmd += " --scenarios " + str(config['N_scen'])

  os.system(runcmd)

################################ GENERATE SEQUENCES ########################################
  
if config['generateSeqs']:

  print("# Generate Seqs")
    
  runcmd = config['mainIgorDir'] + "igor_src/igor"
  runcmd += " -set_wd " + out_dir
  runcmd += " -threads " + str(config['N_cores'])
  #runcmd += " -batch " + batchname
  if not config['useLastEvaluated']:
    if config['SHMmodel']=="none":
      runcmd += " -species " + config['species'] + " -chain heavy_naive"
    elif config['SHMmodel']=="default":
      runcmd += " -species " + config['species'] + " -chain heavy_memory"
    else:
      runcmd += " -set_genomic --V " + genomicV
      if config['chainType']=="HC":
        runcmd += " --D " + genomicD
      runcmd += " --J " + genomicJ
      runcmd += " -set_CDR3_anchors --V " + anchorsV + " --J " + anchorsJ
      runcmd += " -set_custom_model "
      if config['SHMmodel']=="indip":
        runcmd += model_path + "models/model_parms.txt "
        runcmd += model_path + "models/model_marginals.txt"
      elif config['SHMmodel']=="add5mer":
        runcmd += model_path + "supplementary_models/add5mer/V_model/model_parms.txt "
        runcmd += model_path + "supplementary_models/add5mer/V_model/model_marginals.txt"
      elif config['SHMmodel']=="full5mer":
        runcmd += model_path + "supplementary_models/full5mer/V_model/model_parms.txt "
        runcmd += model_path + "supplementary_models/full5mer/V_model/model_marginals.txt"
  elif config['useLastEvaluated']:
    runcmd += " -set_genomic --V " + genomicV
    if config['chainType']=="HC":
      runcmd += " --D " + genomicD
    runcmd += " --J " + genomicJ
    runcmd += " -set_CDR3_anchors --V " + anchorsV + " --J " + anchorsJ
    runcmd += " -set_custom_model "
    runcmd += out_dir + "evaluate/final_parms.txt "
    runcmd += out_dir + "evaluate/final_marginals.txt"
  runcmd += " -generate " + str(config['N_generated'])
  runcmd += " --noerr"

  os.system(runcmd)
