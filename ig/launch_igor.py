#!/usr/bin/python3

# coding: utf-8

# Filename: launch_igor.py (based on old script 'launch_IGoR_IgG_cohortWide.sh')
# Author: Cosimo Lupo, <cosimo.lupo89@gmail.com>
# Last updated: 11 February 2023

import pandas as pd
import os as os
import yaml

with open('config_igor.yaml') as config_file:
  config = yaml.full_load(config_file)

wk_dir = os.path.abspath(config['wk_dir'])
if wk_dir[-1] != '/':
  wk_dir += '/'

input_dir = os.path.abspath(config['input_dir'])
if input_dir[-1] != '/':
  input_dir += '/'

defaultIgorTemplates = os.path.abspath(config['defaultIgorTemplates'])
if defaultIgorTemplates[-1] != '/':
  defaultIgorTemplates += '/'

inferredIgorTemplates = os.path.abspath(config['inferredIgorTemplates'])
if inferredIgorTemplates[-1] != '/':
  inferredIgorTemplates += '/'

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

seq_file = input_dir + config['cohort'] + "_" + config['cellType'] + "/cohortWide_analysis/" + prefix + ".fasta"
out_dir = input_dir + config['cohort'] + "_" + config['cellType'] + "/cohortWide_analysis/IGoR_inferred_" + prefix + "_SHM" + config['SHMmodel'] + "/"
if not os.path.exists(out_dir):
  os.makedirs(out_dir)

template_path = defaultIgorTemplates + config['species'] + "/"
model_path = defaultIgorTemplates  + config['species'] + "/"
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

  runcmd = config['igorExec']
  runcmd += " -set_wd " + out_dir
  runcmd += " -threads " + str(config['N_cores'])
  #runcmd += " -batch " + batchname
  #runcmd += " -set_genomic --V " + config['genomicV']
  #runcmd += " --D " + config['genomicD']
  #runcmd += " --J " + config['genomicJ']
  runcmd += " -run_custom"

  try:
    res = os.system(runcmd)
    if res != 0:
      raise RuntimeError('IGoR error during the model construction.')
  except BaseException as err:
    raise err

  
################################ SEQ READ ########################################

if config['seqRead']:

  print("# Seq Read")

  # build line command
  runcmd = config['igorExec']
  runcmd += " -set_wd " + out_dir
  runcmd += " -threads " + str(config['N_cores'])
  #runcmd += " -batch " + batchname
  runcmd += " -read_seqs " + seq_file
  if config['N_subsampling'] != None and config['N_subsampling'] != 'None':
    runcmd += " -subsample " + str(config['N_subsampling'])

  try:
    res = os.system(runcmd)
    if res != 0:
      raise RuntimeError('IGoR error when reading input sequences.')
  except BaseException as err:
    raise err
  
  
################################ ALIGNMENT ########################################

if config['alignment']:

  print("# Alignment")

  runcmd = config['igorExec']
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

  try:
    res = os.system(runcmd)
    if res != 0:
      raise RuntimeError('IGoR error during the alignment step.')
  except BaseException as err:
    raise err
  
  
################################ FILTERING OF ALIGNMENT ########################################

if config['filter']:

  print("# Filtering")

  runcmd = "mv " + out_dir + "aligns/indexed_sequences.csv " + out_dir + "aligns/indexed_sequences_old.csv";
  runcmd += "; "
  runcmd += "mv " + out_dir + "aligns/V_alignments.csv " + out_dir + "aligns/V_alignments_old.csv";
  runcmd += "; "
  runcmd += "mv " + out_dir + "aligns/D_alignments.csv " + out_dir + "aligns/D_alignments_old.csv";
  runcmd += "; "
  runcmd += "mv " + out_dir + "aligns/J_alignments.csv " + out_dir + "aligns/J_alignments_old.csv";

  try:
    res = os.system(runcmd)
    if res != 0:
      raise RuntimeError('Bash error when copying alignment files during the filtering step.')
  except BaseException as err:
    raise err

  try:
    res = os.system("python filter_IGoR_align.py " + out_dir + " " + str(config['align_thr']));
    if res != 0:
      raise RuntimeError('Python error during the filtering step of alignments.')
  except BaseException as err:
    raise err
  
  
################################ INFER GENERATIVE MODEL ########################################

if config['inferGenModel']:

  print("# Infer Generative Model")

  runcmd = config['igorExec']
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

  try:
    res = os.system(runcmd)
    if res != 0:
      raise RuntimeError('IGoR error during the inference step.')
  except BaseException as err:
    raise err
  
  # Store the inferred model under 'templates' folder

  runcmd = "mkdir -p " + inferredIgorTemplates + prefix + "_SHM" + config['SHMmodel'] + "/"
  runcmd += "; "
  runcmd += "cp " + out_dir + "inference/final_parms.txt" + " " + inferredIgorTemplates + prefix + "_SHM" + config['SHMmodel'] + "/"
  runcmd += "; "
  runcmd += "cp " + out_dir + "inference/final_marginals.txt" + " " + inferredIgorTemplates + prefix + "_SHM" + config['SHMmodel'] + "/"

  try:
    res = os.system(runcmd)
    if res != 0:
      raise RuntimeError('Bash error when copying IGoR inferred model.')
  except BaseException as err:
    raise err

  
################################ EVALUATE GENERATIVE MODEL ########################################

if config['evalGenModel']:

  print("# Evaluate Generative Model")

  runcmd = config['igorExec']
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

  try:
    res = os.system(runcmd)
    if res != 0:
      raise RuntimeError('IGoR error during the evaluation step.')
  except BaseException as err:
    raise err

  
################################ GENERATE SEQUENCES ########################################

if config['generateSeqs']:

  print("# Generate Seqs")

  runcmd = config['igorExec']
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

  try:
    res = os.system(runcmd)
    if res != 0:
      raise RuntimeError('IGoR error during the generation step.')
  except BaseException as err:
    raise err

  
################################ GATHER IGBLAST + IGOR ANALYSIS TOGETHER ########################################

if config['igAnalysisFinalSummary']:

  print("# Gather igBlast and IGoR analysis together")

  if (not os.path.exists(out_dir + 'inference/final_parms.txt')) or (not os.path.exists(out_dir + 'inference/final_marginals.txt')) or (not os.path.exists(out_dir + 'evaluate/final_parms.txt')) or (not os.path.exists(out_dir + 'evaluate/final_marginals.txt')):
    print("\n" + " ### Error! Any of IGoR inference/evaluate final files is not present! ### " + "\n")
    quit()


  ##### ig seqs import #####

  in_file = input_dir + config['cohort'] + "_" + config['cellType'] + "/cohortWide_analysis/" + prefix + ".csv"
  df_IDs = pd.read_csv(in_file, sep=';', keep_default_na=True)
  # df_IDs headers: seq_ID ; seq_nt
  # df_IDs contains some seq_nt duplicates coming from different individuals (but they are very few)
  # To recover IDs, we drop duplicates, so taking the unique seq_nt from the first individual in the list
  df_IDs_unique = df_IDs.drop_duplicates(subset=['seq_nt'], keep='first')

  print(" Total number of sequences in the set: " + str(len(df_IDs.index)) + " of which " + str(len(df_IDs_unique.index)) + " uniques")


  ##### igBlast analysis #####

  in_file = input_dir + config['cohort'] + "_" + config['cellType'] + "/cohortWide_analysis/" + config['cohort'] + "_" + config['cellType'] + "_" + config['chainType'] + "_uniqueSeqs.igBlast_statistics"
  igBlast_df = pd.read_csv(in_file, sep=';', keep_default_na=True)
  igBlast_df = igBlast_df.add_prefix('igBlast_')


  ##### IGoR seqs indexing #####

  # df_indexedSeqs headers: seq_index ; sequence

  in_file = out_dir + "aligns/" + "indexed_sequences.csv"
  df_indexedSeqs = pd.read_csv(in_file, sep=';')
  df_indexedSeqs['seq_index'] = df_indexedSeqs['seq_index'].astype(int)


  ##### IGoR alignment #####

  # df_*_align headers = seq_index ; gene_name ; score ; offset ; insertions ; deletions ; mismatches ; length ; 5_p_align_offset ; 3_p_align_offset

  in_file = out_dir + "aligns/" + "V_alignments.csv"
  df_V_align = pd.read_csv(in_file, sep=';')

  in_file = out_dir + "aligns/" + "D_alignments.csv"
  if config['chainType']=="HC":
    # Import in chunks dramatically reduces the memory usage
    c_size = 100000
    c_list = []
    for chunk in pd.read_csv(in_file, sep=';', chunksize=c_size):
      chunk.sort_values(by='score', ascending=False, inplace=True)
      chunk.drop_duplicates(subset='seq_index', keep='first', inplace=True)
      c_list.append(chunk)
      df_D_align = pd.concat(c_list)
      df_D_align.sort_values(by='score', ascending=False, inplace=True)
      df_D_align.drop_duplicates(subset='seq_index', keep='first', inplace=True)
      df_D_align.reset_index(drop=True, inplace=True)

  in_file = out_dir + "aligns/" + "J_alignments.csv"
  df_J_align = pd.read_csv(in_file, sep=';')


  ##### IGoR output #####

  # df_Pgen headers = seq_index ; Pgen_estimate

  in_file = out_dir + "output/" + "Pgen_counts.csv"
  df_Pgen = pd.read_csv(in_file, sep=';')

  # df_best_scenarios_counts headers = seq_index ; scenario_rank ; scenario_proba_cond_seq ; GeneChoice_V_gene_Undefined_side_prio8_size97 ; GeneChoice_J_gene_Undefined_side_prio7_size7 ; \
  # GeneChoice_D_gene_Undefined_side_prio6_size35 ; Deletion_D_gene_Five_prime_prio6_size44 ; Deletion_D_gene_Three_prime_prio6_size44 ; Deletion_V_gene_Three_prime_prio5_size27 ; \
  # Deletion_J_gene_Five_prime_prio5_size27 ; Insertion_VD_genes_Undefined_side_prio3_size61 ; DinucMarkov_VD_genes_Undefined_side_prio3_size16 ; \
  # Insertion_DJ_gene_Undefined_side_prio2_size61 ; DinucMarkov_DJ_gene_Undefined_side_prio2_size16 ; Mismatches

  in_file = out_dir + "output/" + "best_scenarios_counts.csv"
  df_best_scenarios_counts = pd.read_csv(in_file, sep=';')


  ##### IGoR evaluation #####

  # df_inf_log headers = iteration_n ; seq_processed ; seq_index ; nt_sequence ; n_V_aligns ; n_J_aligns ; seq_likelihood ; seq_mean_n_errors ; seq_n_scenarios ; seq_best_scenario ; time

  in_file = out_dir + "evaluate/" + "inference_logs.txt"
  df_inf_log = pd.read_csv(in_file, sep=';')


  ##### Merge all the dataframes #####

  main_df = df_indexedSeqs
  main_df['seq_index'] = main_df['seq_index'].astype(int)
  main_df.rename(columns = {'sequence': 'IGoR_seq_nt', 'seq_index': 'IGoR_seq_index'}, inplace = True)
  main_df['IGoR_seq_nt'] = main_df['IGoR_seq_nt'].str.upper()
  main_df['seq_ID'] = main_df['IGoR_seq_nt'].map(df_IDs_unique.set_index('seq_nt')['seq_ID']) # for very few cases, same seq_nt comes from different individuals; we pick the first one in the list
  df_V_align_best = df_V_align.loc[df_V_align.groupby(by="seq_index",sort=True)["score"].idxmax()] # old, memory consuming approach
  if config['chainType']=="HC":
    #df_D_align_best = df_D_align.loc[df_D_align.groupby(by="seq_index",sort=True)["score"].idxmax()]
    df_D_align_best = df_D_align # already parsed when uploading in chunks
  df_J_align_best = df_J_align.loc[df_J_align.groupby(by="seq_index",sort=True)["score"].idxmax()]
  #del df_V_align
  #if config['chainType']=="HC":
    #del df_D_align
  #del df_J_align
  #gc.collect() # Should release memory
  #df_V_align = pd.DataFrame()
  #if config['chainType']=="HC":
    #df_D_align = pd.DataFrame()
  #df_J_align = pd.DataFrame()
  #if config['chainType']=="HC":
    #df_D_align.info()
  df_V_align_best['list_mismatches'] = df_V_align_best['mismatches'].apply(lambda x: [int(i) for i in x[1:-1].split(',')] if (len(x)>2) else [])
  df_V_align_best['N_mismatches'] = df_V_align_best['list_mismatches'].apply(lambda x: len(x))
  df_V_align_best['frac_mismatches'] = df_V_align_best['N_mismatches']/df_V_align_best['length']
  main_df['IGoR_V_best_identity'] = main_df['IGoR_seq_index'].map(df_V_align_best.set_index('seq_index')['gene_name'])
  if config['chainType']=="HC":
    main_df['IGoR_D_best_identity'] = main_df['IGoR_seq_index'].map(df_D_align_best.set_index('seq_index')['gene_name'])
  main_df['IGoR_J_best_identity'] = main_df['IGoR_seq_index'].map(df_J_align_best.set_index('seq_index')['gene_name'])
  main_df['IGoR_V_best_score'] = main_df['IGoR_seq_index'].map(df_V_align_best.set_index('seq_index')['score'])
  main_df['IGoR_V_best_offset'] = main_df['IGoR_seq_index'].map(df_V_align_best.set_index('seq_index')['offset'])
  main_df['IGoR_V_best_length'] = main_df['IGoR_seq_index'].map(df_V_align_best.set_index('seq_index')['length'])
  main_df['IGoR_V_best_listMismatches'] = main_df['IGoR_seq_index'].map(df_V_align_best.set_index('seq_index')['list_mismatches'])
  main_df['IGoR_V_best_NofMismatches'] = main_df['IGoR_seq_index'].map(df_V_align_best.set_index('seq_index')['N_mismatches'])
  main_df['IGoR_V_best_fracMismatches'] = main_df['IGoR_seq_index'].map(df_V_align_best.set_index('seq_index')['frac_mismatches'])
  main_df['IGoR_Pgen'] = main_df['IGoR_seq_index'].map(df_Pgen.set_index('seq_index')['Pgen_estimate'])
  df_best_scenarios_counts = df_best_scenarios_counts.query('scenario_rank==1')
  df_best_scenarios_counts['list_mismatches'] = df_best_scenarios_counts['Mismatches'].apply(lambda x: [int(i) for i in x[1:-1].split(',')] if (len(x)>2) else [])
  df_best_scenarios_counts['N_mismatches'] = df_best_scenarios_counts['list_mismatches'].apply(lambda x: len(x))
  df_best_scenarios_counts['N_mismatches'] = df_best_scenarios_counts['N_mismatches'].astype(int)
  main_df['IGoR_bestScenario_condProba'] = main_df['IGoR_seq_index'].map(df_best_scenarios_counts.set_index('seq_index')['scenario_proba_cond_seq'])
  main_df['IGoR_bestScenario_NofMismatches'] = main_df['IGoR_seq_index'].map(df_best_scenarios_counts.set_index('seq_index')['N_mismatches'])
  main_df['IGoR_seq_likelihood'] = main_df['IGoR_seq_index'].map(df_inf_log.set_index('seq_index')['seq_likelihood'])
  main_df['IGoR_averageNofTotalMismatches'] = main_df['IGoR_seq_index'].map(df_inf_log.set_index('seq_index')['seq_mean_n_errors'])
  main_df['IGoR_bestScenario_likelihood'] = main_df['IGoR_seq_index'].map(df_inf_log.set_index('seq_index')['seq_best_scenario'])
  main_df['IGoR_NofScenarios'] = main_df['IGoR_seq_index'].map(df_inf_log.set_index('seq_index')['seq_n_scenarios'])
  main_df['IGoR_hyperMutations_likelihood'] = main_df['IGoR_seq_likelihood'] / main_df['IGoR_Pgen']

  # Merge with igBlast results (how='left', so to only look at sequences in the IGoR dataframe)
  main_df = main_df.merge(igBlast_df, left_on='seq_ID', right_on='igBlast_seq_ID', how='left')
  # Define productivity of a sequence according to frameshift and/or stop codons in CDR3
  # Current criterion is:
  # NP: V_J_frame=="Out-of-frame" | stop_codon_CDR3=="Yes"
  # P: V_J_frame=="In-frame" & stop_codon_CDR3=="No"
  main_df['myDef_productive'] = main_df.apply(lambda row: "Yes" if (row['igBlast_V_J_frame']=="In-frame" and row['igBlast_stop_codon_CDR3']=="No") else "No", axis=1)
  main_df.drop(['igBlast_seq_ID'], axis=1, inplace=True)

  # Export the final dataframe
  out_file = input_dir + config['cohort'] + "_" + config['cellType'] + "/cohortWide_analysis/" + prefix + "_SHM" + config['SHMmodel'] + ".IGoR_summary"

  main_df.to_csv(out_file, sep=';', index=False)
