#!/usr/bin/python3

# coding: utf-8

# Filename: launch_igor.py (based on old script 'launch_IGoR.py')
# Author: Cosimo Lupo, <cosimo.lupo89@gmail.com>
# Last updated: 11 February 2023

import pandas as pd
import numpy as np
from scipy.stats import poisson
import glob as glob
import os as os
import yaml

#sys.path.insert(1, '../libraries/')
from aux_funcs_igor import expand_array, expand_array_2, expand_array_3
from aux_funcs_igor import get_CDR3_len, get_family, parse_hyper_indels

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

# Project specs

seq_file = input_dir + config['input_file_prefix'] + '_' + config['chainType']
if config['filter_by_indels']:
  seq_file += '_noHyperIndels'
if config['edges_trimmed']:
  seq_file += '_trimmed_' + str(config['n_l']) + '_' + str(config['n_r'])
seq_file += '.fasta'

template_path = defaultIgorTemplates + config['species'] + "/"
model_path = defaultIgorTemplates + config['species'] + "/"
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

# SHM reference model (IGoR standard model + cohortWide models inferred on Klein data)

if config['singleModel']:
  cohorts = [config['cohort']]
  producs = [config['produc']]
  SHMmodels = [config['SHMmodel']]
else:
  cohorts = config['cohorts']
  producs = config['producs']
  SHMmodels = config['SHMmodels']

for (cohort,produc,SHMmodel) in [(c,p,s) for c in cohorts for p in producs for s in SHMmodels]:
  
  print('\n\n' + 'Running: ' + config['chainType'] + ' ' + cohort + ' ' + produc + ' ' + SHMmodel + '\n')
  
  out_dir = input_dir + "igor_"
  
  if produc=="standard_NP":
    out_dir += "standard_" + config['chainType'] + "_NP_SHMindip" + "/"
  else:
    model_path = inferredIgorTemplates + cohort + "_IgG_" + config['chainType'] + "_uniqueSeqs_"
    out_dir += "cohortWide_" + cohort + "_" + config['chainType'] + "_"
    if produc=="cohortWide_NP":
      model_path += "NP"
      out_dir += "NP"
    elif produc=="cohortWide_P":
      model_path += "P"
      out_dir += "P"
    else:
      print("wrong model!")
      quit()
    model_path += "_UMIthr_" + str(config['UMI_thr']) + "_trimmed_*_*_noHyperIndels_SHM" + SHMmodel + "/"
    out_dir += "_SHM" + SHMmodel + "/"
    model_paths = glob.glob(model_path)
    try:
      model_path = os.path.abspath(model_paths[0])
      if model_path[-1] != '/':
        model_path += '/'
    except BaseException as err:
      print(err)
  
  if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
  
  if config['chainType']=="HC":
    igor_chain_name = "heavy_memory"
  if config['chainType']=="KC":
    igor_chain_name = "IGK"
  if config['chainType']=="LC":
    igor_chain_name = "IGL"
    
  ################################ SEQ READ ########################################
  
  if config['seqRead']:
    
    print("# Seq Read")
    
    runcmd = config['igorExec']
    runcmd += " -set_wd " + out_dir
    runcmd += " -threads " + str(config['N_cores'])
    #runcmd += " -batch " + batchname
    runcmd += " -read_seqs " + seq_file
    #if config['N_subsampling'] != None and config['N_subsampling'] != 'None':
    #  runcmd += " -subsample " + str(config['N_subsampling'])
    
    os.system(runcmd)
  
  ################################ ALIGNMENT ########################################
  
  if config['alignment']:
    
    print("# Alignment")
    
    runcmd = config['igorExec']
    runcmd += " -set_wd " + out_dir
    runcmd += " -threads " + str(config['N_cores'])
    #runcmd += " -batch " + batchname
    if produc=="standard_NP":
      runcmd += " -species human -chain " + igor_chain_name
    else:
      runcmd += " -set_genomic --V " + genomicV
      if config['chainType']=="HC":
        runcmd += " --D " + genomicD
      runcmd += " --J " + genomicJ
      runcmd += " -set_CDR3_anchors --V " + anchorsV + " --J " + anchorsJ
      runcmd += " -set_custom_model " + model_path + "final_parms.txt " + model_path + "final_marginals.txt"
    runcmd += " -align --all"
    if config['bestGeneOnly']:
      runcmd += " ---best_gene_only true"
    #runcmd += " ---thresh " + str(Sc_thresh)
    runcmd += " ---gap_penalty " + str(config['gap_penalty'])
    #runcmd += " ---matrix " + nuc_matrix
    #runcmd += " --V ---offset_bounds " + str(v_low) + " " + str(v_high)
    #runcmd += " --D ---offset_bounds " + str(d_low) + " " + str(d_high)
    #runcmd += " --J ---offset_bounds " + str(j_low) + " " + str(j_high)
    runcmd += " --V ---thresh " + str(config['V_thresh'])
    if config['chainType']=="HC":
      runcmd += " --D ---thresh " + str(config['D_thresh'])
    runcmd += " --J ---thresh " + str(config['J_thresh'])
    
    os.system(runcmd)
  
  ################################ EVALUATE GENERATIVE MODEL ########################################
  
  if config['evalGenModel']:
    
    print("# Evaluate Generative Model")
    
    runcmd = config['igorExec']
    runcmd += " -set_wd " + out_dir
    runcmd += " -threads " + str(config['N_cores'])
    #runcmd += " -batch " + batchname
    if produc=="standard_NP":
      runcmd += " -species human -chain " + igor_chain_name
    else:
      runcmd += " -set_genomic --V " + genomicV
      if config['chainType']=="HC":
        runcmd += " --D " + genomicD
      runcmd += " --J " + genomicJ
      runcmd += " -set_CDR3_anchors --V " + anchorsV + " --J " + anchorsJ
      runcmd += " -set_custom_model " + model_path + "final_parms.txt " + model_path + "final_marginals.txt"
    runcmd += " -evaluate"
    runcmd += " --L_thresh " + config['L_thresh']
    runcmd += " -output --Pgen"
    runcmd += " --scenarios " + str(config['N_scen'])
    
    os.system(runcmd)
    
    
  ################################ GATHER IGBLAST + IGOR ANALYSIS TOGETHER ########################################
  
  if config['bnabAnalysisFinalSummary']:

    print("# Gather igBlast and IGoR analysis together")

    if (not os.path.exists(out_dir + 'evaluate/final_parms.txt')) or (not os.path.exists(out_dir + 'evaluate/final_marginals.txt')):
      print("\n" + " ### Attention! Any of IGoR evaluate final files is not present! ### " + "\n")

    
    ##### bnabs seqs import #####
    
    in_file = input_dir + config['input_file_prefix'] + "_" + config['chainType'] + ".csv"
    main_df = pd.read_csv(in_file, sep=";", skiprows=1, header=None, names=['bnab_ID','orig_seq_nt'])
    #main_df.rename(columns = {'Antibodies':'bnab_ID', 'Heavy_chain':'orig_seq_nt'}, inplace = True)
    main_df = main_df[main_df['orig_seq_nt'].notnull()].reset_index(drop=True)
    main_df['orig_seq_nt'] = main_df['orig_seq_nt'].str.replace('-','')
    main_df['orig_seq_nt'] = main_df['orig_seq_nt'].str.upper()
    
    in_file = input_dir + config['input_file_prefix'] + "_" + config['chainType'] + "_noHyperIndels.csv"
    noHyperIndels_seqs_df = pd.read_csv(in_file, sep=";", skiprows=1, header=None, names=['bnab_ID','noHyperIndels_seq_nt'])
    #noHyperIndels_seqs_df.rename(columns = {'seq_ID':'bnab_ID', 'noHyperIndels_seq_nt':'noHyperIndels_seq_nt'}, inplace = True)
    noHyperIndels_seqs_df['noHyperIndels_seq_nt'] = noHyperIndels_seqs_df['noHyperIndels_seq_nt'].str.replace('-','')
    noHyperIndels_seqs_df['noHyperIndels_seq_nt'] = noHyperIndels_seqs_df['noHyperIndels_seq_nt'].str.upper()
    
    in_file = input_dir + config['input_file_prefix'] + "_" + config['chainType'] + "_noHyperIndels_trimmed_" + str(config['n_l']) + '_' + str(config['n_r']) + ".csv"
    noHyperIndels_trimmed_seqs_df = pd.read_csv(in_file, sep=";", skiprows=1, header=None, names=['bnab_ID','noHyperIndels_trimmed_seq_nt'])
    #noHyperIndels_trimmed_seqs_df.rename(columns = {'Antibodies':'bnab_ID', 'Heavy_chain':'noHyperIndels_trimmed_seq_nt'}, inplace = True)
    noHyperIndels_trimmed_seqs_df['noHyperIndels_trimmed_seq_nt'] = noHyperIndels_trimmed_seqs_df['noHyperIndels_trimmed_seq_nt'].str.replace('-','')
    noHyperIndels_trimmed_seqs_df['noHyperIndels_trimmed_seq_nt'] = noHyperIndels_trimmed_seqs_df['noHyperIndels_trimmed_seq_nt'].str.upper()
    
    main_df['noHyperIndels_seq_nt'] = main_df['bnab_ID'].map(noHyperIndels_seqs_df.set_index('bnab_ID')['noHyperIndels_seq_nt'])
    main_df['noHyperIndels_trimmed_seq_nt'] = main_df['bnab_ID'].map(noHyperIndels_trimmed_seqs_df.\
                                                                   set_index('bnab_ID')['noHyperIndels_trimmed_seq_nt'])
    
    main_df['orig_seq_nt'] = main_df['orig_seq_nt'].astype(str)
    main_df['orig_seq_nt_length'] = main_df['orig_seq_nt'].apply(len)
    main_df['noHyperIndels_seq_nt'] = main_df['noHyperIndels_seq_nt'].astype(str)
    main_df['noHyperIndels_seq_nt_length'] = main_df['noHyperIndels_seq_nt'].apply(len)
    main_df['noHyperIndels_trimmed_seq_nt'] = main_df['noHyperIndels_trimmed_seq_nt'].astype(str)
    main_df['noHyperIndels_trimmed_seq_nt_length'] = main_df['noHyperIndels_trimmed_seq_nt'].apply(len)
    
    ##### bnab binding site & neutralization data #####
    
    in_file = input_dir + "bnabs.csv"
    features_df = pd.read_csv(in_file, sep=";")
    
    main_df['binding_site'] = main_df['bnab_ID'].map(features_df.set_index('bnab_ID')['binding_site'])
    main_df['AuC'] = main_df['bnab_ID'].map(features_df.set_index('bnab_ID')['AuC'])
    #main_df['AuC'] = [x.replace(',', '.') for x in main_df['Area_under_the_curve']]
    main_df['AuC'] = main_df['AuC'].apply(lambda x: float(x))
    
    ##### igBlast analysis #####
    
    in_file = input_dir + config['input_file_prefix'] + "_" + config['chainType'] + ".igBlast_statistics"
    igBlast_df = pd.read_csv(in_file, sep=';', keep_default_na=True)
    igBlast_df = igBlast_df.add_prefix('igBlast_')
    
    main_df = main_df.merge(igBlast_df, left_on='bnab_ID', right_on='igBlast_bnab_ID', how='inner')
    main_df['myDef_productive'] = main_df.apply(lambda row: "Yes" if (row['igBlast_V_J_frame']=="In-frame" and \
                                                                    row['igBlast_stop_codon_CDR3']=="No") else "No", axis=1)
    main_df.drop(['igBlast_bnab_ID'], axis=1, inplace=True)
    
    main_df['igBlast_hyper_indels_total_L'] = main_df.apply(lambda row: row['igBlast_Total_L_ins']-row['igBlast_Total_L_del'], \
                                                          axis=1)
    main_df['igBlast_V_best_align_effective_length_beforeCDR3'] = main_df.apply(lambda row: \
                                                                              row['igBlast_V_best_align_length_beforeCDR3'] - \
                                                                              row['igBlast_hyper_indels_total_L'], axis=1)
    mu_hyperInsertions = 0.00039
    mu_hyperDeletions = 0.00068
    main_df['likelihood_N_of_hyperins'] = main_df.apply(lambda row: poisson.pmf(row['igBlast_N_ins'],row['igBlast_V_best_align_effective_length_beforeCDR3']*mu_hyperInsertions) if row['igBlast_V_best_align_effective_length_beforeCDR3'] is not np.nan else np.nan, axis=1)
    main_df['likelihood_N_of_hyperdel'] = main_df.apply(lambda row: poisson.pmf(row['igBlast_N_del'],row['igBlast_V_best_align_effective_length_beforeCDR3']*mu_hyperDeletions) if row['igBlast_V_best_align_effective_length_beforeCDR3'] is not np.nan else np.nan, axis=1)
    
    #igBlast_gaps_elaborated['hyperindel_data'] = igBlast_gaps_elaborated['hyperindel_data'].apply(expand_array_3)
    #igBlast_gaps_elaborated['inter_hyperindels_spacings'] = igBlast_gaps_elaborated['inter_hyperindels_spacings'].apply(expand_array_2)
    
    #main_df['igBlast_V_best_align_identity'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('bnab_ID')['V_best'])
    #main_df['igBlast_V_best_align_identity'] = main_df['igBlast_V_best_align_identity'].astype(str)
    #main_df['igBlast_V_best_align_family'] = main_df['igBlast_V_best_align_identity'].apply(get_family)
    #if(chainType=="HC"):
      #main_df['igBlast_D_best_align_identity'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['D_best'])
      #main_df['igBlast_D_best_align_identity'] = main_df['igBlast_D_best_align_identity'].astype(str)
      #main_df['igBlast_D_best_align_family'] = main_df['igBlast_D_best_align_identity'].apply(get_family)
    #main_df['igBlast_J_best_align_identity'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['J_best'])
    #main_df['igBlast_J_best_align_identity'] = main_df['igBlast_J_best_align_identity'].astype(str)
    #main_df['igBlast_J_best_align_family'] = main_df['igBlast_J_best_align_identity'].apply(get_family)
    #main_df['igBlast_stopCodon_presence'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['stop_codon'])
    #main_df['igBlast_V_J_frame'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_J_frame'])
    #main_df['igBlast_productivity'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['productive'])
    #main_df['igBlast_strand'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['strand'])
    #main_df['igBlast_CDR3_identity'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['CDR3'])
    #main_df['igBlast_CDR3_identity'] = main_df['igBlast_CDR3_identity'].astype(str)
    #main_df['igBlast_CDR3_length'] = main_df['igBlast_CDR3_identity'].apply(get_CDR3_len)
    
    #main_df['igBlast_V_best_align_fracMatches'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_frac_of_matches'])
    #main_df['igBlast_V_best_align_length'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_align_length'])
    ##main_df['igBlast_V_best_align_length'] = main_df['igBlast_V_best_align_length'].astype(int)
    #main_df['igBlast_V_best_align_length_beforeCDR3'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_align_length_beforeCDR3'])
    ##main_df['igBlast_V_best_align_length_beforeCDR3'] = main_df['igBlast_V_best_align_length_beforeCDR3'].astype(int)
    #main_df['igBlast_V_best_align_effective_length_beforeCDR3'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_align_effective_length_beforeCDR3'])
    #main_df['igBlast_V_best_align_NofMismatches'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_N_of_mm'])
    ##main_df['igBlast_V_best_align_NofMismatches'] = main_df['igBlast_V_best_align_NofMismatches'].astype(int)
    #main_df['igBlast_V_best_align_fracMismatches'] = main_df['igBlast_V_best_align_NofMismatches'] / main_df['igBlast_V_best_align_effective_length_beforeCDR3']
    #main_df['igBlast_V_best_align_seqStartPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_align_start_seq'])
    #main_df['igBlast_V_best_align_seqEndPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_align_end_seq'])
    #main_df['igBlast_V_best_align_geneStartPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_align_start_gene'])
    #main_df['igBlast_V_best_align_geneEndPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_align_end_gene'])
    #main_df['igBlast_V_best_align_score'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_igBlast_score'])
    
    #if(chainType=="HC"):
      #main_df['igBlast_D_best_align_fracMatches'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['D_best_frac_of_matches'])
      #main_df['igBlast_D_best_align_length'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['D_best_align_length'])
      ##main_df['igBlast_D_best_align_length'] = main_df['igBlast_D_best_align_length'].astype(int)
      #main_df['igBlast_D_best_align_NofMismatches'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['D_best_N_of_mm'])
      ##main_df['igBlast_D_best_align_NofMismatches'] = main_df['igBlast_D_best_align_NofMismatches'].astype(int)
      ##main_df['igBlast_D_best_align_fracMismatches'] = np.where(main_df['igBlast_D_best_align_NofMismatches']=="none","none",main_df['igBlast_D_best_align_NofMismatches'] / main_df['igBlast_D_best_align_length'])
      #main_df['igBlast_D_best_align_seqStartPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['D_best_align_start_seq'])
      ##main_df['igBlast_D_best_align_seqStartPosition'] = main_df['igBlast_D_best_align_seqStartPosition'].astype(int)
      #main_df['igBlast_D_best_align_seqEndPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['D_best_align_end_seq'])
      #main_df['igBlast_D_best_align_geneStartPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['D_best_align_start_gene'])
      #main_df['igBlast_D_best_align_geneEndPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['D_best_align_end_gene'])
      #main_df['igBlast_D_best_align_score'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['D_best_igBlast_score'])
    
    #main_df['igBlast_J_best_align_fracMatches'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['J_best_frac_of_matches'])
    #main_df['igBlast_J_best_align_length'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['J_best_align_length'])
    ##main_df['igBlast_J_best_align_length'] = main_df['igBlast_J_best_align_length'].astype(int)
    #main_df['igBlast_J_best_align_NofMismatches'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['J_best_N_of_mm'])
    ##main_df['igBlast_J_best_align_NofMismatches'] = main_df['igBlast_J_best_align_NofMismatches'].astype(int)
    ##main_df['igBlast_J_best_align_fracMismatches'] = main_df['igBlast_J_best_align_NofMismatches'] / main_df['igBlast_J_best_align_length']
    #main_df['igBlast_J_best_align_seqStartPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['J_best_align_start_seq'])
    ##main_df['igBlast_J_best_align_seqStartPosition'] = main_df['igBlast_J_best_align_seqStartPosition'].astype(int)
    #main_df['igBlast_J_best_align_seqEndPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['J_best_align_end_seq'])
    #main_df['igBlast_J_best_align_geneStartPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['J_best_align_start_gene'])
    #main_df['igBlast_J_best_align_geneEndPosition'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['J_best_align_end_gene'])
    #main_df['igBlast_J_best_align_score'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['J_best_igBlast_score'])
    
    #main_df['igBlast_V_best_align_NofMismatches_fromPairwiseCompar'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_N_of_mm_fromPairwiseCompar'])
    #main_df['igBlast_V_best_align_ListofMismatches_fromPairwiseCompar'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['V_best_list_of_mm_fromPairwiseCompar'])
      
    main_df['igBlast_errors_list'] = main_df['igBlast_errors_list'].apply(lambda x: x.replace('[','').replace(']','').replace(' ','').split(',') if len(x)>2 else [])
    main_df['igBlast_deletions'] = main_df['igBlast_errors_list'].apply(lambda x: parse_hyper_indels(x,"del"))
    main_df['igBlast_insertions'] = main_df['igBlast_errors_list'].apply(lambda x: parse_hyper_indels(x,"ins"))
    #main_df['igBlast_insertions'] = main_df['igBlast_hyper_indels_list'].apply(lambda x: [x[i * 4:(i + 1) * 4] for i in range((len(x) + 4 - 1) // 4 )] if len(x)>0 else [])
    
    #main_df['igBlast_NofHyperInsertions'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['N_of_hyperins'])
    #main_df['igBlast_NofHyperDeletions'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['N_of_hyperdel'])
    #main_df['igBlast_NofHyperIndels'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['N_of_hyperindels'])
    #main_df['igBlast_HyperInsertions_length'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['total_hyperins_length'])
    #main_df['igBlast_HyperDeletions_length'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['total_hyperdel_length'])
    #main_df['igBlast_HyperIndels_length'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['total_hyperindels_length'])
    #main_df['igBlast_HyperIndels_info'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['hyperindel_data'])
    #main_df['igBlast_HyperIndels_interDistance'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['inter_hyperindels_spacings'])
    
    l_hyperInsertions = 3.32
    l_hyperDeletions = 5.28
    
    main_df['igBlast_del_total_likelihood'] = main_df['igBlast_deletions'].apply(lambda x: [np.exp(-i/l_hyperDeletions)/l_hyperDeletions for i in x])
    main_df['igBlast_del_total_likelihood'] = main_df.apply(lambda row: np.prod(row['igBlast_del_total_likelihood']) if row['igBlast_N_del']>0 else 1, axis=1)
    main_df['igBlast_ins_total_likelihood'] = main_df['igBlast_insertions'].apply(lambda x: [np.exp(-i/l_hyperInsertions)/l_hyperInsertions for i in x])
    main_df['igBlast_ins_total_likelihood'] = main_df.apply(lambda row: np.prod(row['igBlast_ins_total_likelihood']) if row['igBlast_N_ins']>0 else 1, axis=1)
    
    main_df['igBlast_hyper_indels_total_likelihood'] = main_df.apply(lambda row: row['likelihood_N_of_hyperins']*row['likelihood_N_of_hyperdel']*row['igBlast_ins_total_likelihood']*row['igBlast_del_total_likelihood'], axis=1)
    
    #igBlast_gaps_expanded['likelihood_hyperindel_length'] = np.where(igBlast_gaps_expanded['hyperindel_type']=="ins",np.exp(-igBlast_gaps_expanded['hyperindel_length']/l_hyperInsertions)/l_hyperInsertions,np.exp(-igBlast_gaps_expanded['hyperindel_length']/l_hyperDeletions)/l_hyperDeletions)
    
    #igBlast_gaps_expanded_tmp = igBlast_gaps_expanded.groupby(['seq_id'], as_index=False)['likelihood_hyperindel_length'].prod()
    
    #igBlast_gaps_elaborated['likelihood_hyperindel_total_length'] = igBlast_gaps_elaborated['seq_id'].map(igBlast_gaps_expanded_tmp.set_index('seq_id')['likelihood_hyperindel_length'])
    #igBlast_gaps_elaborated['likelihood_hyperindel_total'] = igBlast_gaps_elaborated['likelihood_N_of_hyperins'] * igBlast_gaps_elaborated['likelihood_N_of_hyperdel'] * np.where(igBlast_gaps_elaborated['N_of_hyperindels']>0,igBlast_gaps_elaborated['likelihood_hyperindel_total_length'],1)
    
    #main_df['igBlast_HyperIndels_total_likelihood'] = main_df['bnab_ID'].map(igBlast_gaps_elaborated.set_index('seq_id')['likelihood_hyperindel_total'])
    
    #key_correlations = ['seq_id','V_best','D_best','J_best','stop_codon',\
    #                    'V_J_frame','productive','strand','CDR3','V_align_length',\
    #                    'V_align_length_beforeCDR3','N_of_mm','align_start_seq','align_end_seq','align_start_gene',\
    #                    'align_end_gene','igBlast_score','N_of_hyperindels','hyperindel_3prime_type','hyperindel_5prime_type',\
    #                    'inter_hyperindels_spacing','hyperindel_3prime_length','hyperindel_5prime_length'];
    
    #igBlast_gaps_correlations = pd.read_csv(wk_dir + "all_bnabs_seqs.gaps_correlations", sep='\s+', header=None, names=key_correlations, low_memory=False, keep_default_na=False);
    
    # inserire qui info da gaps correlations
    
    
    ##### IGoR seqs indexing #####
    
    #IGoR_indexing = pd.read_csv(wk_dir + sub_dir + "aligns/" + "indexed_sequences.csv", sep=";")
    #IGoR_indexing['seq_index'] = IGoR_indexing['seq_index'].astype(int)
    
    in_file = out_dir + "aligns/" + "indexed_sequences.csv"
    IGoR_noHyperIndels_indexing = pd.read_csv(in_file, sep=";")
    IGoR_noHyperIndels_indexing['seq_index'] = IGoR_noHyperIndels_indexing['seq_index'].astype(int)
    
    #main_df['IGoR_seq_index'] = main_df['orig_seq_nt'].map(IGoR_indexing.set_index('sequence')['seq_index'])
    main_df['IGoR_noHyperIndels_seq_index'] = main_df['noHyperIndels_trimmed_seq_nt'].map(IGoR_noHyperIndels_indexing.set_index('sequence')['seq_index'])
    
    ##### IGoR alignment #####
    
    #IGoR_V_align = pd.read_csv(sub_dir + "aligns/" + "V_alignments.csv", sep=";")
    #IGoR_V_align_best = IGoR_V_align.loc[IGoR_V_align.groupby(by="seq_index",sort=True)["score"].idxmax()];
    #IGoR_V_align_best['mismatches'] = IGoR_V_align_best['mismatches'].apply(expand_array);
    #IGoR_V_align_best['N_mismatches'] = IGoR_V_align_best['mismatches'].apply(len);
    #IGoR_V_align_best['frac_mismatches'] = IGoR_V_align_best['N_mismatches']/IGoR_V_align_best['length'];
    
    in_file = out_dir + "aligns/" + "V_alignments.csv"
    IGoR_noHyperIndels_V_align = pd.read_csv(in_file, sep=";")
    IGoR_noHyperIndels_V_align_best = IGoR_noHyperIndels_V_align.loc[IGoR_noHyperIndels_V_align.groupby(by="seq_index",sort=True)["score"].idxmax()];
    IGoR_noHyperIndels_V_align_best['mismatches'] = IGoR_noHyperIndels_V_align_best['mismatches'].apply(expand_array);
    IGoR_noHyperIndels_V_align_best['N_mismatches'] = IGoR_noHyperIndels_V_align_best['mismatches'].apply(len);
    IGoR_noHyperIndels_V_align_best['frac_mismatches'] = IGoR_noHyperIndels_V_align_best['N_mismatches']/IGoR_noHyperIndels_V_align_best['length'];
    
    #main_df['IGoR_V_best_align_identity'] = main_df['IGoR_seq_index'].map(IGoR_V_align_best.set_index('seq_index')['gene_name'])
    main_df['IGoR_noHyperIndels_V_best_align_identity'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_V_align_best.set_index('seq_index')['gene_name'])
    
    #main_df['IGoR_V_best_align_score'] = main_df['IGoR_seq_index'].map(IGoR_V_align_best.set_index('seq_index')['score'])
    main_df['IGoR_noHyperIndels_V_best_align_score'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_V_align_best.set_index('seq_index')['score'])
    
    #main_df['IGoR_V_best_align_offset'] = main_df['IGoR_seq_index'].map(IGoR_V_align_best.set_index('seq_index')['offset'])
    main_df['IGoR_noHyperIndels_V_best_align_offset'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_V_align_best.set_index('seq_index')['offset'])
    
    #main_df['IGoR_V_best_align_length'] = main_df['IGoR_seq_index'].map(IGoR_V_align_best.set_index('seq_index')['length'])
    main_df['IGoR_noHyperIndels_V_best_align_length'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_V_align_best.set_index('seq_index')['length'])
    
    #main_df['IGoR_V_best_align_mismatches'] = main_df['IGoR_seq_index'].map(IGoR_V_align_best.set_index('seq_index')['mismatches'])
    main_df['IGoR_noHyperIndels_V_best_align_mismatches'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_V_align_best.set_index('seq_index')['mismatches'])
    
    #main_df['IGoR_V_best_align_NofMismatches'] = main_df['IGoR_seq_index'].map(IGoR_V_align_best.set_index('seq_index')['N_mismatches'])
    main_df['IGoR_noHyperIndels_V_best_align_NofMismatches'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_V_align_best.set_index('seq_index')['N_mismatches'])
    
    #main_df['IGoR_V_best_align_fracMismatches'] = main_df['IGoR_seq_index'].map(IGoR_V_align_best.set_index('seq_index')['frac_mismatches'])
    main_df['IGoR_noHyperIndels_V_best_align_fracMismatches'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_V_align_best.set_index('seq_index')['frac_mismatches'])
  
    if config['chainType']=="HC":
      #IGoR_D_align = pd.read_csv(sub_dir + "aligns/" + "D_alignments.csv", sep=";")
      #IGoR_D_align_best = IGoR_D_align.loc[IGoR_D_align.groupby(by="seq_index",sort=True)["score"].idxmax()];
      #IGoR_D_align_best['mismatches'] = IGoR_D_align_best['mismatches'].apply(expand_array);
      #IGoR_D_align_best['N_mismatches'] = IGoR_D_align_best['mismatches'].apply(len);
      #IGoR_D_align_best['frac_mismatches'] = IGoR_D_align_best['N_mismatches']/IGoR_D_align_best['length'];
      
      in_file = out_dir + "aligns/" + "D_alignments.csv"
      IGoR_noHyperIndels_D_align = pd.read_csv(in_file, sep=";")
      IGoR_noHyperIndels_D_align_best = IGoR_noHyperIndels_D_align.loc[IGoR_noHyperIndels_D_align.groupby(by="seq_index",sort=True)["score"].idxmax()];
      IGoR_noHyperIndels_D_align_best['mismatches'] = IGoR_noHyperIndels_D_align_best['mismatches'].apply(expand_array);
      IGoR_noHyperIndels_D_align_best['N_mismatches'] = IGoR_noHyperIndels_D_align_best['mismatches'].apply(len);
      IGoR_noHyperIndels_D_align_best['frac_mismatches'] = IGoR_noHyperIndels_D_align_best['N_mismatches']/IGoR_noHyperIndels_D_align_best['length'];
      
      #main_df['IGoR_D_best_align_identity'] = main_df['IGoR_seq_index'].map(IGoR_D_align_best.set_index('seq_index')['gene_name'])
      main_df['IGoR_noHyperIndels_D_best_align_identity'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_D_align_best.set_index('seq_index')['gene_name'])
      
      #main_df['IGoR_D_best_align_score'] = main_df['IGoR_seq_index'].map(IGoR_D_align_best.set_index('seq_index')['score'])
      main_df['IGoR_noHyperIndels_D_best_align_score'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_D_align_best.set_index('seq_index')['score'])
      
      #main_df['IGoR_D_best_align_offset'] = main_df['IGoR_seq_index'].map(IGoR_D_align_best.set_index('seq_index')['offset'])
      main_df['IGoR_noHyperIndels_D_best_align_offset'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_D_align_best.set_index('seq_index')['offset'])
      
      #main_df['IGoR_D_best_align_length'] = main_df['IGoR_seq_index'].map(IGoR_D_align_best.set_index('seq_index')['length'])
      main_df['IGoR_noHyperIndels_D_best_align_length'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_D_align_best.set_index('seq_index')['length'])
      
      #main_df['IGoR_D_best_align_mismatches'] = main_df['IGoR_seq_index'].map(IGoR_D_align_best.set_index('seq_index')['mismatches'])
      main_df['IGoR_noHyperIndels_D_best_align_mismatches'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_D_align_best.set_index('seq_index')['mismatches'])
      
      #main_df['IGoR_D_best_align_NofMismatches'] = main_df['IGoR_seq_index'].map(IGoR_D_align_best.set_index('seq_index')['N_mismatches'])
      main_df['IGoR_noHyperIndels_D_best_align_NofMismatches'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_D_align_best.set_index('seq_index')['N_mismatches'])
      
      #main_df['IGoR_D_best_align_fracMismatches'] = main_df['IGoR_seq_index'].map(IGoR_D_align_best.set_index('seq_index')['frac_mismatches'])
      main_df['IGoR_noHyperIndels_D_best_align_fracMismatches'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_D_align_best.set_index('seq_index')['frac_mismatches'])
    
    #IGoR_J_align = pd.read_csv(sub_dir + "aligns/" + "J_alignments.csv", sep=";")
    #IGoR_J_align_best = IGoR_J_align.loc[IGoR_J_align.groupby(by="seq_index",sort=True)["score"].idxmax()];
    #IGoR_J_align_best['mismatches'] = IGoR_J_align_best['mismatches'].apply(expand_array);
    #IGoR_J_align_best['N_mismatches'] = IGoR_J_align_best['mismatches'].apply(len);
    #IGoR_J_align_best['frac_mismatches'] = IGoR_J_align_best['N_mismatches']/IGoR_J_align_best['length'];
    
    in_file = out_dir + "aligns/" + "J_alignments.csv"
    IGoR_noHyperIndels_J_align = pd.read_csv(in_file, sep=";")
    IGoR_noHyperIndels_J_align_best = IGoR_noHyperIndels_J_align.loc[IGoR_noHyperIndels_J_align.groupby(by="seq_index",sort=True)["score"].idxmax()];
    IGoR_noHyperIndels_J_align_best['mismatches'] = IGoR_noHyperIndels_J_align_best['mismatches'].apply(expand_array);
    IGoR_noHyperIndels_J_align_best['N_mismatches'] = IGoR_noHyperIndels_J_align_best['mismatches'].apply(len);
    IGoR_noHyperIndels_J_align_best['frac_mismatches'] = IGoR_noHyperIndels_J_align_best['N_mismatches']/IGoR_noHyperIndels_J_align_best['length'];
    
    #main_df['IGoR_J_best_align_identity'] = main_df['IGoR_seq_index'].map(IGoR_J_align_best.set_index('seq_index')['gene_name'])
    main_df['IGoR_noHyperIndels_J_best_align_identity'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_J_align_best.set_index('seq_index')['gene_name'])
    
    #main_df['IGoR_J_best_align_score'] = main_df['IGoR_seq_index'].map(IGoR_J_align_best.set_index('seq_index')['score'])
    main_df['IGoR_noHyperIndels_J_best_align_score'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_J_align_best.set_index('seq_index')['score'])
    
    #main_df['IGoR_J_best_align_offset'] = main_df['IGoR_seq_index'].map(IGoR_J_align_best.set_index('seq_index')['offset'])
    main_df['IGoR_noHyperIndels_J_best_align_offset'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_J_align_best.set_index('seq_index')['offset'])
    
    #main_df['IGoR_J_best_align_length'] = main_df['IGoR_seq_index'].map(IGoR_J_align_best.set_index('seq_index')['length'])
    main_df['IGoR_noHyperIndels_J_best_align_length'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_J_align_best.set_index('seq_index')['length'])
    
    #main_df['IGoR_J_best_align_mismatches'] = main_df['IGoR_seq_index'].map(IGoR_J_align_best.set_index('seq_index')['mismatches'])
    main_df['IGoR_noHyperIndels_J_best_align_mismatches'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_J_align_best.set_index('seq_index')['mismatches'])
    
    #main_df['IGoR_J_best_align_NofMismatches'] = main_df['IGoR_seq_index'].map(IGoR_J_align_best.set_index('seq_index')['N_mismatches'])
    main_df['IGoR_noHyperIndels_J_best_align_NofMismatches'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_J_align_best.set_index('seq_index')['N_mismatches'])
    
    #main_df['IGoR_J_best_align_fracMismatches'] = main_df['IGoR_seq_index'].map(IGoR_J_align_best.set_index('seq_index')['frac_mismatches'])
    main_df['IGoR_noHyperIndels_J_best_align_fracMismatches'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_J_align_best.set_index('seq_index')['frac_mismatches'])
    
    
    ##### IGoR output #####
    
    in_file = out_dir + "output/" + "Pgen_counts.csv"
    #IGoR_Pgen = pd.read_csv(wk_dir + "all_bnabs_output/" + "Pgen_counts.csv", sep=";")
    IGoR_noHyperIndels_Pgen = pd.read_csv(in_file, sep=";")
    
    #main_df['IGoR_Pgen'] = main_df['IGoR_seq_index'].map(IGoR_Pgen.set_index('seq_index')['Pgen_estimate'])
    main_df['IGoR_noHyperIndels_Pgen'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_Pgen.set_index('seq_index')['Pgen_estimate'])
    
    in_file = out_dir + "output/" + "best_scenarios_counts.csv"
    #IGoR_BestScen = pd.read_csv(wk_dir + "all_bnabs_output/" + "best_scenarios_counts.csv", sep=";")
    IGoR_noHyperIndels_BestScen = pd.read_csv(in_file, sep=";")
    
    #IGoR_BestScen = IGoR_BestScen.query('scenario_rank==1')
    IGoR_noHyperIndels_BestScen =IGoR_noHyperIndels_BestScen.query('scenario_rank==1')
    
    #IGoR_BestScen['mismatches'] = IGoR_BestScen['Mismatches'].apply(expand_array);
    #IGoR_BestScen['N_mismatches'] = IGoR_BestScen['mismatches'].apply(len);
    #IGoR_BestScen['N_mismatches'] = IGoR_BestScen['N_mismatches'].astype(int)
    
    IGoR_noHyperIndels_BestScen['mismatches'] = IGoR_noHyperIndels_BestScen['Mismatches'].apply(expand_array);
    IGoR_noHyperIndels_BestScen['N_mismatches'] = IGoR_noHyperIndels_BestScen['mismatches'].apply(len);
    IGoR_noHyperIndels_BestScen['N_mismatches'] = IGoR_noHyperIndels_BestScen['N_mismatches'].astype(int)
    
    #main_df['IGoR_bestScenario_condProba'] = main_df['IGoR_seq_index'].map(IGoR_BestScen.set_index('seq_index')['scenario_proba_cond_seq'])
    main_df['IGoR_noHyperIndels_bestScenario_condProba'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_BestScen.set_index('seq_index')['scenario_proba_cond_seq'])
    
    #main_df['IGoR_bestScenario_NofMismatches'] = main_df['IGoR_seq_index'].map(IGoR_BestScen.set_index('seq_index')['N_mismatches'])
    main_df['IGoR_noHyperIndels_bestScenario_NofMismatches'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_BestScen.set_index('seq_index')['N_mismatches'])
    
    ##### IGoR evaluation #####
    
    in_file = out_dir + "evaluate/" + "inference_logs.txt"
    #IGoR_evaluate = pd.read_csv(wk_dir + "all_bnabs_evaluate/" + "inference_logs.txt", sep=";")
    IGoR_noHyperIndels_evaluate = pd.read_csv(in_file, sep=";")
    
    #main_df['IGoR_seq_likelihood'] = main_df['IGoR_seq_index'].map(IGoR_evaluate.set_index('seq_index')['seq_likelihood'])
    main_df['IGoR_noHyperIndels_seq_likelihood'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_evaluate.set_index('seq_index')['seq_likelihood'])
    
    #main_df['IGoR_averageNofTotalMismatches'] = main_df['IGoR_seq_index'].map(IGoR_evaluate.set_index('seq_index')['seq_mean_n_errors'])
    main_df['IGoR_noHyperIndels_averageNofTotalMismatches'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_evaluate.set_index('seq_index')['seq_mean_n_errors'])
    
    #main_df['IGoR_bestScenario_likelihood'] = main_df['IGoR_seq_index'].map(IGoR_evaluate.set_index('seq_index')['seq_best_scenario'])
    main_df['IGoR_noHyperIndels_bestScenario_likelihood'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_evaluate.set_index('seq_index')['seq_best_scenario'])
    
    #main_df['IGoR_NofScenarios'] = main_df['IGoR_seq_index'].map(IGoR_evaluate.set_index('seq_index')['seq_n_scenarios'])
    main_df['IGoR_noHyperIndels_NofScenarios'] = main_df['IGoR_noHyperIndels_seq_index'].map(IGoR_noHyperIndels_evaluate.set_index('seq_index')['seq_n_scenarios'])
    
    #main_df['IGoR_hyperMutations_likelihood'] = main_df['IGoR_seq_likelihood'] / main_df['IGoR_Pgen']
    main_df['IGoR_noHyperIndels_hyperMutations_likelihood'] = main_df['IGoR_noHyperIndels_seq_likelihood'] / main_df['IGoR_noHyperIndels_Pgen']
    
    ##### print for check #####
    
    #with pd.option_context('display.max_rows', 999):
    #        print(main_df)
    
    ##### export the dataframe #####
    
    out_file = out_dir + config['input_file_prefix'] + "_" + out_dir.split('/')[-2] + '.df'
    main_df.to_csv(out_file, index=False, sep=';')
    
    # old name of this sub-folder: IGoR_analysis
    if not os.path.isdir(wk_dir + 'igor_bnabs_summary'):
      os.mkdir(wk_dir + 'igor_bnabs_summary')
    
    out_file = wk_dir + 'igor_bnabs_summary/' + config['input_file_prefix'] + "_" + out_dir.split('/')[-2] + '.df'
    main_df.to_csv(out_file, index=False, sep=';')
