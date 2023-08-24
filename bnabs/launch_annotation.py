#!/usr/bin/python3

# coding: utf-8

# Filename: launch_annotation.py
# Author: Cosimo Lupo, <cosimo.lupo89@gmail.com>
# Last updated: 10 February 2023

import numpy as np
import pandas as pd
import glob as glob
import os as os
import shutil
import multiprocessing as mp
pd.set_option('display.max_columns',100)
import Bio
from Bio import SeqIO
import datetime
import re
import yaml

#sys.path.insert(1, '../libraries/')
from aux_funcs_annotation import natural_keys, make_csv_from_fasta, make_fasta_from_csv
from aux_funcs_annotation import run_igBlast, parse_igBlast, revert_seq, reconstruct_gapped_seqs, write_anchored_seqs

########## settings ##########

with open('config_annotation.yaml') as config_file:
  config = yaml.full_load(config_file)

wkdir = os.path.abspath(config['wkdir'])
if wkdir[-1] != '/':
  wkdir += '/'
inputdir = os.path.abspath(config['inputdir'])
if inputdir[-1] != '/':
  inputdir += '/'

chainType_dict = {"heavy": "HC", "kappa": "KC", "lambda": "LC"}
chainType_shortHand = chainType_dict[config['chainType']]

########## main ##########

print('\n' + '\033[1m' + " ***** Pre-processing, annotation through igBlast and post-processing ***** " + '\033[0m' + '\n')

# Step 1: check which format seqs have to be read and (if needed) produce the other format

if config['preProcessFiles']:

  print("\nStep 1: check which format seqs have to be read and (if needed) produce the other format")

  filenames = []

  # check csv
  csv_files = glob.glob(inputdir + config['input_file_prefix'] + '_' + chainType_shortHand + '.csv')
  if(len(csv_files)>0):
    for j,fullfilename in enumerate(csv_files):
      filename = fullfilename.split('/')[-1].split('.')[0]
      filenames.append(filename)

  # check fasta
  fasta_files = glob.glob(inputdir + config['input_file_prefix'] + '_' + chainType_shortHand + '.fasta')
  if(len(fasta_files)>0):
    for j,fullfilename in enumerate(fasta_files):
      filename = fullfilename.split('/')[-1].split('.')[0]
      filenames.append(filename)

  # make unique
  filenames = list(set(filenames))
  filenames.sort(key=natural_keys)

  # produce the missing format (if any)
  for i,filename in enumerate(filenames):
    fullfilename = inputdir + filename

    if(os.path.isfile(fullfilename + '.csv')==True and os.path.isfile(fullfilename + '.fasta')==False):
      try:
        make_fasta_from_csv(fullfilename + '.csv', headers=True, sep=';')
      except BaseException as err:
        raise err

    elif(os.path.isfile(fullfilename + '.csv')==False and os.path.isfile(fullfilename + '.fasta')==True):
      try:
        make_csv_from_fasta(fullfilename + '.fasta', headers=['bnab_ID','raw_seq_nt'], sep=';')
      except BaseException as err:
        raise err

# Step 2: run igBlast

if config['runIgBlast']:

  print("\nStep 2: run igBlast")

  fasta_files = glob.glob(inputdir + config['input_file_prefix'] + '_' + chainType_shortHand + '.fasta')

  if(len(fasta_files)>0):
    fasta_files.sort(key=natural_keys)
    for i,fullfilename in enumerate(fasta_files):
      in_file = fullfilename.split('.')[0] + '.fasta'
      out_file = fullfilename.split('.')[0] + '.igBlast_raw_output'

      t1 = datetime.datetime.now()
      try:
        run_igBlast(in_file, config['species'], config['chainType'])
      except BaseException as err:
        raise ValueError(err)
      t2 = datetime.datetime.now()
      print('        igBlast running time:', t2-t1)

# Step 3: parse the raw output from igBlast

if config['parseIgBlast']:

  print("\nStep 3: parse the raw output from igBlast")

  igBlast_raw_output_files = glob.glob(inputdir + config['input_file_prefix'] + '_' + chainType_shortHand + '.igBlast_raw_output')

  if(len(igBlast_raw_output_files)>0):
    igBlast_raw_output_files.sort(key=natural_keys)
    for i,fullfilename in enumerate(igBlast_raw_output_files):
      in_file = fullfilename.split('.')[0] + '.igBlast_raw_output'
      out_file = fullfilename.split('.')[0] + '.igBlast_statistics'

      t1 = datetime.datetime.now()

      try:
        df = parse_igBlast(in_file, config['chainType'], requireJ=True)
      except BaseException as err:
        raise err

      # Quality filtering
      #df = df[df['V_best_align_length_beforeCDR3']>=config['V_min_len']]    # V gene should align at least for V_min_len nt
      df = df[df['strand']=="+"]    # Only sequences read in the correct direction

      # mapping of nt sequences through IDs (included primers and C segment, if any)
      seqs_file = fullfilename.split('.')[0] + '.csv'
      raw_seqs = pd.read_csv(seqs_file, sep=';', low_memory=False, keep_default_na=True, index_col=False)
      if config['chainType']=="heavy":
        raw_seqs = raw_seqs.rename(columns={"heavy_seq_nt": "raw_seq_nt"})
      elif config['chainType']=="kappa" or config['chainType']=="lambda":
        raw_seqs = raw_seqs.rename(columns={"light_seq_nt": "raw_seq_nt"})
      raw_seqs['bnab_ID'] = raw_seqs['bnab_ID'].astype('str')
      df = df.rename(columns={"seq_ID": "bnab_ID"})
      df['raw_seq_nt'] = df['bnab_ID'].map(raw_seqs.set_index('bnab_ID')['raw_seq_nt'])
      df['seq_nt'] = df.apply(lambda row: row['raw_seq_nt'][row['V_best_align_start_seq']-1:row['J_best_align_end_seq']], axis=1)
      df['seq_nt_len'] = df['seq_nt'].apply(lambda x: len(x))

      # Reconstruct gapped query and germlines (useful for Natanael's trees)
      blast_database = config['blast_database'] + config['species'] + "/"
      if config['chainType']=="heavy":
        Vdatabase = blast_database + "BCR_Heavy/forIgBlast_IGHV_Homo_Sapiens_F.fasta"
        Ddatabase = blast_database + "BCR_Heavy/forIgBlast_IGHD_Homo_Sapiens_F.fasta"
        Jdatabase = blast_database + "BCR_Heavy/forIgBlast_IGHJ_Homo_Sapiens_F.fasta"
      elif config['chainType']=="kappa":
        Vdatabase = blast_database + "BCR_Kappa/forIgBlast_IGKV_Homo_Sapiens_F.fasta"
        Ddatabase = blast_database + "BCR_Heavy/forIgBlast_IGHD_Homo_Sapiens_F.fasta"
        Jdatabase = blast_database + "BCR_Kappa/forIgBlast_IGKJ_Homo_Sapiens_F.fasta"
      elif config['chainType']=="lambda":
        Vdatabase = blast_database + "BCR_Lambda/forIgBlast_IGLV_Homo_Sapiens_F.fasta"
        Ddatabase = blast_database + "BCR_Heavy/forIgBlast_IGHD_Homo_Sapiens_F.fasta"
        Jdatabase = blast_database + "BCR_Lambda/forIgBlast_IGLJ_Homo_Sapiens_F.fasta"
      V_dict = {}
      with open(Vdatabase, 'r') as temp_f:
        for line in temp_f:
          line = line.strip()
          if(line!=''):
            if(line[0]=='>'):
              V_name = line[1:]
            else:
              V_string = line.upper()
              V_dict[V_name] = V_string
      df['temp'] = df.apply(lambda row: reconstruct_gapped_seqs(row), axis=1)
      df['gapped_query'] = df['temp'].apply(lambda x: x.split(';')[0])
      df['gapped_germline'] = df['temp'].apply(lambda x: x.split(';')[1])

      # Revert indels in the sequence
      df['reverted_seq_nt'] = df.apply(lambda row: revert_seq(row), axis=1)

      # Drop unnecessary data
      df = df.drop(['raw_seq_nt', 'V_best_aligned_query', 'V_best_aligned_germline', \
                    'J_best_aligned_query', 'J_best_aligned_germline','temp'], axis=1)

      # Export on file
      df.to_csv(out_file, index=False, sep=';')

      t2 = datetime.datetime.now()
      print('        parsing & formatting time:', t2-t1)

# Step 4: sort NP-P into new files and gather cohortwide data

if config['produceFinalFiles']:

  print("\nStep 4: sort NP-P into new files and gather cohortwide data")

  igBlast_statistics_files = glob.glob(inputdir + config['input_file_prefix'] + '_' + chainType_shortHand + '.igBlast_statistics')

  if(len(igBlast_statistics_files)>0):
    igBlast_statistics_files.sort(key=natural_keys)
    for i,fullfilename in enumerate(igBlast_statistics_files):
      in_file = fullfilename.split('.')[0] + '.igBlast_statistics'

      # Open .igBlast_statistics files
      df = pd.read_csv(in_file, sep=';', low_memory=False, keep_default_na=True, index_col=False)
      df['cellType'] = config['cellType']

      # Trim sequences at the two edges
      df['seq_nt_trimmed'] = df['seq_nt'].apply(lambda x: x[config['n_l']:-config['n_r']])
      df['reverted_seq_nt_trimmed'] = df['reverted_seq_nt'].apply(lambda x: x[config['n_l']:-config['n_r']])

      #df = df.drop_duplicates(subset=['reverted_seq_nt'],keep='first').reset_index(drop=True)    # Drop seq_nt duplicates
      #if config['onlyAnnotatedCDR3']:
      #  df = df[df['CDR3_nt']==df['CDR3_nt']]

      out_file = "/".join(fullfilename.split('/')[:-1]) + '/' + fullfilename.split('/')[-1].split('.')[0] + '_noHyperIndels.csv'
      df[['bnab_ID','reverted_seq_nt']].to_csv(out_file, index=False, sep=';')
      try:
        make_fasta_from_csv(out_file.split('.')[0] + '.csv', headers=True, sep=';')
      except BaseException as err:
        raise err

      out_file = "/".join(fullfilename.split('/')[:-1]) + '/' + fullfilename.split('/')[-1].split('.')[0] + '_noHyperIndels_trimmed_' + str(config['n_l']) + '_' + str(config['n_r']) + '.csv'
      df[['bnab_ID','reverted_seq_nt_trimmed']].to_csv(out_file, index=False, sep=';')
      try:
        make_fasta_from_csv(out_file.split('.')[0] + '.csv', headers=True, sep=';')
      except BaseException as err:
        raise err
