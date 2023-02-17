#!/usr/bin/python3

# coding: utf-8

# Filename: launch_annotation.py
# Author: Cosimo Lupo, <cosimo.lupo89@gmail.com>
# Last updated: 10 February 2023

import numpy as np
import pandas as pd
import glob
import os
import sys
import shutil
import multiprocessing as mp
pd.set_option('display.max_columns',100)
from scipy.stats import poisson
from scipy.stats import binom
from scipy.optimize import curve_fit
from scipy.optimize import minimize
#from pathlib import Path
import Bio
from Bio import SeqIO
import datetime
import re
import yaml

#sys.path.insert(1, '/home/lupo/')

from aux_funcs_annotation import natural_keys, make_csv_from_fasta, make_fasta_from_csv
from aux_funcs_annotation import run_igBlast, parse_igBlast, revert_seq, reconstruct_gapped_seqs, write_anchored_seqs

########## settings ##########

with open('config_annotation.yaml') as config_file:
  config = yaml.full_load(config_file)

wk_dir = os.path.abspath(config['wk_dir'])
if wk_dir[-1] != '/':
  wk_dir += '/'
input_dir = os.path.abspath(config['input_dir'])
if input_dir[-1] != '/':
  input_dir += '/'

cohort_dir = input_dir + config['cohort'] + '_' + config['cellType'] + '/'
sample_dirs = glob.glob(cohort_dir + 'BZR*')

chainType_dict = {"heavy": "HC", "kappa": "KC", "lambda": "LC"}
chainType_shortHand = chainType_dict[config['chainType']]

########## main ##########

print('\n' + '\033[1m' + " ***** Pre-processing, annotation through igBlast and post-processing ***** " + '\033[0m' + '\n')

# Step 1: check which format seqs have to be read and (if needed) produce the other format

if config['preProcessFiles']:

  print("\nStep 1: check which format seqs have to be read and (if needed) produce the other format")

  sample_dirs.sort(key=natural_keys)
  for i,sample_dir in enumerate(sample_dirs):

    sample = sample_dir.split('/')[-1]
    print("    sample: " + sample + "  [" + config['cohort'] + " " + config['cellType'] + " " + config['chainType'] + "]")
    filenames = []

    # check csv
    csv_files = glob.glob(sample_dir + '/' + '*' + chainType_shortHand + '.csv')
    if(len(csv_files)>0):
      for j,fullfilename in enumerate(csv_files):
        filename = fullfilename.split('/')[-1].split('.')[0]
        filenames.append(filename)

    # check fasta
    fasta_files = glob.glob(sample_dir + '/' + '*' + chainType_shortHand + '.fasta')
    if(len(fasta_files)>0):
      for j,fullfilename in enumerate(fasta_files):
        filename = fullfilename.split('/')[-1].split('.')[0]
        filenames.append(filename)

    # make unique
    filenames = list(set(filenames))
    filenames.sort(key=natural_keys)

    # produce the missing format (if any)
    for i,filename in enumerate(filenames):
      fullfilename = sample_dir + '/' + filename

      if(os.path.isfile(fullfilename + '.csv')==True and os.path.isfile(fullfilename + '.fasta')==False):
        try:
          make_fasta_from_csv(fullfilename + '.csv', headers=True, sep=';')
        except BaseException as err:
          print(err)

      elif(os.path.isfile(fullfilename + '.csv')==False and os.path.isfile(fullfilename + '.fasta')==True):
        try:
          make_csv_from_fasta(fullfilename + '.fasta', headers=['seq_ID','raw_seq_nt'], sep=';')
        except BaseException as err:
          print(err)

# Step 2: run igBlast

if config['runIgBlast']:

  print("\nStep 2: run igBlast")

  sample_dirs.sort(key=natural_keys)
  for i,sample_dir in enumerate(sample_dirs):
    sample = sample_dir.split('/')[-1]
    print("    sample: " + sample + "  [" + config['cohort'] + " " + config['cellType'] + " " + config['chainType'] + "]")
    fasta_files = glob.glob(sample_dir + '/' + '*' + chainType_shortHand + '.fasta')
    if(len(fasta_files)>0):
      fasta_files.sort(key=natural_keys)
      for j,fullfilename in enumerate(fasta_files):
        in_file = fullfilename.split('.')[0] + '.fasta'
        out_file = fullfilename.split('.')[0] + '.igBlast_raw_output'

        t1 = datetime.datetime.now()
        try:
          run_igBlast(in_file, config['species'], config['chainType'])
        except BaseException as err:
          raise err
        t2 = datetime.datetime.now()
        print('        igBlast running time:', t2-t1)

# Step 3: parse the raw output from igBlast

if config['parseIgBlast']:

  print("\nStep 3: parse the raw output from igBlast")

  sample_dirs.sort(key=natural_keys)
  for i,sample_dir in enumerate(sample_dirs):
    sample = sample_dir.split('/')[-1]
    print("    sample: " + sample + "  [" + config['cohort'] + " " + config['cellType'] + " " + config['chainType'] + "]")
    igBlast_raw_output_files = glob.glob(sample_dir + '/' + '*' + chainType_shortHand + '.igBlast_raw_output')
    if(len(igBlast_raw_output_files)>0):
      igBlast_raw_output_files.sort(key=natural_keys)
      for j,fullfilename in enumerate(igBlast_raw_output_files):
        in_file = fullfilename.split('.')[0] + '.igBlast_raw_output'
        out_file = fullfilename.split('.')[0] + '.igBlast_statistics'

        t1 = datetime.datetime.now()

        try:
          df = parse_igBlast(in_file, config['chainType'], requireJ=True)
        except BaseException as err:
          print(err)

        # UMI counts
        df['UMI_counts'] = df['seq_ID'].apply(lambda x: int(x.split(':')[3][1:]))

        # Quality filtering
        df = df[df['V_best_align_length_beforeCDR3']>=config['V_min_len']]    # V gene should align at least for V_min_len nt
        df = df[df['strand']=="+"]    # Only sequences read in the correct direction
        #df = df[df['UMI_counts']>config['UMI_thr']]    # Only UMIs counted at least UMI_thr times

        # mapping of nt sequences through IDs (included primers and C segment, if any)
        seqs_file = fullfilename.split('.')[0] + '.csv'
        raw_seqs = pd.read_csv(seqs_file, sep=';', low_memory=False, keep_default_na=True, index_col=False)
        raw_seqs['seq_ID'] = raw_seqs['seq_ID'].astype('str')
        df['raw_seq_nt'] = df['seq_ID'].map(raw_seqs.set_index('seq_ID')['raw_seq_nt'])
        df['seq_nt'] = df.apply(lambda row: row['raw_seq_nt'][row['V_best_align_start_seq']-1:row['J_best_align_end_seq']], axis=1)
        df['seq_nt_len'] = df['seq_nt'].apply(lambda x: len(x))

        # Clonal abundance
        counts = df.groupby('seq_nt').count().reset_index()
        counts = counts[['seq_nt','seq_ID']]
        counts.rename(columns={'seq_ID': 'seq_nt_counts'}, inplace=True)
        df['seq_nt_counts'] = df['seq_nt'].map(counts.set_index('seq_nt')['seq_nt_counts'])
        df = df.drop_duplicates(subset=['seq_nt'],keep='first')    # Drop seq_nt duplicates

        # Reconstruct gapped query and germlines (useful for Natanael's trees)
        '''
        if(config['chainType']=="heavy"):
          Vdatabase = blast_database + "BCR_Heavy/forIgBlast_IGHV_Homo_Sapiens_F.fasta"
          Ddatabase = blast_database + "BCR_Heavy/forIgBlast_IGHD_Homo_Sapiens_F.fasta"
          Jdatabase = blast_database + "BCR_Heavy/forIgBlast_IGHJ_Homo_Sapiens_F.fasta"
        elif(config['chainType']=="kappa"):
          Vdatabase = blast_database + "BCR_Kappa/forIgBlast_IGKV_Homo_Sapiens_F.fasta"
          Ddatabase = blast_database + "BCR_Heavy/forIgBlast_IGHD_Homo_Sapiens_F.fasta"
          Jdatabase = blast_database + "BCR_Kappa/forIgBlast_IGKJ_Homo_Sapiens_F.fasta"
        elif(config['chainType']=="lambda"):
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
        '''
        df['temp'] = df.apply(lambda row: reconstruct_gapped_seqs(row), axis=1)
        df['gapped_query'] = df['temp'].apply(lambda x: x.split(';')[0])
        df['gapped_germline'] = df['temp'].apply(lambda x: x.split(';')[1])

        # Revert indels in the sequence
        df['reverted_seq_nt'] = df.apply(lambda row: revert_seq(row), axis=1)

        # Drop unnecessary data
        df = df.drop(['raw_seq_nt', 'V_best_aligned_query', 'V_best_aligned_germline', \
                      'J_best_aligned_query', 'J_best_aligned_germline', 'temp'], axis=1)

        # Export on file
        df.to_csv(out_file, index=False, sep=';')

        t2 = datetime.datetime.now()
        print('        parsing & formatting time:', t2-t1)

# Step 4: sort NP-P into new files and gather cohortwide data

if config['sortAndCohort']:

  print("\nStep 4: sort NP-P into new files and gather cohortwide data")


  # Check if the cohort-wide already exists (otherwise create it)
  if not os.path.exists(cohort_dir + 'cohortWide_analysis'):
    os.makedirs(cohort_dir + 'cohortWide_analysis')


  ##### gathering together initial sequences #####

  df_cohort = pd.DataFrame();

  sample_dirs.sort(key=natural_keys)
  for i,sample_dir in enumerate(sample_dirs):

    sample = sample_dir.split('/')[-1]

    seqs_files = glob.glob(sample_dir + '/' + '*' + chainType_shortHand + '.csv')
    seqs_files.sort(key=natural_keys)

    if(len(seqs_files)>0):
      for j,fullfilename in enumerate(seqs_files):
        in_file = fullfilename.split('.')[0] + '.csv'
        df = pd.read_csv(in_file, sep=';', low_memory=False, keep_default_na=True, index_col=False)

        # Append to the cohort dataframe
        df_cohort = pd.concat([df_cohort,df]).reset_index(drop=True)

  # df_cohort headers: seq_ID ; raw_seq_nt
  # df_cohort contains some raw_seq_nt duplicates coming from different individuals (but they are very few)
  # So we drop duplicates, taking the unique raw_seq_nt from the first individual in the list
  df_cohort = df_cohort.drop_duplicates(subset=['raw_seq_nt'], keep='first')

  # Write the cohort-wide .csv file with all the sequences together
  df_cohort.to_csv(cohort_dir + 'cohortWide_analysis/' + config['cohort'] + '_' + config['cellType'] + '_' + chainType_shortHand + '.csv', index=False, sep=';')


  ##### gathering together .igBlast_statistics files and then sorting sequences #####

  df_cohort = pd.DataFrame();

  sample_dirs.sort(key=natural_keys)
  for i,sample_dir in enumerate(sample_dirs):

    sample = sample_dir.split('/')[-1]
    print("    sample: " + sample + "  [" + config['cohort'] + " " + config['cellType'] + " " + config['chainType'] + ", ", end="")

    igBlast_statistics_files = glob.glob(sample_dir + '/' + '*' + chainType_shortHand + '.igBlast_statistics')
    igBlast_statistics_files.sort(key=natural_keys)

    if(len(igBlast_statistics_files)>0):
      for j,fullfilename in enumerate(igBlast_statistics_files):
        in_file = fullfilename.split('.')[0] + '.igBlast_statistics'

        # Open .igBlast_statistics files
        df = pd.read_csv(in_file, sep=';', low_memory=False, keep_default_na=True, index_col=False)
        df['cohort'] = config['cohort']
        df['cellType'] = config['cellType']
        df['sample'] = sample

        # Drop duplicate sequences
        df = df.drop_duplicates(subset=['seq_nt'],keep='first').reset_index(drop=True)    # Drop seq_nt duplicates

        # Trim sequences at the two edges
        df['trimmed_seq_nt'] = df['seq_nt'].apply(lambda x: x[config['n_l']:-config['n_r']])

        # Keep only sequences above UMI count threshold
        if config['UMIfilter']:
          df = df[df['UMI_counts']>=config['UMI_thr']]

        # Keep only sequences with annotated CDR3
        if config['onlyAnnotatedCDR3']:
          df = df.query('CDR3_nt==CDR3_nt')

        # Append to the cohort dataframe
        df_cohort = pd.concat([df_cohort,df]).reset_index(drop=True)

        print("data_size: " + str(len(df)) + "]")

        for filterOutHyperIndels in [False,True]:
          for trimEdges in [False,True]:

            # Filter P/NP sequences
            # Here I'm using the following definition:
            # - P:  Both V,J are in-frame and no stop codon in the CDR3
            # - NP: Either V,J are out-frame or there is a stop codon in the CDR3
            #df_NP = df.loc[(df['V_J_frame']=="Out-of-frame") | (df['stop_codon_CDR3']=="Yes")]
            df_NP = df[df['V_J_frame']=="Out-of-frame"]
            df_P = df[(df['V_J_frame']=="In-frame") & (df['stop_codon_CDR3']=="No")]

            # Export sequences on file
            out_file_NP = fullfilename.split('.')[0] + '_uniqueSeqs_NP'
            out_file_P = fullfilename.split('.')[0] + '_uniqueSeqs_P'
            if config['UMIfilter']:
              out_file_NP += '_UMIthr_' + str(config['UMI_thr'])
              out_file_P += '_UMIthr_' + str(config['UMI_thr'])
            if(trimEdges):
              out_file_NP += '_trimmed_' + str(config['n_l']) + '_' + str(config['n_r'])
              out_file_P += '_trimmed_' + str(config['n_l']) + '_' + str(config['n_r'])
              keys2 = ['seq_ID','trimmed_seq_nt','V_best_identity','V_best_align_start_seq','V_best_align_end_seq','V_best_align_start_gene','V_best_align_end_gene','N_indels']
              df_NP = df_NP[keys2]
              df_P = df_P[keys2]
              df_NP = df_NP.rename(columns={'trimmed_seq_nt': 'seq_nt'})
              df_P = df_P.rename(columns={'trimmed_seq_nt': 'seq_nt'})
            else:
              keys2 = ['seq_ID','seq_nt','V_best_identity','V_best_align_start_seq','V_best_align_end_seq','V_best_align_start_gene','V_best_align_end_gene','N_indels']
              df_NP = df_NP[keys2]
              df_P = df_P[keys2]
            if(filterOutHyperIndels):
              out_file_NP += '_noHyperIndels'
              out_file_P += '_noHyperIndels'
              df_NP = df_NP.query('N_indels==0')
              df_P = df_P.query('N_indels==0')

            # Write csv
            df_NP[['seq_ID','seq_nt']].to_csv(out_file_NP + '.csv', index=False, sep=';')
            df_P[['seq_ID','seq_nt']].to_csv(out_file_P + '.csv', index=False, sep=';')

            # Write fasta
            ofile = open(out_file_NP + '.fasta', "w")
            for ii in range(len(df_NP)):
              ofile.write(">" + df_NP.iloc[ii,0] + "\n" + df_NP.iloc[ii,1] + "\n")
            ofile.close()
            ofile = open(out_file_P + '.fasta', "w")
            for ii in range(len(df_P)):
              ofile.write(">" + df_P.iloc[ii,0] + "\n" + df_P.iloc[ii,1] + "\n")
            ofile.close()

            if(filterOutHyperIndels==False and trimEdges==False):

              # Formatted csv for indels software
              # When numbering starting from 0:
              #  - initial position is given by 'start' index - 1
              #  - final position (not included) is given by 'end' index, so that last included is given by 'end' index - 1

              anchored_subset_size = 100000

              for produc in ['NP','P']:

                out_file = (out_file_NP if (produc == 'NP') else out_file_P)
                out_file += '_anchored.csv'
                out_f = open(out_file, "w")
                out_f.write("seq_ID;aligned_seq_nt;V_best;V_best_start;V_best_end\n")
                df2 = (df_NP if (produc == 'NP') else df_P)
                for r,row in df2.iterrows():
                  out_f.write(write_anchored_seqs(row) + "\n")
                out_f.close()

                if(len(df2)>anchored_subset_size):
                  out_file = (out_file_NP if (produc == 'NP') else out_file_P)
                  out_file += '_anchored_subset100K.csv'
                  out_f = open(out_file, "w")
                  out_f.write("seq_ID;aligned_seq_nt;V_best;V_best_start;V_best_end\n")
                  df2 = df2.sample(n=anchored_subset_size).sort_index().reset_index(drop=True)
                  for r,row in df2.iterrows():
                    out_f.write(write_anchored_seqs(row) + "\n")
                  out_f.close()

  # Export the cohort-wide .csv and .igBlast_statistics files
  df_cohort[['seq_ID','seq_nt']].to_csv(cohort_dir + 'cohortWide_analysis/' + config['cohort'] + '_' + config['cellType'] + '_' + chainType_shortHand + '_uniqueSeqs.csv', index=False, sep=';')
  df_cohort.to_csv(cohort_dir + 'cohortWide_analysis/' + config['cohort'] + '_' + config['cellType'] + '_' + chainType_shortHand + '_uniqueSeqs.igBlast_statistics', index=False, sep=';')

  for filterOutHyperIndels in [False,True]:
    for trimEdges in [False,True]:

      # Filter P/NP sequences
      # Here I'm using the following definition:
      # - P:  Both V,J are in-frame and no stop codon in the CDR3
      # - NP: Either V,J are out-frame or there is a stop codon in the CDR3
      #df_cohort_NP = df_cohort.loc[(df_cohort['V_J_frame']=="Out-of-frame") | (df_cohort['stop_codon_CDR3']=="Yes")]
      df_cohort_NP = df_cohort[df_cohort['V_J_frame']=="Out-of-frame"]
      df_cohort_P = df_cohort[(df_cohort['V_J_frame']=="In-frame") & (df_cohort['stop_codon_CDR3']=="No")]

      # Export sequences on file
      out_file = cohort_dir + 'cohortWide_analysis/' + config['cohort'] + '_' + config['cellType'] + '_' + chainType_shortHand
      out_file_NP = out_file + '_uniqueSeqs_NP'
      out_file_P = out_file + '_uniqueSeqs_P'
      if config['UMIfilter']:
        out_file_NP += '_UMIthr_' + str(config['UMI_thr'])
        out_file_P += '_UMIthr_' + str(config['UMI_thr'])
      if(trimEdges):
        out_file_NP += '_trimmed_' + str(config['n_l']) + '_' + str(config['n_r'])
        out_file_P += '_trimmed_' + str(config['n_l']) + '_' + str(config['n_r'])
        keys2 = ['seq_ID','trimmed_seq_nt','V_best_identity','V_best_align_start_seq','V_best_align_end_seq','V_best_align_start_gene','V_best_align_end_gene','N_indels']
        df_cohort_NP = df_cohort_NP[keys2]
        df_cohort_P = df_cohort_P[keys2]
        df_cohort_NP = df_cohort_NP.rename(columns={'trimmed_seq_nt': 'seq_nt'})
        df_cohort_P = df_cohort_P.rename(columns={'trimmed_seq_nt': 'seq_nt'})
      else:
        keys2 = ['seq_ID','seq_nt','V_best_identity','V_best_align_start_seq','V_best_align_end_seq','V_best_align_start_gene','V_best_align_end_gene','N_indels']
        df_cohort_NP = df_cohort_NP[keys2]
        df_cohort_P = df_cohort_P[keys2]
      if(filterOutHyperIndels):
        out_file_NP += '_noHyperIndels'
        out_file_P += '_noHyperIndels'
        df_cohort = df_cohort.query('N_indels==0')

      # Write csv
      df_cohort_NP.to_csv(out_file_NP + '.csv', index=False, sep=';')
      df_cohort_P.to_csv(out_file_P + '.csv', index=False, sep=';')

      # Write fasta
      ofile = open(out_file_NP + '.fasta', "w")
      for ii in range(len(df_cohort_NP)):
        ofile.write(">" + df_cohort_NP.iloc[ii,0] + "\n" + df_cohort_NP.iloc[ii,1] + "\n")
      ofile.close()
      ofile = open(out_file_P + '.fasta', "w")
      for ii in range(len(df_cohort_P)):
        ofile.write(">" + df_cohort_P.iloc[ii,0] + "\n" + df_cohort_P.iloc[ii,1] + "\n")
      ofile.close()

      if(filterOutHyperIndels==False and trimEdges==False):

        # Formatted csv for indels software
        # When numbering starting from 0:
        #  - initial position is given by 'start' index - 1
        #  - final position (not included) is given by 'end' index, so that last included is given by 'end' index - 1

        anchored_subset_size = 100000

        for produc in ['NP','P']:

          out_file = (out_file_NP if (produc == 'NP') else out_file_P)
          out_file += '_anchored.csv'
          out_f = open(out_file, "w")
          out_f.write("seq_ID;aligned_seq_nt;V_best;V_best_start;V_best_end\n")
          df2 = (df_cohort_NP if (produc == 'NP') else df_cohort_P)
          for r,row in df2.iterrows():
            out_f.write(write_anchored_seqs(row) + "\n")
          out_f.close()

          if(len(df2)>anchored_subset_size):
            out_file = (out_file_NP if (produc == 'NP') else out_file_P)
            out_file += '_anchored_subset100K.csv'
            out_f = open(out_file, "w")
            out_f.write("seq_ID;aligned_seq_nt;V_best;V_best_start;V_best_end\n")
            df2 = df2.sample(n=anchored_subset_size).sort_index().reset_index(drop=True)
            for r,row in df2.iterrows():
              out_f.write(write_anchored_seqs(row) + "\n")
            out_f.close()

quit()
