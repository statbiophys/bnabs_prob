#!/usr/bin/python3

# coding: utf-8

# Filename: compute_stats.py
# Author: Cosimo Lupo, <cosimo.lupo89@gmail.com>
# Last updated: 14 July 2023

import pandas as pd
import numpy as np
import glob
import os
import random

produc = "P"
SHMmodel = "add5mer"
cohorts = ["healthy_control","hiv1","hcv",\
           "hiv1_non_neutralizers","hiv1_weak_neutralizers",\
           "hiv1_intermediate_neutralizers","hiv1_top_neutralizers",\
           "hiv1_art_off","hiv1_art_on"]
T_bootstrap = 500

source_ig_data = "./"

if not os.path.isdir(source_ig_data + "stats"):
    os.mkdir(source_ig_data + "stats")

print()

print('%s sequences, %s SHM model' % (produc, SHMmodel))

for chain_type in ["HC","KC","LC"]:
    
    print(' > %s chain' % chain_type)
    
    if chain_type=="HC":
        dx = 1
        x_Pgen = np.arange(-60,0+0.5*dx,dx)
        dx = 3
        x_Pshm = np.arange(-120,0+0.5*dx,dx)
    elif chain_type=="KC" or chain_type=="LC":
        dx = 0.75
        x_Pgen = np.arange(-20,0+0.5*dx,dx)
        dx = 2
        x_Pshm = np.arange(-100,0+0.5*dx,dx)
    
    Pgen_df_to_export = pd.DataFrame()
    Pshm_df_to_export = pd.DataFrame()
    
    for c,cohort in enumerate(cohorts):
        
        print('   * cohort %s' % cohort, end='')
        
        filenames = source_ig_data + "datasets/" + cohort + "_IgG/cohortWide_analysis/"
        filenames += cohort + "_IgG_" + chain_type  + "_uniqueSeqs_" + produc + "_UMIthr_3_trimmed_*_*_noHyperIndels_SHM" + SHMmodel + ".IGoR_summary"
        
        if len(glob.glob(filenames))==1:
            
            df = pd.read_csv(glob.glob(filenames)[0], sep=';', low_memory=False)
            df = df.query('IGoR_Pgen==IGoR_Pgen').reset_index(drop=True)
            df = df.query('IGoR_hyperMutations_likelihood==IGoR_hyperMutations_likelihood').reset_index(drop=True)
            print(' [N=%d]' % (len(df)))
            
            # compute mean and median values for each donor
            
            if cohort in ["healthy_control","hiv1","hcv"]:
                
                samples = sorted(list(set(df['igBlast_sample'])))
                nan_list = [np.nan for i in range(len(samples))]
                sampleAver_df_to_export = pd.DataFrame({'sample': samples, \
                                                        'mean_P_gen': nan_list, 'median_P_gen': nan_list, \
                                                        'mean_P_SHM': nan_list, 'median_P_SHM': nan_list \
                                                       })
                sampleAver_df_to_export.set_index('sample', inplace=True)
                
                for s,sample in enumerate(samples):
                    sample_df = df.query('igBlast_sample==@sample').reset_index(drop=True)
                    
                    sample_data = [np.log10(_) for _ in sample_df['IGoR_Pgen']]
                    sampleAver_df_to_export.at[sample, 'mean_P_gen'] = np.mean(sample_data)
                    sampleAver_df_to_export.at[sample, 'median_P_gen'] = np.median(sample_data)
                    
                    sample_data = [np.log10(_) for _ in sample_df['IGoR_hyperMutations_likelihood']]
                    sampleAver_df_to_export.at[sample, 'mean_P_SHM'] = np.mean(sample_data)
                    sampleAver_df_to_export.at[sample, 'median_P_SHM'] = np.median(sample_data)
                
                sampleAver_df_to_export = sampleAver_df_to_export.reset_index(drop=False)
                sampleAver_df_to_export.to_csv(source_ig_data + 'stats/sampleAveragedValues_Pgen_and_Pshm_' + cohort + '_IgG_' + chain_type + '_' + produc + '_SHM' + SHMmodel + '.csv', index=False, sep=';')
            
            # use bootstrap to estimate mean distribution and error
            
            Pgen_data = [np.log10(_) for _ in df['IGoR_Pgen']]
            Pshm_data = [np.log10(_) for _ in df['IGoR_hyperMutations_likelihood']]
            
            # P_gen
            L = len(Pgen_data)
            a = np.zeros(len(x_Pgen)-1)
            a2 = np.zeros(len(x_Pgen)-1)
            count = 0
            for r in range(T_bootstrap):
                temp_data = [random.choice(Pgen_data) for _ in range(L)]
                hist, bin_edges = np.histogram(temp_data, bins=x_Pgen, density=True)
                a += hist
                a2 += hist*hist
                count += 1
            mean_hist = a/count
            std_hist = np.sqrt((a2/count)-(a/count)**2)
            if(c==0):
                Pgen_df_to_export['P_gen'] = bin_edges[:-1]+0.5*(bin_edges[1]-bin_edges[0])
            Pgen_df_to_export[cohort+'_bootstrap_mean'] = mean_hist
            Pgen_df_to_export[cohort+'_bootstrap_std'] = std_hist
            
            # P_SHM
            L = len(Pshm_data)
            a = np.zeros(len(x_Pshm)-1)
            a2 = np.zeros(len(x_Pshm)-1)
            count = 0
            for r in range(T_bootstrap):
                temp_data = [random.choice(Pshm_data) for _ in range(L)]
                hist, bin_edges = np.histogram(temp_data, bins=x_Pshm, density=True)
                a += hist
                a2 += hist*hist
                count += 1
            mean_hist = a/count
            std_hist = np.sqrt((a2/count)-(a/count)**2)
            if(c==0):
                Pshm_df_to_export['P_SHM'] = bin_edges[:-1]+0.5*(bin_edges[1]-bin_edges[0])
            Pshm_df_to_export[cohort+'_bootstrap_mean'] = mean_hist
            Pshm_df_to_export[cohort+'_bootstrap_std'] = std_hist
            
        else:
            
            print('Ambiguous result for \'' + cohort + '\' cohort...')
    
    Pgen_df_to_export.to_csv(source_ig_data + 'stats/Pgen_distr_' + chain_type + '_' + produc + '_SHM' + SHMmodel + '.csv2', index=False, sep=';')
    Pshm_df_to_export.to_csv(source_ig_data + 'stats/Pshm_distr_' + chain_type + '_' + produc + '_SHM' + SHMmodel + '.csv2', index=False, sep=';')
    
    print()
