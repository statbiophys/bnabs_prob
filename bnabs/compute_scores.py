#!/usr/bin/python3

# coding: utf-8

# Filename: compute_scores.py
# Author: Cosimo Lupo, <cosimo.lupo89@gmail.com>
# Last updated: 13 July 2023

import pandas as pd
import numpy as np
import glob as glob
import os as os
import yaml
from sklearn import linear_model
from scipy.stats import linregress
clf = linear_model.LinearRegression()

# Bulding scores

chainType = "LC"
produc = "cohortWide_P"
SHMmodel = "full5mer"
cohorts = ["healthy_control","hiv1","hcv",\
           "hiv1_non_neutralizers","hiv1_weak_neutralizers",\
           "hiv1_intermediate_neutralizers","hiv1_top_neutralizers",\
           "hiv1_art_off","hiv1_art_on"]
log_AuC = True
useRegularizedMarginals = True

source_bnabs_data = "./"
prefix_name = "bnabs_seqs_igor_cohortWide_"

print()

for fit_model in [1,2,7]:
    
    print('Fit model n. %d' % fit_model)
    
    for c,cohort in enumerate(cohorts):
        
        print(' > cohort: ' + cohort, end = '')
        
        in_file = prefix_name + cohort + "_" + chainType + "_"
        if produc=="cohortWide_NP":
            in_file += "NP"
        elif produc=="cohortWide_P":
            in_file += "P"
        in_file += "_SHM" + SHMmodel
        if useRegularizedMarginals:
            in_file += "_regularized"
        in_file += ".df"
        
        df = pd.read_csv(source_bnabs_data + "igor_bnabs_summary/" + in_file, sep=';')
        N_before = len(df)
        df = df.query('IGoR_noHyperIndels_seq_likelihood==IGoR_noHyperIndels_seq_likelihood')
        df = df.query('IGoR_noHyperIndels_seq_likelihood>0')
        df = df.query('igBlast_CDR3_nt==igBlast_CDR3_nt')
        df = df.query('AuC>10')
        N_after = len(df)
        
        print(' [%d/%d bnabs]' % (N_after,N_before))
        
        if c==0:
            if log_AuC:
                df['log10_AuC'] = df['AuC'].apply(lambda x: np.log10(x))
                main_df = df[['bnab_ID','log10_AuC']]
            else:
                main_df = df[['bnab_ID','AuC']]
        
        if fit_model==1:
            # Model n. 1: only P_gen
            X = [np.log10(df['IGoR_noHyperIndels_Pgen'])]
            fit_model_eq = '$\mathcal{S} \equiv c_1\,log_{10}{\,P_{gen}}$'
        elif fit_model==2:
            # Model n. 2: only P_SHM
            X = [np.log10(df['IGoR_noHyperIndels_hyperMutations_likelihood'])]
            fit_model_eq = '$\mathcal{S} \equiv c_1\,log_{10}{\,P_{SHM}}$'
        elif fit_model==7:
            # Model n. 7: P_gen & P_SHM
            X = [np.log10(df['IGoR_noHyperIndels_Pgen']),\
                 np.log10(df['IGoR_noHyperIndels_hyperMutations_likelihood'])]
            fit_model_eq = '$\mathcal{S} \equiv c_1\,log_{10}{\,P_{gen}}+c_2\,log_{10}{\,P_{SHM}}$'
        elif fit_model==21:
            X = [np.log10(df['IGoR_noHyperIndels_Pgen']),\
                 np.log10(df['IGoR_noHyperIndels_hyperMutations_likelihood']),\
                 df['igBlast_N_indels']]
            fit_model_eq = '$\mathcal{S} \equiv c_1\,log_{10}{\,P_{gen}}+c_2\,log_{10}{\,P_{SHM}}+c_3\,N_I$'
        
        X = np.reshape(np.array(X).T,(len(X[0]),len(X)))
        if log_AuC:
            Y = np.log10(df['AuC'])
        else:
            Y = df['AuC']
        clf.fit(X,Y)
        Z = np.dot(X,clf.coef_)
        
        features = [",".join([str(x) for x in X[i]]) for i in range(len(X))]
        coeffs = ",".join([str(x) for x in clf.coef_])
        
        slope, intercept, r_value, p_value, std_err = linregress(Z,Y)
        print('     r = %.3g ,  p = %.3g' % (r_value,p_value))
        
        df_out = pd.DataFrame({'bnab_ID': np.array(df['bnab_ID']), \
                               'features': features, \
                               'fit_coeffs': [coeffs for i in range(len(X))], \
                               'fit_score': Z, \
                               'fit_slope': [slope for i in range(len(X))], \
                               'fit_intercept': [intercept for i in range(len(X))], \
                               'fit_r_value': [r_value for i in range(len(X))], \
                               'fit_p_value': [p_value for i in range(len(X))] \
                              })
        
        for col in df_out.columns.values.tolist():
            if(col!='bnab_ID'):
                df_out.rename(columns={col: col+'_'+cohort}, inplace=True)
        
        main_df = main_df.merge(df_out, how='left', on='bnab_ID')
        main_df.fillna(np.nan, inplace=True)
    
    if not os.path.isdir(source_bnabs_data + "neutr_score_prediction"):
        os.mkdir(source_bnabs_data + "neutr_score_prediction")
    
    score_file_name = source_bnabs_data + "neutr_score_prediction/"
    if log_AuC:
        score_file_name += "logAuC"
    else:
        score_file_name += "AuC"
    score_file_name += "_score_" + chainType + "_" + produc + "_SHM" + SHMmodel + "_fitModel_" + ("%02d" % fit_model)
    if useRegularizedMarginals:
        score_file_name += "_regularized"
    score_file_name += ".csv"
    main_df.to_csv(score_file_name, index=False, sep=';')
    
    print()
