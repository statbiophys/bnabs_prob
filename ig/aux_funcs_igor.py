#!/usr/bin/python3

# coding: utf-8

# Filename: aux_funcs_igor.py
# Author: Cosimo Lupo, <cosimo.lupo89@gmail.com>
# Last updated: 03 July 2023

import numpy as np
import os
#import sys
#import pandas as pd
#import Bio
#from Bio import SeqIO
#import re
import yaml
import itertools
import scipy

########## settings ##########

with open('config_igor.yaml') as config_file:
    config = yaml.full_load(config_file)

#blast_database = config['blast_database']

# Apparently, os.system() is deprecated now
# and subprocess modules have to used.
# More info at:
# https://docs.python.org/3/library/subprocess.html

########## defs ##########

def read_marginals_from_file(in_file_marginals, in_file_parms):
  
  has_Dgene = False
  
  marginals = {}
  
  with open(in_file_marginals, "r") as f:
    rows = f.readlines()
    rows = [_.strip() for _ in rows]
    
    for r,row in enumerate(rows):
      
      if row=='@v_choice':
        N_V = int(rows[r+1][5:-1])
        P_V = [float(_) for _ in rows[r+3][1:].split(',')]
        P_V = np.array(P_V)
        P_V /= np.sum(P_V)
      
      if row=='@j_choice':
        N_J = int(rows[r+1][5:-1].split(',')[1])
        P_JgivenV = np.zeros((N_V,N_J))
        for v in range(N_V):
          x = np.array([float(_) for _ in rows[r+3+2*v][1:].split(',')])
          for j in range(N_J):
            P_JgivenV[v,j] = x[j]
      
      if row=='@d_gene':
        has_Dgene = True
        N_D = int(rows[r+1][5:-1].split(',')[2])
        P_DgivenVJ = np.zeros((N_V,N_D,N_J))
        for v in range(N_V):
          for j in range(N_J):
            x = np.array([float(_) for _ in rows[r+3+2*(v*N_J+j)][1:].split(',')])
            for d in range(N_D):
              P_DgivenVJ[v,d,j] = x[d]
      
      if row=='@v_3_del':
        L_V3del = int(rows[r+1][5:-1].split(',')[1])
        P_V3delgivenV = np.zeros((N_V,L_V3del))
        for v in range(N_V):
          x = np.array([float(_) for _ in rows[r+3+2*v][1:].split(',')])
          for n in range(L_V3del):
            P_V3delgivenV[v,n] = x[n]
      
      if row=='@j_5_del':
        L_J5del = int(rows[r+1][5:-1].split(',')[1])
        P_J5delgivenJ = np.zeros((N_J,L_J5del))
        for j in range(N_J):
          x = np.array([float(_) for _ in rows[r+3+2*j][1:].split(',')])
          for n in range(L_J5del):
            P_J5delgivenJ[j,n] = x[n]
      
      if row=='@d_5_del':
        has_Dgene = True
        L_D5del = int(rows[r+1][5:-1].split(',')[1])
        P_D5delgivenD = np.zeros((N_D,L_D5del))
        for d in range(N_D):
          x = np.array([float(_) for _ in rows[r+3+2*d][1:].split(',')])
          for n in range(L_D5del):
            P_D5delgivenD[d,n] = x[n]
      
      if row=='@d_3_del':
        has_Dgene = True
        L_D3del = int(rows[r+1][5:-1].split(',')[2])
        P_D3delgivenDD5del = np.zeros((N_D,L_D5del,L_D3del))
        for d in range(N_D):
          for n in range(L_D5del):
            x = np.array([float(_) for _ in rows[r+3+2*(d*L_D5del+n)][1:].split(',')])
            for m in range(L_D3del):
              P_D3delgivenDD5del[d,n,m] = x[m]
      
      if row=='@vd_ins':
        has_Dgene = True
        N_VDins = int(rows[r+1][5:-1])
        P_VDins = [float(_) for _ in rows[r+3][1:].split(',')]
        P_VDins = np.array(P_VDins)
        P_VDins /= np.sum(P_VDins)
      
      if row=='@dj_ins':
        has_Dgene = True
        N_DJins = int(rows[r+1][5:-1])
        P_DJins = [float(_) for _ in rows[r+3][1:].split(',')]
        P_DJins = np.array(P_DJins)
        P_DJins /= np.sum(P_DJins)
      
      if row=='@vj_ins':
        has_Dgene = False
        N_VJins = int(rows[r+1][5:-1])
        P_VJins = [float(_) for _ in rows[r+3][1:].split(',')]
        P_VJins = np.array(P_VJins)
        P_VJins /= np.sum(P_VJins)
    
    P_VJ = np.zeros((N_V,N_J))
    for v in range(N_V):
      for j in range(N_J):
        P_VJ[v,j] = P_JgivenV[v,j]*P_V[v]
    P_VJ /= np.sum(P_VJ)
    
    P_J = np.zeros(N_J)
    for j in range(N_J):
      P_J[j] = np.sum(P_VJ[:,j])
    
    P_VV3del = np.zeros((N_V,L_V3del))
    for v in range(N_V):
      for n in range(L_V3del):
        P_VV3del[v,n] = P_V3delgivenV[v,n]*P_V[v]
    P_VV3del /= np.sum(P_VV3del)
    
    P_JJ5del = np.zeros((N_J,L_J5del))
    for j in range(N_J):
      for n in range(L_J5del):
        P_JJ5del[j,n] = P_J5delgivenJ[j,n]*P_J[j]
    P_JJ5del /= np.sum(P_JJ5del)
    
    if has_Dgene:
      
      P_VDJ = np.zeros((N_V,N_D,N_J))
      for v in range(N_V):
        for d in range(N_D):
          for j in range(N_J):
            P_VDJ[v,d,j] = P_DgivenVJ[v,d,j]*P_VJ[v,j]
      P_VDJ /= np.sum(P_VDJ)
      
      P_D = np.zeros(N_D)
      for d in range(N_D):
        P_D[d] = np.sum(P_VDJ[:,d,:])
      
      P_DD5del = np.zeros((N_D,L_D5del))
      for d in range(N_D):
        for n in range(L_D5del):
          P_DD5del[d,n] = P_D5delgivenD[d,n]*P_D[d]
      P_DD5del /= np.sum(P_DD5del)
      
      P_DD5delD3del = np.zeros((N_D,L_D5del,L_D3del))
      for d in range(N_D):
        for n in range(L_D5del):
          for m in range(L_D3del):
            P_DD5delD3del[d,n,m] = P_D3delgivenDD5del[d,n,m]*P_DD5del[d,n]
      P_DD5delD3del /= np.sum(P_DD5delD3del)
    
    marginals['V'] = P_V
    marginals['J'] = P_J
    marginals['VJ'] = P_VJ
    marginals['JgivenV'] = P_JgivenV
    marginals['VV3del'] = P_VV3del
    marginals['V3delgivenV'] = P_V3delgivenV
    marginals['JJ5del'] = P_JJ5del
    marginals['J5delgivenJ'] = P_J5delgivenJ
    if has_Dgene:
      marginals['D'] = P_D
      marginals['VDJ'] = P_VDJ
      marginals['DgivenVJ'] = P_DgivenVJ
      marginals['DD5del'] = P_DD5del
      marginals['DD5delD3del'] = P_DD5delD3del
      marginals['D5delgivenD'] = P_D5delgivenD
      marginals['D3delgivenDD5del'] = P_D3delgivenDD5del
      marginals['VDins'] = P_VDins
      marginals['DJins'] = P_DJins
    else:
      marginals['VJins'] = P_VJins
    
  with open(in_file_parms, "r") as f:
    rows = f.readlines()
    rows = [_.strip() for _ in rows]
    
    weighted_aver_Vlen = 0
    
    for r,row in enumerate(rows):
      
      if row=='@Event_list':
        for v in range(N_V):
          L = len(rows[r+2+v][1:].split(';')[1])
          i = int(rows[r+2+v][1:].split(';')[2])
          weighted_aver_Vlen += L*P_V[i]
          
      if row=='@ErrorRate':
        if 'HypermutationfullNmererrorrate' in rows[r+1]:
          Nmer_size = int(rows[r+1][1:].split(';')[1])
          P_Nmer = [float(_) for _ in rows[r+2].split(';')]
          marginals['Nmer'] = P_Nmer
          marginals['weighted_aver_Vlen'] = weighted_aver_Vlen
  
  return marginals

def regularize_joint_marginals(marginals, N_seqs):
  
  regularized_marginals = {}
  
  if 'VDJ' in marginals.keys():
    # full VDJ
    has_Dgene = True
    N_V, N_D, N_J = np.shape(marginals['VDJ'])
    
    regularized_marginals['VDJ'] = np.zeros((N_V,N_D,N_J))
    for v in range(N_V):
      for d in range(N_D):
        for j in range(N_J):
          regularized_marginals['VDJ'][v,d,j] = (marginals['VDJ'][v,d,j]*N_seqs + 1) / (N_seqs + N_V*N_D*N_J)
    
    regularized_marginals['V'] = np.zeros(N_V)
    for v in range(N_V):
      regularized_marginals['V'][v] = np.sum(regularized_marginals['VDJ'][v,:,:])
    
    regularized_marginals['D'] = np.zeros(N_D)
    for d in range(N_D):
      regularized_marginals['D'][d] = np.sum(regularized_marginals['VDJ'][:,d,:])
    
    regularized_marginals['J'] = np.zeros(N_J)
    for j in range(N_J):
      regularized_marginals['J'][j] = np.sum(regularized_marginals['VDJ'][:,:,j])
    
    regularized_marginals['VJ'] = np.zeros((N_V,N_J))
    for v in range(N_V):
      for j in range(N_J):
        regularized_marginals['VJ'][v,j] = np.sum(regularized_marginals['VDJ'][v,:,j])
    
    N_D, L_D5del, L_D3del = np.shape(marginals['DD5delD3del'])
    regularized_marginals['DD5delD3del'] = np.zeros((N_D,L_D5del,L_D3del))
    for d in range(N_D):
      for n in range(L_D5del):
        for m in range(L_D3del):
          regularized_marginals['DD5delD3del'][d,n,m] = (marginals['DD5delD3del'][d,n,m]*N_seqs + 1) / (N_seqs + N_D*L_D5del*L_D3del)
    
    regularized_marginals['DD5del'] = np.zeros((N_D,L_D5del))
    for d in range(N_D):
      for n in range(L_D5del):
        regularized_marginals['DD5del'][d,n] = np.sum(regularized_marginals['DD5delD3del'][d,n,:])
    
    regularized_marginals['D3delgivenDD5del'] = np.zeros((N_D,L_D5del,L_D3del))
    for d in range(N_D):
      for n in range(L_D5del):
        for m in range(L_D3del):
          if regularized_marginals['DD5delD3del'][d,n,m]>0:
            regularized_marginals['D3delgivenDD5del'][d,n,m] = regularized_marginals['DD5delD3del'][d,n,m]/regularized_marginals['DD5del'][d,n]
    
    regularized_marginals['D5delgivenD'] = np.zeros((N_D,L_D5del))
    for d in range(N_D):
      for n in range(L_D5del):
        if regularized_marginals['DD5del'][d,n]>0:
          regularized_marginals['D5delgivenD'][d,n] = regularized_marginals['DD5del'][d,n]/regularized_marginals['D'][d]
    
    N_VDins = len(marginals['VDins'])
    regularized_marginals['VDins'] = np.zeros(N_VDins)
    for n in range(N_VDins):
      regularized_marginals['VDins'][n] = (marginals['VDins'][n]*N_seqs + 1) / (N_seqs + N_VDins)
    
    N_DJins = len(marginals['DJins'])
    regularized_marginals['DJins'] = np.zeros(N_DJins)
    for n in range(N_DJins):
      regularized_marginals['DJins'][n] = (marginals['DJins'][n]*N_seqs + 1) / (N_seqs + N_DJins)
    
    regularized_marginals['DgivenVJ'] = np.zeros((N_V,N_D,N_J))
    for v in range(N_V):
      for d in range(N_D):
        for j in range(N_J):
          if regularized_marginals['VDJ'][v,d,j]>0:
            regularized_marginals['DgivenVJ'][v,d,j] = regularized_marginals['VDJ'][v,d,j]/regularized_marginals['VJ'][v,j]
    
  else:
    # only V and J
    has_Dgene = False
    N_V, N_J = np.shape(marginals['VJ'])
    
    regularized_marginals['VJ'] = np.zeros((N_V,N_J))
    for v in range(N_V):
      for j in range(N_J):
        regularized_marginals['VJ'][v,j] = (marginals['VJ'][v,j]*N_seqs + 1) / (N_seqs + N_V*N_J)
    
    regularized_marginals['V'] = np.zeros(N_V)
    for v in range(N_V):
      regularized_marginals['V'][v] = np.sum(regularized_marginals['VJ'][v,:])
    
    regularized_marginals['J'] = np.zeros(N_J)
    for j in range(N_J):
      regularized_marginals['J'][j] = np.sum(regularized_marginals['VJ'][:,j])
    
    N_VJins = len(marginals['VJins'])
    regularized_marginals['VJins'] = np.zeros(N_VJins)
    for n in range(N_VJins):
      regularized_marginals['VJins'][n] = (marginals['VJins'][n]*N_seqs + 1) / (N_seqs + N_VJins)
  
  regularized_marginals['JgivenV'] = np.zeros((N_V,N_J))
  for v in range(N_V):
    for j in range(N_J):
      if regularized_marginals['VJ'][v,j]>0:
        regularized_marginals['JgivenV'][v,j] = regularized_marginals['VJ'][v,j]/regularized_marginals['V'][v]
  
  N_V, L_V3del = np.shape(marginals['VV3del'])
  regularized_marginals['VV3del'] = np.zeros((N_V,L_V3del))
  for v in range(N_V):
      for n in range(L_V3del):
        regularized_marginals['VV3del'][v,n] = (marginals['VV3del'][v,n]*N_seqs + 1) / (N_seqs + N_V*L_V3del)
  
  regularized_marginals['V3delgivenV'] = np.zeros((N_V,L_V3del))
  for v in range(N_V):
    for n in range(L_V3del):
      if regularized_marginals['VV3del'][v,n]>0:
        regularized_marginals['V3delgivenV'][v,n] = regularized_marginals['VV3del'][v,n]/regularized_marginals['V'][v]
  
  N_J, L_J5del = np.shape(marginals['JJ5del'])
  regularized_marginals['JJ5del'] = np.zeros((N_J,L_J5del))
  for j in range(N_J):
      for n in range(L_J5del):
        regularized_marginals['JJ5del'][j,n] = (marginals['JJ5del'][j,n]*N_seqs + 1) / (N_seqs + N_J*L_J5del)
  
  regularized_marginals['J5delgivenJ'] = np.zeros((N_J,L_J5del))
  for j in range(N_J):
    for n in range(L_J5del):
      if regularized_marginals['JJ5del'][j,n]>0:
        regularized_marginals['J5delgivenJ'][j,n] = regularized_marginals['JJ5del'][j,n]/regularized_marginals['J'][j]
  
  if 'Nmer' in marginals.keys():
    # regularization v1
    """
    S = np.sum(marginals['Nmer'])
    N = len(marginals['Nmer'])
    L = marginals['weighted_aver_Vlen']
    regularized_marginals['Nmer'] = np.zeros(N)
    for n in range(N):
      regularized_marginals['Nmer'][n] = (marginals['Nmer'][n]*N_seqs*L + 1) / (N_seqs*L + N/S)
    """
    # regularization v2
    S = np.sum(marginals['Nmer'])
    N = len(marginals['Nmer'])
    alphabet = ['A','C','G','T']
    Nmer_size = int(np.log10(N)/np.log10(len(alphabet)))
    allNmers = [''.join(i) for i in itertools.product(alphabet, repeat=Nmer_size)]
    eps = 1e-3
    regularized_marginals['Nmer'] = np.zeros(N)
    for n in range(N):
      neighbours = [i for (i,_) in enumerate(allNmers) if scipy.spatial.distance.hamming(list(allNmers[n]),list(_))==1/Nmer_size]
      regularized_marginals['Nmer'][n] = (1-eps)*marginals['Nmer'][n] + eps*np.mean([marginals['Nmer'][i] for i in neighbours])
    
  return regularized_marginals

def get_all_marginals(joint_marginal):
  
  if len(np.shape(joint_marginal))==3:
    # full VDJ
    has_Dgene = True
    N_V, N_D, N_J = np.shape(joint_marginal)
    P_VDJ = joint_marginal
    P_V = np.zeros(N_V)
    for v in range(N_V):
      P_V[v] = np.sum(joint_marginal[v,:,:])
    P_D = np.zeros(N_D)
    for d in range(N_D):
      P_D[d] = np.sum(joint_marginal[:,d,:])
    P_J = np.zeros(N_J)
    for j in range(N_J):
      P_J[j] = np.sum(joint_marginal[:,:,j])
    P_VJ = np.zeros((N_V,N_J))
    for v in range(N_V):
      for j in range(N_J):
        P_VJ[v,j] = np.sum(joint_marginal[v,:,j])
    P_JgivenV = np.zeros((N_V,N_J))
    for v in range(N_V):
      for j in range(N_J):
        if P_VJ[v,j]>0:
          P_JgivenV[v,j] = P_VJ[v,j]/P_V[v]
    P_DgivenVJ = np.zeros((N_V,N_D,N_J))
    for v in range(N_V):
      for d in range(N_D):
        for j in range(N_J):
          if joint_marginal[v,d,j]>0:
            P_DgivenVJ[v,d,j] = joint_marginal[v,d,j]/P_VJ[v,j]
  
  elif len(np.shape(joint_marginal))==2:
    # only V and J
    has_Dgene = False
    N_V, N_J = np.shape(joint_marginal)
    P_VJ = joint_marginal
    P_V = np.zeros(N_V)
    for v in range(N_V):
      P_V[v] = np.sum(joint_marginal[v,:])
    P_J = np.zeros(N_J)
    for j in range(N_J):
      P_J[j] = np.sum(joint_marginal[:,j])
    P_JgivenV = np.zeros((N_V,N_J))
    for v in range(N_V):
      for j in range(N_J):
        if joint_marginal[v,j]>0:
          P_JgivenV[v,j] = joint_marginal[v,j]/P_V[v]
  
  # sanity checks
  if np.abs(np.sum(P_V)-1)>1e-6:
    print('Error in the normalization of P_V!')
  if has_Dgene:
    if np.abs(np.sum(P_D)-1)>1e-6:
      print('Error in the normalization of P_D!')
  if np.abs(np.sum(P_J)-1)>1e-6:
    print('Error in the normalization of P_J!')
  if np.abs(np.sum(P_VJ)-1)>1e-6:
    print('Error in the normalization of P_VJ!')
  if has_Dgene:
    if np.abs(np.sum(P_VDJ)-1)>1e-6:
      print('Error in the normalization of P_VDJ!')
  if np.abs(np.sum(P_JgivenV[0,:])-1)>1e-6:
    print('Error in the normalization of P_JgivenV!')
  if has_Dgene:
    if np.abs(np.sum(P_DgivenVJ[0,:,0])-1)>1e-6:
      print('Error in the normalization of P_DgivenVJ!')
  
  if has_Dgene:
    return P_V, P_D, P_J, P_VJ, P_JgivenV, P_DgivenVJ
  else:
    return P_V, P_J, P_JgivenV

def write_marginals_on_file(in_file_marginals, in_file_parms, N_seqs):
  
  original_marginals = read_marginals_from_file(in_file_marginals, in_file_parms)
  regularized_marginals = regularize_joint_marginals(original_marginals, N_seqs)
  
  #original_joint_marginal = read_marginals(in_file_marginals)
  """
  if len(np.shape(original_joint_marginal))==3:
    # full VDJ
    has_Dgene = True
    N_V, N_D, N_J = np.shape(original_joint_marginal)
    P_VDJ = original_joint_marginal
    P_VDJ_reg = regularize_joint_marginals(original_joint_marginal, N_seqs)
    P_V_reg, P_D_reg, P_J_reg, P_VJ_reg, P_JgivenV_reg, P_DgivenVJ_reg = get_all_marginals(P_VDJ_reg)
    
  elif len(np.shape(original_joint_marginal))==2:
    # only V and J
    has_Dgene = False
    N_V, N_J = np.shape(original_joint_marginal)
    P_VJ = original_joint_marginal
    P_VJ_reg = regularize_joint_marginals(original_joint_marginal, N_seqs)
    P_V_reg, P_J_reg, P_JgivenV_reg = get_all_marginals(P_VJ_reg)
  """
  
  # Updating marginals file
  
  with open(in_file_marginals, "r") as f:
    original_rows = f.readlines()
    original_rows = [_.strip() for _ in original_rows]
  
  regularized_rows = original_rows
  
  for r,row in enumerate(original_rows):
    
    if row=='@v_choice':
      N_V = len(regularized_marginals['V'])
      regularized_rows[r+3] = '%' + ','.join(['%.6g' % _ for _ in regularized_marginals['V']])
    
    if row=='@j_choice':
      N_J = len(regularized_marginals['J'])
      for v in range(N_V):
        regularized_rows[r+3+2*v] = '%' + ','.join(['%.6g' % _ for _ in regularized_marginals['JgivenV'][v,:]])
    
    if row=='@d_gene':
      N_D = len(regularized_marginals['D'])
      for v in range(N_V):
        for j in range(N_J):
          regularized_rows[r+3+2*(v*N_J+j)] = '%' + ','.join(['%.6g' % _ for _ in regularized_marginals['DgivenVJ'][v,:,j]])
    
    if row=='@v_3_del':
      L_V3del = np.shape(regularized_marginals['VV3del'])[1]
      for v in range(N_V):
        regularized_rows[r+3+2*v] = '%' + ','.join(['%.6g' % _ for _ in regularized_marginals['V3delgivenV'][v,:]])
    
    if row=='@j_5_del':
      L_J5del = np.shape(regularized_marginals['JJ5del'])[1]
      for j in range(N_J):
        regularized_rows[r+3+2*j] = '%' + ','.join(['%.6g' % _ for _ in regularized_marginals['J5delgivenJ'][j,:]])
    
    if row=='@d_5_del':
      L_D5del = np.shape(regularized_marginals['DD5del'])[1]
      for d in range(N_D):
        regularized_rows[r+3+2*d] = '%' + ','.join(['%.6g' % _ for _ in regularized_marginals['D5delgivenD'][d,:]])
    
    if row=='@d_3_del':
      L_D3del = np.shape(regularized_marginals['DD5delD3del'])[2]
      for d in range(N_D):
        for n in range(L_D5del):
          regularized_rows[r+3+2*(d*L_D5del+n)] = '%' + ','.join(['%.6g' % _ for _ in regularized_marginals['D3delgivenDD5del'][d,n,:]])
    
    if row=='@vd_ins':
      regularized_rows[r+3] = '%' + ','.join(['%.6g' % _ for _ in regularized_marginals['VDins']])
    
    if row=='@dj_ins':
      regularized_rows[r+3] = '%' + ','.join(['%.6g' % _ for _ in regularized_marginals['DJins']])
    
    if row=='@vj_ins':
      regularized_rows[r+3] = '%' + ','.join(['%.6g' % _ for _ in regularized_marginals['VJins']])
    
    if regularized_rows[r][-1]!='\n':
      regularized_rows[r] += '\n'
  
  out_file_marginals = os.path.dirname(in_file_marginals) + '/' + os.path.splitext(os.path.basename(in_file_marginals))[0] + '_regularized' + os.path.splitext(os.path.basename(in_file_marginals))[1]
  
  with open(out_file_marginals, "w") as f:
    f.writelines(regularized_rows)
  
  # Updating parms file
  
  with open(in_file_parms, "r") as f:
    original_rows = f.readlines()
    original_rows = [_.strip() for _ in original_rows]
  
  regularized_rows = original_rows
  
  for r,row in enumerate(original_rows):
    
    if row=='@ErrorRate':
      if ('HypermutationfullNmererrorrate' in original_rows[r+1]) and ('Nmer' in original_marginals.keys()):
        regularized_rows[r+2] = ';'.join(['%.6g' % _ for _ in regularized_marginals['Nmer']])
  
    if regularized_rows[r][-1]!='\n':
      regularized_rows[r] += '\n'
  
  out_file_parms = os.path.dirname(in_file_parms) + '/' + os.path.splitext(os.path.basename(in_file_parms))[0] + '_regularized' + os.path.splitext(os.path.basename(in_file_parms))[1]
  
  with open(out_file_parms, "w") as f:
    f.writelines(regularized_rows)
  
  return
