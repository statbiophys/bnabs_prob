#!/usr/bin/python3

# coding: utf-8

# Filename: aux_funcs_igor.py
# Author: Cosimo Lupo, <cosimo.lupo89@gmail.com>
# Last updated: 10 February 2023

def expand_array(strvar):
  if(len(strvar)>2):
    return [int(xstr) for xstr in strvar[1:-1].split(',')]
  else:
    return []

def expand_array_2(strvar):
  if(len(strvar)>2):
    return [int(xstr) for xstr in strvar[1:-1].split('|')]
  else:
    return []

def expand_array_3(strvar):
  if(len(strvar)>2):
    return [xstr for xstr in strvar[1:-1].split('|')]
  else:
    return []

def get_CDR3_len(CDR3):
  if(CDR3!="N/A"):
    return len(CDR3)
  else:
    return 0

def get_family(gene):
  return (gene.split(',')[0]).split('*')[0]

def parse_hyper_indels(errors_list,indel_type):
  if(indel_type=="del"):
    if(len(errors_list)>0):
      v = [y for y in errors_list if y[0]=="d"]
      if(len(v)>0):
        return [len(y.split('|')[3]) for y in v]
      else:
        return []
    else:
      return []
  elif(indel_type=="ins"):
    if(len(errors_list)>0):
      v = [y for y in errors_list if y[0]=="i"]
      if(len(v)>0):
        return [len(y.split('|')[3]) for y in v]
      else:
        return []
    else:
      return []
