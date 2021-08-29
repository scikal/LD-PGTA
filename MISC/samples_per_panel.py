#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 22:46:22 2021

@author: ariad
"""
from operator import itemgetter

with open('ALL_panel.txt', 'r') as f:
    SAMPLES = {line.strip() for line in f}

with open('igsr_samples.tsv', 'r') as f:
    f.readline()   # Bite off the header
    DATA = (line.strip().split('\t') for line in f)
    RELEVANT = sorted((i for i in DATA if i[0] in SAMPLES and '30x' in i[-1]), key=itemgetter(5,3,1,0))
    
    
header = "sample population group sex\n"
C = {}
with open("EUR_panel.samples", "w") as C['EUR'], open("EAS_panel.samples", "w") as C['EAS'], open("SAS_panel.samples", "w") as C['SAS'], open("AMR_panel.samples", "w") as C['AMR'], open("AFR_panel.samples", "w") as C['AFR']:
    for superpopulation_code in ('EUR','EAS','SAS','AMR','AFR'):
        C[superpopulation_code].write(header)
    for sample_name,sex,biosample_ID,population_code,population_name,superpopulation_code,superpopulation_name,population_elastic_ID,data_collections in RELEVANT:
        C[superpopulation_code].write(f"{sample_name:s} {population_code:s} {superpopulation_code:s} {1 if sex=='male' else 2:d}\n")
             
with open("EUR_panel.txt", "w") as C['EUR'], open("EAS_panel.txt", "w") as C['EAS'], open("SAS_panel.txt", "w") as C['SAS'], open("AMR_panel.txt", "w") as C['AMR'], open("AFR_panel.txt", "w") as C['AFR']:
    for sample_name,sex,biosample_ID,population_code,population_name,superpopulation_code,superpopulation_name,population_elastic_ID,data_collections in RELEVANT:
        C[superpopulation_code].write(f"{sample_name:s}\n")
