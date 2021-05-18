#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  2 01:37:18 2021

@author: ariad
"""
work_dir0 = '/home/ariad/Dropbox/postdoc_JHU/LD-PGTA_ecosystem/build_reference_panel/samples_per_panel/'
REF = {}
for sp in ('EUR','AFR','SAS','EAS','AMR'):
    with open(work_dir0 + f'{sp:s}_panel.txt','rt') as f:
        REF[sp] = [i.strip('\n') for i in f]

def check(*x):
    for sp in ('EUR','AFR','SAS','EAS','AMR','ERROR'):
        #if sp=='ERROR':
        #    print(x)
        if all(i in REF[sp] for i in x): 
            break
    return sp

import glob, shutil
for sp in ('EUR','AFR','SAS','EAS','AMR'):
    work_dir = f'/mybox/simulations2/results_mixed_{sp:s}/'
    filenames = glob.glob(work_dir + '*.p.bz2')
    for source in filenames:
        a = source.rpartition('/')[-1].split('.')
        if a[1]=='monosomy':
            A = a[-4][:-1]
            sp1=check(A)
            if sp1!=sp:
                print(source)
                dest = f'/mybox/simulations2/results_mixed_{sp1:s}/'+source.rpartition('/')[-1] 
                shutil.move(source, dest)                
        elif a[1]=='SPH':
            A = a[-4][:-1]
            B = a[-5][:-1]
            sp1 = check(A,B)
            if sp1!=sp:
                print(source)
                dest = f'/mybox/simulations2/results_mixed_{sp1:s}/'+source.rpartition('/')[-1] 
                shutil.move(source, dest)                
        elif a[1]=='BPH':
            A = a[-4][:-1]
            B = a[-5][:-1]
            C = a[-6][:-1]
            sp1 = check(A,B,C)
            if sp1!=sp:
                print(source)
                dest = f'/mybox/simulations2/results_mixed_{sp1:s}/'+source.rpartition('/')[-1] 
                shutil.move(source, dest)                
        elif a[1]=='transitions':
            A = a[5][:-1]
            B = a[6][:-1]
            C = a[7][:-1]
            sp1 = check(A,B,C)
            if sp1!=sp:
                print(source)
                dest = f'/mybox/simulations2/results_mixed_{sp1:s}/'+source.rpartition('/')[-1] 
                shutil.move(source, dest)                
        
        
        
        
    