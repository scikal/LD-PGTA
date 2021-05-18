#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 10:08:27 2021

@author: ariad
"""
import glob
import bz2, gzip, pickle
filenames = glob.glob("/mybox/simulations/**/simulated*obs.p.bz2", recursive=True)

for filename in filenames:
    try:
        print(filename)
        Open = {'bz2': bz2.open, 'gzip': gzip.open}.get(filename.rpartition('.')[-1], open)
            
        with Open(filename, 'rb') as f:
            obs_tab = pickle.load(f)
            info = pickle.load(f)
        
        if len(obs_tab[0])==4:
            obs_tab2 = tuple((a,c,d) for a,b,c,d in obs_tab)
            
        
        with Open(filename, "wb") as f:
            pickle.dump(obs_tab2, f, protocol=4)
            pickle.dump(info, f, protocol=4)
    except:
        print('failed!')
        
    