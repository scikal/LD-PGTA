#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""" Creates statistical models for the case that two haplotypes are associated 
    with one reference panel and the third haplotype with a different reference
    panel. """

from itertools import product
from collections import defaultdict
from math import gcd as greatest_common_divisor
from pickle import dump
from time import time
from bz2 import BZ2File
import argparse, sys

def ENGINE(number_of_reads,degeneracies):
    """ Generates polysomy statistical models for n-reads, based of a list of
    the degeneracy of each homolog. """

    degeneracies_dict = {i:w for i,w in enumerate(degeneracies) if w>0}
    model = defaultdict(int)
    for sequence in product(degeneracies_dict, repeat=number_of_reads):
        haplotypes = defaultdict(list,{0:[]})
        weight = 1
        for read_ind,hap in enumerate(sequence):
            haplotypes[hap].append(read_ind)
            weight *=  degeneracies_dict[hap]
        key = tuple(tuple(indices) for indices in haplotypes.values())
        model[key] += weight
    return model

def representationA(model):
     """ Represents the model that is returned from the function ENGINE. """
     result = ''
     for partition,weight in model.items():
        if result!='':
            result += '+'
        l = [''.join((chr(read+65) for read in hap)) for hap in partition  if len(hap)] 
        sandwitch = ''.join(prefix+hap+')' for prefix,hap in zip(('g(','f(','f('),l))
        if weight!=1:
            result += f'{weight:d}'
        result += f'{sandwitch:s}'
     return result

def representationB(model):
    """ An alternative representation of the model. """
    result = ''
    for partition,weight in model.items():
        if result!='':
            result += '+'
        l = [''.join((chr(Aa+read_ind) for read_ind in hap)) for Aa,hap in zip((65,97,97),partition) if len(hap)] 
        if weight!=1:
            result += f'{weight:d}*'
        result +=  f"{'*'.join(l):s}"
    return result


def BPH(number_of_reads):
    """ Builds a statistical model for n-reads under the BPH scenario. """

    degeneracies = (1, 1, 1)
    model = ENGINE(number_of_reads,degeneracies)
    return model

def SPH(number_of_reads):
    """ Builds a statistical model for n-reads under the SPH scenario. """

    degeneracies = (1, 2)
    model = ENGINE(number_of_reads,degeneracies)
    return model

def DISOMY(number_of_reads):
    """ Builds a statistical model for n-reads under the diploidy scenario. """

    degeneracies = (1, 1)
    model = ENGINE(number_of_reads,degeneracies)
    return model

def MONOSOMY(number_of_reads):
    """ Builds a statistical model for n-reads under the monosomy scenario. """

    ### return {1: {(1,1): int(number_of_reads*'1',2)}}
    degeneracies = (1,)
    model = ENGINE(number_of_reads,degeneracies)
    return model


if __name__ == "__main__":
    print("--- Two reads ---")
    print("Disomy: ",representationB(DISOMY(2)))
    print("SPH:",representationB(SPH(2)))
    print("BPH:",representationB(BPH(2)))
    print("--- Three reads ---")
    print("Disomy:",representationB(DISOMY(3)))
    print("SPH:",representationB(SPH(3)))
    print("BPH:",representationB(BPH(3)))
    print("--- Four reads ---")
    print("Disomy:",representationB(DISOMY(4)))
    print("SPH:",representationB(SPH(4)))
    print("BPH:",representationB(BPH(4)))
    sys.exit(0)

else:
    print("The module 'verify' was imported.")
