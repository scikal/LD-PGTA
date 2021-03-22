#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LIKELIHOODS_CALCULATOR_HOMOGENOUES

Given reads that originated form the same genomic window and a reference panel
of a single population, the likelihood of observed reads under four scenarios,
namely, monosomy, disomy, SPH and BPH is calculated.

BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.

Daniel Ariad (daniel@ariad.org)
Dec 21, 2020
"""

import pickle, os, sys, bz2

from functools import reduce
from operator import and_, itemgetter
from itertools import combinations

try:
    from gmpy2 import popcount
except ModuleNotFoundError:
    print('caution: the module gmpy2 is missing.')
    def popcount(x):
        """ Counts non-zero bits in positive integer. """
        return bin(x).count('1')

class examine_homogeneous:
    """ Based on two IMPUTE2 arrays, which contain the legend and haplotypes,
    and a dictionary with statisitcal models (models_dict), it allows to
    calculate the likelihoods of observed alleles under various statistical
    models (monosomy, disomy, SPH and BPH). """
   
    
    def __init__(self, obs_tab, leg_tab, hap_tab, sam_tab, models_dict, number_of_haplotypes):
        """ Initialize the attributes of the class. """
        self.models_dict = models_dict
        self.hap_dict = self.build_hap_dict(obs_tab, leg_tab, hap_tab, number_of_haplotypes)
        self.total_number_of_haplotypes_in_reference_panel = number_of_haplotypes

    def build_hap_dict(self, obs_tab, leg_tab, hap_tab, number_of_haplotypes):
        """ Returns a dictionary that lists chromosome positions of SNPs and gives
            their relevent row from haplotypes table. The row is stored as a tuple
            of booleans, where True represents the observed allele. We denote
            the returned dictionary as the reference panel."""
    
        hap_dict = dict()
        mismatches = 0
        
        b = (1 << number_of_haplotypes) - 1 #### equivalent to int('1'*number_of_haplotypes,2)
        
        for (pos, ind, read_id, base) in obs_tab:
            chr_id, pos2, ref, alt = leg_tab[ind]
            if pos!=pos2:
                raise Exception('error: the line numbers in obs_tab refer to the wrong chromosome positions in leg_tab.')
            if base==alt:
                hap_dict[(pos,base)] = hap_tab[ind]
            elif base==ref:
                hap_dict[(pos,base)] = hap_tab[ind] ^ b ### ^b flips all bits of the binary number, hap_tab[ind] using bitwise xor operator. 
            else:
                mismatches += 1
    
        print('%.2f%% of the reads matched known alleles.' % (100*(1-mismatches/len(obs_tab))))
    
        return hap_dict

    def intrenal_hap_dict(self, *alleles):
        """ This function allows treatment of alleles and haplotypes on an
        equal footing. This is done in three steps: (1) All the alleles and
        haplotypes are enumerated. (2) For each given haplotype, tuples in the
        reference panel, associated with the haplotype's alleles, are
        intersected. (3) A dictionary that lists all the alleles and haplotypes
        by their index is returned. The dictionary gives for each allele and
        haplotype their associated tuple and intersected tuple, respectively. """
        hap = dict()

        for i, X in enumerate(alleles):
            if type(X[0])==tuple: #Checks if X is a tuple/list of alleles.
                n = len(X)
                if n==1:
                    hap[1 << i] = self.hap_dict[X[0]]
                elif n==2:
                    hap[1 << i] = self.hap_dict[X[0]] & self.hap_dict[X[1]]
                else:
                    hap[1 << i] = reduce(and_,itemgetter(*X)(self.hap_dict))

            elif type(X[0])==int: #Checks if X is a single allele.
                hap[1 << i] = self.hap_dict[X]
            else:
                raise Exception('error: joint_frequencies only accepts alleles and tuple/list of alleles.')

        return hap

    def joint_frequencies_combo(self, *alleles, normalize):
        """ Based on the reference panel, it calculates joint frequencies of
            observed alleles. The function arguments are alleles, that is, 
            tuples of position and base, e.g., (100,'T'), (123, 'A') and 
            (386, 'C'). Each allele is enumerated according to the order it
            was received by the function. The function returns a dictionary that
            lists all the possible subgroups of the given alleles. Each key in
            the dictionary is a tuple of intergers that are lexicographically
            sorted. Moreover, each integer within the keys corresponds to an
            enumerated allele. For each subgroup of alleles the dictionary
            gives the joint frequencies in a given population. The function
            arguments can also include haplotypes, that is, tuples of alleles;
            Haplotypes are treated in the same manner as alleles. """

        hap = self.intrenal_hap_dict(*alleles)

        result = {c: popcount(A) for c,A in hap.items() }

        for C in combinations(hap, 2):
            result[C[0]|C[1]] = popcount(hap[C[0]]&hap[C[1]])

        for C in combinations(hap, 3):
            result[C[0]|C[1]|C[2]] = popcount(hap[C[0]]&hap[C[1]]&hap[C[2]])

        for r in range(4,len(alleles)):
            for C in combinations(hap, r):
                result[sum(C)] = popcount(reduce(and_,itemgetter(*C)(hap)))

        if len(alleles)>=4:
            result[sum(hap.keys())] = popcount(reduce(and_,hap.values()))

        if normalize:
            result = {k: v/self.total_number_of_haplotypes_in_reference_panel for k,v in result.items()}
            
        return result
    
    def likelihoods(self, *alleles):
        """ Calculates the likelihood to observe a set with alleles
        and haplotypes under four scenarios, namely, monosomy, disomy, SPH
        and BPH. """
        
        F = self.joint_frequencies_combo(*alleles, normalize=False)
        N = self.total_number_of_haplotypes_in_reference_panel #Divide values by N to normalize the joint frequencies.
        models = self.models_dict[len(alleles)]

        ### BPH ###
        (((A0, A1),((B0,),)),) = models['BPH'][1].items()
        BPH = F[B0] * A0 / ( A1 * N )

        BPH += sum( sum(F[B0] * F[B1] for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in models['BPH'][2].items()) / N**2

        if len(alleles)>2:
            BPH += sum( sum(F[B0] * sum(F[B1] * F[B2] for (B1, B2) in C[B0]) for B0 in C) * A0 / A1
                       for (A0, A1), C in models['BPH'][3].items()) / N**3

        ### SPH ###
        (((A0, A1),((B0,),)),) = models['SPH'][1].items()
        SPH = F[B0] * A0 / ( A1 * N ) 

        SPH += sum( sum(F[B0] * F[B1] for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in models['SPH'][2].items()) / N**2

        
        ### DIPLOIDY ###
        (((A0, A1),((B0,),)),) = models['DISOMY'][1].items()
        DISOMY = F[B0] * A0 / ( A1 * N ) 

        DISOMY += sum( sum(F[B0] * F[B1] for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in models['DISOMY'][2].items()) / N**2

        ### MONOSOMY ###
        ((B0,),) = models['MONOSOMY'][1][(1,1)]
        MONOSOMY = F[B0] / N 
        #MONOSOMY = F[int(N*'1',2)] / N 

        result = (MONOSOMY, DISOMY, SPH, BPH)
        return result
    
    def likelihoods2(self, *alleles):
        """ Calculates the likelihood to observe two alleles/haplotypes
        under four scenarios, namely, monosomy, disomy, SPH and BPH. """
        
        F = self.joint_frequencies_combo(*alleles, normalize=True)
        a, b, ab = F[1], F[2], F[3]
        BPH = (ab+2*a*b)/3 #The likelihood of three unmatched haplotypes.
        SPH = (5*ab+4*a*b)/9 #The likelihood of two identical haplotypes out three.
        DISOMY = (ab+a*b)/2 #The likelihood of diploidy.
        MONOSOMY = ab #The likelihood of monosomy.
        return MONOSOMY, DISOMY, SPH, BPH

    def likelihoods3(self, *alleles):
        """ Calculates the likelihood to observe three alleles/haplotypes
        under four scenarios, namely, monosomy, disomy, SPH and BPH. """
        
        F = self.joint_frequencies_combo(*alleles, normalize=True)
        a, b, c, ab, ac, bc, abc = F[1], F[2], F[4], F[3], F[5], F[6], F[7]
        BPH = (abc+2*(ab*c+ac*b+bc*a+a*b*c))/9 #The likelihood of three unmatched haplotypes.
        SPH = abc/3+2*(ab*c+ac*b+bc*a)/9  #The likelihood of two identical haplotypes out three.
        DISOMY = (abc+ab*c+ac*b+bc*a)/4 #The likelihood of diploidy.
        MONOSOMY = abc #The likelihood of monosomy.
        return MONOSOMY, DISOMY, SPH, BPH
    
    def likelihoods4(self, *alleles):
        """ Calculates the likelihood to observe four alleles/haplotypes
        under four scenarios, namely, monosomy, disomy, SPH and BPH. """
        
        F = self.joint_frequencies_combo(*alleles, normalize=True)
        a, b, c, d = F[1], F[2], F[4], F[8],
        ab, ac, ad, bc, bd, cd = F[3], F[5], F[9], F[6], F[10], F[12],
        abc, abd, acd, bcd = F[7], F[11], F[13], F[14]
        abcd = F[15]
        BPH = (abcd+2*(ab*c*d+a*bd*c+a*bc*d+ac*b*d+a*b*cd+ad*b*c+abc*d+a*bcd+acd*b+abd*c+ab*cd+ad*bc+ac*bd))/27  #The likelihood of three unmatched haplotypes.
        SPH = (17*abcd+10*(abc*d+bcd*a+acd*b+abd*c)+8*(ab*cd+ad*bc+ac*bd))/81  #The likelihood of two identical haplotypes out three.
        DISOMY = (abcd+abc*d+bcd*a+acd*b+abd*c+ab*cd+ad*bc+ac*bd)/8 #The likelihood of diploidy.
        MONOSOMY = abcd #The likelihood of monosomy.
        return MONOSOMY, DISOMY, SPH, BPH
        
def wrapper_of_examine_for_debugging(obs_filename,leg_filename,hap_filename,models_filename):
    """ Wrapper function of the class examine_homogeneous. It receives an observations
    file, IMPUTE2 legend file, IMPUTE2 haplotypes file, and a file with four
    statistical models. Based on the given data it creates and returns an
    instance of the class. """

    from MAKE_OBS_TAB import read_impute2

    if not os.path.isfile(obs_filename): raise Exception('Error: OBS file does not exist.')
    if not os.path.isfile(leg_filename): raise Exception('Error: LEGEND file does not exist.')
    if not os.path.isfile(hap_filename): raise Exception('Error: HAP file does not exist.')
    if not os.path.isfile(models_filename): raise Exception('Error: MODELS file does not exist.')

    leg_tab = read_impute2(leg_filename, filetype='leg')
    hap_tab, number_of_haplotypes = read_impute2(hap_filename, filetype='hap')
    with open(obs_filename, 'rb') as f:
        obs_tab = pickle.load(f)
        #info = pickle.load(f)

    load_model = bz2.BZ2File if models_filename[-4:]=='pbz2' else open
    with load_model(models_filename, 'rb') as f:
        models_dict = pickle.load(f)

    return examine_homogeneous(obs_tab, leg_tab, hap_tab, None, models_dict, number_of_haplotypes)

if __name__ != "__main__":
    print('The module LIKELIHOODS_CALCULATOR was imported.')
else:
    print('The module LIKELIHOODS_CALCULATOR was invoked directly')
    sys.exit(0)

###############################   END OF FILE   ###############################



"""
if __name__ != "__main__":
    print("The module LIKELIHOODS_CALCULATOR was imported.")
else:
    print("Executed when invoked directly")
    #sys.exit(0)
    import time, random
    a = time.time()
    obs_filename = 'results_EUR/mixed3haploids.X0.50.HG02035B.HG00451A.HG02513A.chr21.recomb.1.00.obs.p'
    hap_filename = '../build_reference_panel/EUR_panel.hg38.BCFtools/chr21_EUR_panel.hap.gz'
    leg_filename = '../build_reference_panel/EUR_panel.hg38.BCFtools/chr21_EUR_panel.legend.gz'
    models_filename = 'MODELS/MODELS16.p'
    
    A = wrapper_of_examine_for_debugging(obs_filename,leg_filename,hap_filename,models_filename)

    alleles = tuple(A.hap_dict.keys())

    frequencies = lambda *x: {bin(a)[2:]:b for a,b in A.joint_frequencies_combo(*x,normalize=True).items()}

    likelihoods = A.likelihoods
    likelihoods2 = A.likelihoods2
    likelihoods3 = A.likelihoods3
    likelihoods4 = A.likelihoods4

    x = random.randrange(len(alleles)) #123
    haplotypes = (alleles[x:x+4],alleles[x+4:x+8],alleles[x+8:x+12],alleles[x+12:x+16])

    print('-----joint_frequencies_combo-----')    
    print(frequencies(alleles[x+0]))
    print(frequencies(alleles[x:x+4]))
    print('-----likelihoods4-haplotypes-----')
    print(haplotypes)
    print(frequencies(*haplotypes))
    print(likelihoods(*haplotypes))
    print(likelihoods4(*haplotypes))
    print('-----likelihoods2-----')
    print(alleles[x:x+2])
    print(frequencies(*alleles[x:x+2]))
    print(likelihoods(*alleles[x:x+2]))
    print(likelihoods2(*alleles[x:x+2]))
    print('-----likelihoods3-----')
    print(alleles[x:x+3])
    print(frequencies(*alleles[x:x+3]))
    print(likelihoods(*alleles[x:x+3]))
    print(likelihoods3(*alleles[x:x+3]))
    print('-----likelihoods4-----')
    print(alleles[x:x+4])
    print(frequencies(*alleles[x:x+4]))
    print(likelihoods(*alleles[x:x+4]))
    print(likelihoods4(*alleles[x:x+4]))

    b = time.time()

    print('Done in %.3f sec.' % ((b-a)))

"""