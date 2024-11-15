#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BALANCED_ROC_CURVES

Based on simulated data that was created by MIX_HAPLOIDS and analyzed by
ANEUPLOIDY_TEST, balanced ROC (Receiver Operating Characteristic) curves
for predicting BPH (both parental homologs) and SPH (single parental homologs)
are created.

The balanced ROC curve is a plot of BTPR (Balanced True Positive Rate) vs. BFPR
(Balanced False Positive Rate), where BTPR=0.5*[TPR(BPH) + TPR(SPH)] and B
FPR=0.5*[FPR(BPH) + FPR(SPH)].

The genome is divided into bins. For each simulated data, the mean LLR and the
standard deviation for a bin are calculated. A positive (negative) prediction
is made if the bounds of the confidence interval lie on the positive (negative)
side of the number line. In each bin, the z-score is varied and the balanced
positive and negative rates are calculated for each value it takes.

Daniel Ariad (daniel@ariad.org)
Dec 30, 2021
"""

import pickle, bz2, gzip, sys, os, argparse
from math import log
from operator import attrgetter
from collections import defaultdict, namedtuple

leg_tuple = namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
sam_tuple = namedtuple('sam_tuple', ('sample_id', 'group1', 'group2', 'sex')) #Encodes the rows of the samples table
obs_tuple = namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table

def LLR(y,x):
    """ Calculates the logarithm of y over x and deals with edge cases. """
    if x and y:
        result = log(y/x)
    elif x and not y:
        result = -1.23456789
    elif not x and y:
        result = +1.23456789
    elif not x and not y:
        result = 0
    else:
        result = None
    return result

def chr_length(chr_id):
    """ Return the chromosome length for a given chromosome, based on the reference genome hg38."""
    #The data of chromosome length was taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
    length_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                  'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
                  'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                  'chr16': 90338345, 'chr17':  83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                  'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
    return length_dict[chr_id]

def mean_and_std_of_mean_of_rnd_var(A):
    """ Calculates the mean and population standard deviation of the mean of random variables.
        Each row of A represents a random variable, with observations in the columns."""
    if type(A)==dict:
        A = tuple(tuple(i) for i in A.values()) 
    
    M, N = len(A), len(A[0])
    mu = sum(sum(likelihoods_in_window)/N for likelihoods_in_window in A)
    arg = ((sum(sampled_likelihoods) - mu)**2 for sampled_likelihoods in zip(*A))
    std = (sum(arg) / (N - 1))**.5 / M
    mean = mu / M
    return mean, std

def load_likelihoods(filename):
    """ Loads from a file a dictionary that lists genomic windows that contain
    at least two reads and gives the bootstrap distribution of the
    log-likelihood ratios (LLRs). """
    
    Open = {'bz2': bz2.open, 'gzip': gzip.open}.get(filename.rpartition('.')[-1], open)
    try:
        with Open(filename, 'rb') as f:
            likelihoods = pickle.load(f)
            info = pickle.load(f)
        return likelihoods, info
    except Exception as e:
        print(filename,e)
        return None


def show_info(filename, info, pairs):
    S = info['statistics']
    ancestral_makeup = ", ".join("{:.1f}% {}".format(100*v, k) for k, v in info['ancestral_makeup'].items()) if type(info['ancestral_makeup'])==dict else ', '.join(info['ancestral_makeup'])
    print('\nFilename: %s' % filename)
    print('\nSummary statistics')
    print('------------------')    
    print('Chromosome ID: %s, Depth: %.2f.' % (info['chr_id'],info['depth']))
    print('Number of genomic windows: %d, Mean and standard error of genomic window size: %d, %d.' % (S.get('num_of_windows',0),S.get('window_size_mean',0),S.get('window_size_std',0)))
    print('Mean and standard error of meaningful reads per genomic window: %.1f, %.1f.' % (S.get('reads_mean',0), S.get('reads_std',0)))
    print('Ancestral makeup: %s, Fraction of alleles matched to the reference panel: %.3f.' % (ancestral_makeup, info['statistics']['matched_alleles']))

    for pair in pairs:
        if 'LLRs_per_chromosome' in S and tuple(pair) in S['LLRs_per_chromosome']:
            L = S['LLRs_per_chromosome'][tuple(pair)]
            print(f"--- Chromosome-wide LLR between {pair[0]:s} and {pair[1]:s} ----")
            print(f"Mean LLR: {L['mean_of_mean']:.3f}, Standard error of the mean LLR: {L['std_of_mean']:.3f}")
            print(f"Fraction of genomic windows with a negative LLR: {L['fraction_of_negative_LLRs']:.3f}")
        

def bin_genomic_windows(windows,chr_id,num_of_bins):
    """ Lists the bins and gives the genomic windows that they contain. """
    bin_size = chr_length(chr_id) / num_of_bins
    result = {}
    j = 0
    
    for i in range(num_of_bins): ### All bins before the first the genomic window are filled with Nones.
        if sum(windows[0])/2 < (i+1)*bin_size:
            break
        result[i/num_of_bins,(i+1)/num_of_bins] = None
    
    for k,(a,b) in enumerate(windows):
        if not bin_size*i <= (a+b)/2 < bin_size*(i+1):
            result[i/num_of_bins,(i+1)/num_of_bins] = (j,k)
            j = k
            for i in range(i+1,num_of_bins): #Proceed to the next non-empty bin; Empty bins are filled with Nones.
                if (a+b)/2 < (i+1)*bin_size:
                    break
                result[i/num_of_bins,(i+1)/num_of_bins] = None
    
    for i in range(i,num_of_bins): ### All bins after the last the genomic window are filled with Nones.
        result[i/num_of_bins,(i+1)/num_of_bins] = (j,k) if j != k else None
        j = k 
    return result

def binning(LLRs_per_window,info,num_of_bins):
    """ Genomic windows are distributed into bins. The LLRs in a genomic windows
    are regarded as samples of a random variable. Within each bin, we calculate
    the mean and population standard deviation of the mean of random variables. 
    The boundaries of the bins as well as the mean LLR and the standard-error
    per bin are returned. """
             
    #K,M,V = tuple(LLR_stat.keys()), *zip(*LLR_stat.values())
    list_of_windows = [*LLRs_per_window.keys()]
    bins = bin_genomic_windows(list_of_windows, info['chr_id'], num_of_bins)
    X = [*bins]
    
    LLR_matrix = [*LLRs_per_window.values()]
    Y, E = [], []
    for C in bins.values():
        if C:
            mean, std = mean_and_std_of_mean_of_rnd_var(LLR_matrix[C[0]:C[1]])
        else:
            mean, std = None, None
        
        Y.append(mean)
        E.append(std)
    
    return X,Y,E

def collect_data(criteria, num_of_bins, work_dir):
    """ Iterates over all the data files in the folder and creates a dictionary
    the list all the files that fit the criteria and gives their analysis
    (via the bucketing function). """

    import glob
    
    scenarios = criteria['scenarios']
    del criteria['scenarios']
    
    filenames = glob.glob(work_dir.rstrip('/') +'/'+ "simulated*LLR.p*")
    result = {scenarios[0]: defaultdict(dict), scenarios[1]: defaultdict(dict)}
    for filename in filenames:
        try:
            likelihoods, info = load_likelihoods(filename)
            if likelihoods==None: continue
            subinfo = {x: info.get(x,None) for x in criteria.keys()}
            scenario = filename.split('.',2)[1]
            
            ###show_info(filename, info, scenarios)
            ###print(subinfo)
            ###print(criteria)
            ###print(criteria==subinfo)
            ###print(scenario, scenarios)
            if subinfo==criteria and scenario in scenarios:
                chr_id = info['chr_id']
                                                                    
                LLRs = {window: tuple(LLR(attrgetter(scenarios[0])(l), attrgetter(scenarios[1])(l)) for l in likelihoods_in_window)
                                   for window,likelihoods_in_window in likelihoods.items()} 
                
                result[scenario][chr_id][filename] = {bucket: {'mean': mean, 'std': std}
                                                   for (bucket,mean,std) in zip(*binning(LLRs,info,num_of_bins))}
        except Exception as err:
            print(err)
    return result

def prediction_rates(data, thresholds, positive='both'):
    """ Creates a nested dictionary that lists chromosomes, buckets and
    thresholds. For each entry the dictionary gives the false and true positive
    rates. """

    label_B, label_A = data.keys()
    B, A = data.values()
    
    prediction_rates = {}
    for chr_id in A:
        prediction_rates[chr_id] = {}
        some_filename = next(iter(A[chr_id]))
        ###print(len(A[chr_id]),len(B[chr_id]))
        for bucket in A[chr_id][some_filename]:
            prediction_rates[chr_id][bucket] = {}
            for z in thresholds:
                true_B = [file[bucket]['mean'] > z * file[bucket]['std'] for file in B[chr_id].values() if file[bucket]['mean']!=None]
                false_B = [file[bucket]['mean'] > z * file[bucket]['std'] for file in A[chr_id].values() if file[bucket]['mean']!=None]

                false_A = [file[bucket]['mean'] < -z * file[bucket]['std'] for file in B[chr_id].values() if file[bucket]['mean']!=None]
                true_A = [file[bucket]['mean'] < -z * file[bucket]['std'] for file in A[chr_id].values() if file[bucket]['mean']!=None]

                if true_B==[] or false_B==[]:
                    break
                elif positive == 'both':
                    TPR = 0.5 * (true_B.count(1)/len(true_B) + true_A.count(1)/len(true_A))
                    FPR = 0.5 * (false_B.count(1)/len(false_B) + false_A.count(1)/len(false_A))
                elif positive == label_A:
                    TPR = true_A.count(1)/len(true_A)
                    FPR = false_A.count(1)/len(false_A)
                elif positive == label_B:
                    TPR = true_B.count(1)/len(true_B)
                    FPR = false_B.count(1)/len(false_B)
                else:
                    break

                prediction_rates[chr_id][bucket][z] = (FPR,TPR) ### <--- Structure of the nested dictionary

    return prediction_rates

def main(work_dir,output_filename,number_of_bins,criteria,compress):
    """ Creates Balanced ROC curves for predicting aneuploidy. """
    
    assert os.path.isdir(work_dir), 'The path to the directory that contains simulated data does not exist.'
    Z = [i/33 for i in [*range(-1800,-300,300)]+[*range(-300,300)]+[*range(300,1800,300)]]
    data = collect_data(criteria.copy(), number_of_bins, work_dir)
    R = prediction_rates(data, thresholds = Z)
    N = {scenario: {chr_id: len(data[scenario][chr_id]) for chr_id in data[scenario]} for scenario in data}
    
    Open = {'bz2': bz2.open, 'gz': gzip.open}.get(compress, open)
    ext = ('.'+compress) * (compress in ('bz2','gz'))
    with Open(output_filename + ext, "wb") as f:
        pickle.dump(criteria, f, protocol=4) #Simulation info
        pickle.dump(N, f, protocol=4) #Statistical info
        pickle.dump(R, f, protocol=4) #BTPR vs. BFPR for each bin along the genome
        
    return N, R

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description='Based on simulated data that was created by MIX_HAPLOIDS and analyzed by'
                'DETECT_CROSSOVERS, Balanced ROC curves for predicting aneuploidy are created as the z-score varies.')
    parser.add_argument('path_to_simulated_data', type=str, metavar='PATH_TO_SIMULATED_DATA', 
                        help='Path of a directory that contains LLR files with simulated scenarios.')
    parser.add_argument('output_filename', type=str, metavar='OUTPUT_FILENAME',
                        help='The output filename.')
    parser.add_argument('-n', '--number-of-bins', type=int, metavar='INT',  default='15',
                        help='The genome is divided into bins and for each bin a ROC curve is calculated. Default value is 15.')
    parser.add_argument('-c', '--compress', metavar='gz/bz2/unc', type=str, default='unc',  choices=['gz','bz2','unc'],
                        help='Output compressed via gzip, bzip2 or uncompressed. Default is uncompressed.')
    parser.add_argument('-o', '--scenarios', type=str, nargs=2,
                        metavar='BPH/SPH/disomy/monosomy', default=['BPH', 'SPH'], choices=['BPH','SPH','disomy','monosomy'],
                        help="Two simulated scenarios for which a balanced ROC curve would be created. The default is \"BPH SPH\". ")
    parser.add_argument('-a', '--ancestral-makeup', metavar='STR', nargs='+',
                        help='Apply a criterion for the ancestral makeup: \n'
                             'a. For non-admixtures the argument consists a single superpopulation, e.g., EUR. \n'
                             'b. For recent admixtures the argument consists two superpopulations, e.g., EUR EAS. \n'
                             'c. For distant admixtures the argument consists of the superpoplations and their proportions, e.g, EUR 0.8 EAS 0.1 SAS 0.1. \n')
    parser.add_argument('-i', '--chr-id', type=str, metavar='STR',
                        help='Apply a criterion for the chromosome number, e.g., chrX.')
    parser.add_argument('-w', '--window-size', type=int, metavar='INT', 
                        help='Apply a criterion for the size of the genomic window.')
    parser.add_argument('-s', '--subsamples', type=int, metavar='INT', 
                        help='Apply a criterion for the number of subsamples per genomic window.')
    parser.add_argument('-m', '--min-reads', type=int, metavar='INT',
                        help='Apply a criterion for the minimal number of reads per homolog, admitting non-zero score.')
    parser.add_argument('-M', '--max-reads', type=int, metavar='INT',
                        help='Apply a criterion for the maximal number of sampled reads per bootstrap iteration.')
    parser.add_argument('-f', '--min-HF', type=float, metavar='FLOAT',
                        help='Apply a criterion for the minimal haplotype frequency.')
    parser.add_argument('-S', '--min-score', type=int, metavar='INT',
                        help='Apply a criterion for the minimal score of a read.')
    parser.add_argument('-l', '--read-length', type=int, metavar='INT',  
                        help='Apply a criterion for the number of base pairs in read.')
    parser.add_argument('-d', '--depth', type=float, metavar='FLOAT',  
                        help='Apply a criterion for the depth of coverage.')
    
    args = vars(parser.parse_args())
    
    #print(args)
    criteria = {}
    criteria['scenarios'] =  args['scenarios']
    if args['depth']!=None: criteria['depth'] = args['depth']
    if args['window_size']!=None: criteria['window_size'] = args['window_size']
    if args['min_reads']!=None: criteria['min_reads'] = args['min_reads']
    if args['max_reads']!=None: criteria['max_reads'] = args['max_reads']
    if args['min_score']!=None: criteria['min_score'] = args['min_score']
    if args['min_HF']!=None: criteria['min_HF'] = args['min_HF']
    if args['chr_id']!=None: criteria['chr_id'] = args['chr_id']
    if args['read_length']!=None: criteria['read_length'] = args['read_length']
   
    if args['ancestral_makeup']!=None:
        strings_even = all(i.isalpha() for i in args['ancestral_makeup'][0::2])
        strings_odd = all(i.isalpha() for i in args['ancestral_makeup'][1::2])
        floats_odd = all(i.replace('.','',1).isdigit() for i in args['ancestral_makeup'][1::2])
 
        if strings_even and strings_odd and len(args['ancestral_makeup']) in (1,2):
            criteria['ancestry'] = set(args['ancestral_makeup'])
        elif  strings_even and floats_odd:
            criteria['ancestry'] = dict(zip(args['ancestral_makeup'][0::2],map(float,args['ancestral_makeup'][1::2])))
        else:
            print('caution: unsupported argument was supplied with --ancestral-makeup. Thus, this criterion would be ignored.')
    
    print('- Path to simulated data:',args['path_to_simulated_data'],'\n- Output filename:',args['output_filename'],'\n- Number of bins:',args['number_of_bins'])
    print('- Criteria:')
    print("\n".join("\t{}. {}:\t{}".format(i, k, v) for i, (k, v) in enumerate(criteria.items(),start=1)))
    N,R = main(args['path_to_simulated_data'],args['output_filename'],args['number_of_bins'],criteria,args['compress'])
    print('- Simulated files that fulfilled the criteria:\n',N)
    sys.exit(0)
else:
    print("The module BALANCED_ROC_CURVE was imported.")

### END OF FILE ###
