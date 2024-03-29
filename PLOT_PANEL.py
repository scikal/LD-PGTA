#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PLOT_PANEL

Plots log-likelihood ratio vs. chromosomal position from a LLR file.

May 20, 2020
"""

import pickle, bz2, gzip, collections, math
from statistics import mean, variance, median
from math import log
from operator import attrgetter, itemgetter
from itertools import accumulate
import argparse, sys

likelihoods_tuple = collections.namedtuple('likelihoods_tuple', ('monosomy', 'disomy', 'SPH', 'BPH'))

def chr_length(chr_id):
    """ Return the chromosome length for a given chromosome, based on the reference genome hg38.""" 
    #The data of chromosome length was taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
    
    length_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                  'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
                  'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                  'chr16': 90338345, 'chr17':  83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                  'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
        
    return length_dict[f'chr{chr_id:s}' if chr_id.isdigit() else chr_id]

def mean_and_var(x):
    """ Calculates the mean and variance. """
    cache = tuple(x)
    m = mean(cache)
    var = variance(cache, xbar=m)
    return m, var

def mean_and_std_of_mean_of_rnd_var(A):
    """ Calculates the mean and sample standard deviation of the mean of random variables.
        Each row of A represents a random variable, with observations in the columns."""
    if type(A)==dict:
        A = tuple(tuple(i) for i in A.values()) 
    
    M, N = len(A), len(A[0]) ### M is the number of random variables, while N is the number of samples.
    mu = sum(sum(likelihoods_in_window) for likelihoods_in_window in A) / N
    arg = ((sum(sampled_likelihoods) - mu)**2 for sampled_likelihoods in zip(*A))
    std = (sum(arg) / (N - 1))**.5 / M
    mean = mu / M
    return mean, std

def median_and_mad_of_mean_of_rnd_var(A):
    """ Calculates the median of GW's mean and median absolute deviation of the GW's mean.
        Each row of A corresponds to a random variable/genomic window, while each column 
        corresponds to a different observation/sample. """
    if type(A)==dict:
        A = tuple(tuple(i) for i in A.values()) 
    
    M, N = len(A), len(A[0]) ### M is the number of random variables, while N is the number of samples.
    m = median(sum(sampled_likelihoods) for sampled_likelihoods in zip(*A))
    mad = median(abs(sum(sampled_likelihoods) - m) for sampled_likelihoods in zip(*A)) 
    mu = m / M
    adjust_mad =  mad / ( 0.6744897501960817 * M) 
    return mu, adjust_mad

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

def load_likelihoods(filename):
    """ Loads from a file a dictionary that lists genomic windows that contain
    at least two reads and gives the bootstrap distribution of the 
    log-likelihood ratios (LLRs). """
    
    Open = {'bz2': bz2.open, 'gzip': gzip.open}.get(filename.rpartition('.')[-1], open)
    
    with Open(filename, 'rb') as f:
        likelihoods = pickle.load(f)
        info = pickle.load(f)
    return likelihoods, info

def show_info(filename, info, pairs=(('BPH','SPH'),)):
    S = info['statistics']
    ancestral_makeup = ", ".join("{:.1f}% {}".format(100*v, k) for k, v in info['ancestral_makeup'].items()) if type(info['ancestral_makeup'])==dict else ', '.join(info['ancestral_makeup'])
    matched_alleles = ", ".join("{}: {:.1f}%".format(k,100*v) for k, v in info['statistics']['matched_alleles'].items())

    print('\nFilename: %s' % filename)
    print('\nSummary statistics')
    print('------------------')    
    print('Chromosome ID: %s, Depth: %.2f.' % (info['chr_id'],info['depth']))
    print('Number of genomic windows: %d, Mean and standard error of genomic window size: %d, %d.' % (S.get('num_of_windows',0),S.get('window_size_mean',0),S.get('window_size_std',0)))
    print('Mean and standard error of meaningful reads per genomic window: %.1f, %.1f.' % (S.get('reads_mean',0), S.get('reads_std',0)))
    print('Ancestral makeup: %s, Fraction of alleles matched to the reference panel: %s.' % (ancestral_makeup, matched_alleles))

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
            #####mean, std = median_and_mad_of_mean_of_rnd_var(LLR_matrix[C[0]:C[1]])
        else:
            mean, std = None, None
        
        Y.append(mean)
        E.append(std)
    
    return X,Y,E

def detect_crossovers_v2(genomic_windows, mean_of_LLRs, variance_of_LLRs, z_score=1.96, lookahead=20):
    """ Detecting crossovers by indetifying transitions between BPH and SPH 
        regions. """
        
    crossovers = {}
    #Scan the chromosome in the 5'-to-3' direction to find crossovers.
    x_coord = tuple(0.5*(a+b) for a,b in genomic_windows)
    acc_means = tuple(accumulate(mean_of_LLRs))
    acc_vars = tuple(accumulate(variance_of_LLRs))
    triple = tuple(zip(x_coord, acc_means, acc_vars))
    
    #Maxima and minima candidates are temporarily stored in mx and mn, respectively.
    mx, mn, last_ind, last_extremum = None, None, 0, 0
        
    def recover_skipped_extremum(extremum_type,last_ind,ind):
        """ Recovers skipped extremum """
        ### extremum_typ is the max (minimum) function if a maximum (minimum) was skipped.
        ### last_ind is the index of the last detected extremum.
        ### ind is the index of the most recent detected extremum.
        M0, M2, V0, V2 = acc_means[last_ind],acc_means[ind], acc_vars[last_ind], acc_vars[ind]
        Z1, X1, M1, V1 = extremum_type((((M1-M0)/(V1-V0)**.5-(M2-M1)/(V2-V1)**.5, X1, M1, V1) 
                              for X1,M1,V1 in triple[last_ind+5:ind-4] if V0!=V1!=V2), key=itemgetter(0))
        kappa = min(abs(M2-M1)/(V2-V1)**.5,abs(M0-M1)/(V1-V0)**.5)
        print(f'Recovering skipped {extremum_type.__name__:s} point:',(Z1, X1, M1, V1))
        return {X1: kappa}
        
    for index, (x, y, v) in enumerate(triple):
        
        if  mx==None or y > mx:
            mx_index,mx_pos,mx,mx_var = index,x,y,v
        
        if  mn==None or y < mn:
            mn_index,mn_pos,mn,mn_var = index,x,y,v

        if mx!=None and 0 < (mx-y)-z_score*(v-mx_var)**.5 and index-mx_index>=lookahead: #maximum point detected
            for x2, y2, v2 in triple[max(mx_index-lookahead,last_ind):last_ind:-1]:
                if 0 < (mx-y2)-z_score*(mx_var-v2)**.5:
                   
                    if last_extremum == +1: #Recovering skipped minimum point.
                        crossovers.update(recover_skipped_extremum(min,last_ind,mx_index)) 

                    kappa = min((mx-y2)/(mx_var-v2)**.5,(mx-y)/(v-mx_var)**.5)
                    crossovers[mx_pos] = kappa
                    mx, mn, last_ind, last_extremum = None, None, mx_index, +1 # set algorithm to find the next minimum
                    break
                
        if mn!=None and 0 < (y-mn)-z_score*(v-mn_var)**.5  and index-mn_index>=lookahead: #minimum point detected
            for x2, y2, v2 in triple[max(mn_index-lookahead,last_ind):last_ind:-1]:
                if  0 < (y2-mn)-z_score*(mn_var-v2)**.5:
                    
                    if last_extremum == -1: #Recovering skipped maximum point.
                        crossovers.update(recover_skipped_extremum(max,last_ind,mn_index))
                        
                    kappa = min((y2-mn)/(mn_var-v2)**.5,(y-mn)/(v-mn_var)**.5)
                    crossovers[mn_pos] = kappa
                    
                    mx, mn, last_ind, last_extremum = None, None, mn_index, -1 # set algorithm to find the next maxima
                    break
                
    return crossovers

def capitalize(x):
    return x[0].upper() + x[1:]
    
def panel_plot(DATA,**kwargs):
    """ Creates a multi-panel figure. For each numbered chromosome, a figure 
        depicts the log-likelihood ratio vs. chromosomal position for BPH over
        SPH. """
    
    import matplotlib as mpl
    mpl.rcParams.update({'figure.max_open_warning': 0})
    
    
    scale = kwargs.get('scale', 0.5)
    save = kwargs.get('save', '')
    z_score = kwargs.get('z_score', 1.96)
    lookahead = kwargs.get('lookahead', 30)

    bin_size = kwargs.get('bin_size', 4000000)
    pairs = kwargs.get('pairs', (('BPH','disomy'),('disomy','SPH'),('SPH','monosomy')))
    
    fs=28 * scale
    columns = 6
    rows = math.ceil(len(DATA)/columns)
        
    if save!='':
            mpl.use('Agg')
    else:
        #['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']
        mpl.use('Qt5Agg')
    
    
    import matplotlib.pyplot as plt
    num_of_bins = {info['chr_id']: chr_length(info['chr_id'])//bin_size for (likelihoods,info) in DATA.values()}

    colors = {frozenset(('BPH','disomy')):(177/255,122/255,162/255),
              frozenset(('disomy','SPH')):(242/255,142/255,44/255),
              frozenset(('SPH','monosomy')):(239/255,106/255,92/255),
              frozenset(('disomy','monosomy')):(104/255,162/255,183/255),
              frozenset(('BPH','SPH')):(104/255,162/255,104/255)}

    if len(DATA)>columns:
        fig,axs = plt.subplots(rows ,columns, sharex='col', sharey='row', figsize=(6.666 * columns * scale, 5.625 * rows * scale))
        fig.subplots_adjust(left=0.05, bottom=0.1, right=.99, top=(0.92 if kwargs.get('title',None) else 0.96), wspace=None, hspace=None)
    else:
        fig,axs = plt.subplots(rows ,columns, sharex='none', sharey='row', figsize=( 6.666 * columns * scale, 1.25 * 5.625 * rows * scale))
        fig.subplots_adjust(left=0.05, bottom=0.3, right=.99, top=(0.82 if kwargs.get('title',None) else 0.86), wspace=None, hspace=None)
    
    
    AX = [i for j in axs for i in j] if len(DATA)>columns else axs
    
    H = {}
    YMAX = [0]*len(DATA)
    crossovers = {}
    for a,b in pairs:
        for g,(ax1,(likelihoods,info)) in enumerate(zip(AX,DATA.values())):
    

            LLRs = {window: tuple(LLR(attrgetter(a)(l), attrgetter(b)(l)) for l in likelihoods_in_window)
                               for window,likelihoods_in_window in likelihoods.items()} 
                                                               
            X,Y,E = binning(LLRs,info,num_of_bins[info['chr_id']])
            Y = [(y if y else 0) for y in Y]
            E = [(z_score*e if e else 0) for e in E]
            
            T = [(x[1]+x[0])/2 for x in X]                
            steps_x = [X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]]
            steps_y = [i for i in Y for j in (1,2)]
            H[a,b] = ax1.plot(steps_x, steps_y, label=f'{capitalize(a):s} vs. {capitalize(b):s}' ,color=colors[frozenset((a,b))], linewidth=2, zorder=10, scalex=True, scaley=True, alpha=0.8)
            
            P = [(x[1]-x[0])/2 for x in X]                
            ax1.errorbar(T, Y, xerr = P, ecolor=colors[frozenset((a,b))],marker=None, ls='none',alpha=1, zorder=13, linewidth=5*scale) 
            ax1.errorbar(T, Y, yerr = E, ecolor='black',marker=None, ls='none',alpha=0.2, zorder=15, linewidth=4*scale) 
            
            yabsmax = max(map(abs,Y))
            
            if pairs==(('BPH','SPH'),) or pairs==(('SPH','BPH'),):
                genomic_windows = info['statistics']['LLRs_per_genomic_window'][('BPH','SPH')]
                mean_of_LLRs = [*map(itemgetter(0),genomic_windows.values())]
                variance_of_LLRs = [*map(itemgetter(1),genomic_windows.values())]
                unnormalized_crossovers = detect_crossovers_v2(genomic_windows, mean_of_LLRs, variance_of_LLRs, z_score=z_score, lookahead=lookahead)
                l = chr_length(info['chr_id'])
                crossovers[g] = [pos/l for pos in unnormalized_crossovers] #Normalize position according to the chromosome length.
                
            YMAX[g] = yabsmax if YMAX[g]< yabsmax else YMAX[g]

    for g,(ax1,(identifier,(likelihoods,info))) in enumerate(zip(AX,DATA.items())):
        mean_genomic_window_size = info['statistics']['window_size_mean']/chr_length(info['chr_id']) 
        ymax = max(YMAX[columns*(g//columns):columns*(g//columns+1)])
        ax1.errorbar( 0.88-mean_genomic_window_size, -0.76*ymax,marker=None, ls='none', xerr=25*mean_genomic_window_size, linewidth=2*scale, color='k', capsize=4*scale, zorder=20)
        ax1.text(     0.88-mean_genomic_window_size, -0.82*ymax, '25 GW',  horizontalalignment='center', verticalalignment='top',fontsize=2*fs//3, zorder=20)
        ax1.plot([0,1],[0,0],color='black', ls='dotted',alpha=0.7,zorder=0, linewidth=2*scale, scalex=False, scaley=False)
        ax1.set_title(identifier,fontsize=fs)

    for g,ax1 in enumerate(AX[:len(DATA)]):
        ymax = max(YMAX[columns*(g//columns):columns*(g//columns+1)])
        ax1.set_ylim((-1.01*ymax,+1.01*ymax))
        ax1.set_xlim((0,1)) 
        
        #Replace ticks along the x-axis 
        X_ticks = [i/10 for i in range(0,11,2)]
        X_labels = [('%g' % j) for j in X_ticks] 
        ax1.set_xticks(X_ticks)
        ax1.set_xticklabels(X_labels)
        
        ax1.tick_params(axis='x', labelsize=fs) 
        ax1.tick_params(axis='y', labelsize=fs)
        ax1.xaxis.set_tick_params(width=2*scale)
        ax1.yaxis.set_tick_params(width=2*scale)
        ###ax1.grid(color='black', linestyle='-.', linewidth=1,alpha=0.5)
        for axis in ['top','bottom','left','right']:
            ax1.spines[axis].set_linewidth(2*scale)
            
        if pairs==(('BPH','SPH'),) or pairs==(('SPH','BPH'),):
            for i in crossovers[g]:
                ax1.plot([i,i],[-1.01*ymax,1.01*ymax],color='purple', ls='dotted',alpha=0.7,zorder=19, linewidth=2*scale, scalex=False, scaley=False)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    
    plt.xlabel('Chromosomal position', fontsize=fs,labelpad=23*scale)
    plt.ylabel('Log-likelihood ratio', fontsize=fs,labelpad=45*scale)        
    
    for l in range(1,len(AX)-len(DATA)+1):
        AX[-l].tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False, width=0)
        for axis in ['top','bottom','left','right']:
            AX[-l].spines[axis].set_visible(False)
        AX[-l].xaxis.set_tick_params(labelbottom=True)
    if len(H)>1:
        fig.legend(handles=[i[0] for i in H.values()], title='', loc='lower right', ncol=len(H), fancybox=True,fontsize=fs) # bbox_to_anchor=(.5, .45)
    
    if kwargs.get('title',None): 
        fig.suptitle(kwargs['title'], fontsize=int(1.2*fs), color='black', fontweight="bold") #, y=1.01
            
        
    if save!='':
        print('Saving plot...')
        #####plt.tight_layout()
        extension = 'svg'
        plt.savefig('.'.join([save,extension]), format=extension) # bbox_inches='tight'
        plt.close(fig)
    else:
       #plt.tight_layout()
       plt.show() 

def single_plot(likelihoods,info,**kwargs):
    """ Creates a figure  depicts the log-likelihood ratio vs. chromosomal
        position for (a) BPH over disomy, (b) disomy over SPH and (c) SPH over 
        monosomy. """
        
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    scale = kwargs.get('scale', 1)
    z_score = kwargs.get('z_score', 1.96)
    bin_size = kwargs.get('bin_size', 4000000)
    save = kwargs.get('save', '')
    pairs = kwargs.get('pairs', (('BPH','disomy'),('disomy','SPH'),('SPH','monosomy')))
    
    if save!='':
        mpl.use('Agg')
    else:
        #['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']
        mpl.use('Qt5Agg')
    
    
    num_of_bins = {info['chr_id']: chr_length(info['chr_id'])//bin_size}

    
    fs = 24 * scale
    
    colors = {frozenset(('BPH','disomy')):(177/255,122/255,162/255),
              frozenset(('disomy','SPH')):(242/255,142/255,44/255),
              frozenset(('SPH','monosomy')):(239/255,106/255,92/255),
              frozenset(('disomy','monosomy')):(104/255,162/255,183/255),
              frozenset(('BPH','SPH')):(104/255,162/255,104/255)}
    
    LLRs = {(i,j): 
            {window: tuple(LLR(attrgetter(i)(l), attrgetter(j)(l)) for l in likelihoods_in_window)
                       for window,likelihoods_in_window in likelihoods.items()} 
                                                        for i,j in pairs}
        
    fig,(ax1)=plt.subplots(1,1, figsize=(16 * scale, 9 * scale))
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    H = {}
    for p,LLRs_per_genomic_window in LLRs.items():
        X,Y,E = binning(LLRs_per_genomic_window,info,num_of_bins[info['chr_id']])
        Y = [(y if y else 0) for y in Y]
        E = [(z_score*e if e else 0) for e in E]        
        T = [(x[1]+x[0])/2 for x in X]            
        
        ###ax1.plot([X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]],[i for i in Y for j in (1,2)], label=f'{p[0]:s} vs. {p[1]:s}',color=colors[p],linewidth=2)
        
        steps_x = [X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]]
        steps_y = [i for i in Y for j in (1,2)]
        H[p] = ax1.plot(steps_x, steps_y, label=f'{capitalize(p[0]):s} vs. {capitalize(p[1]):s}' ,color=colors[frozenset(p)], linewidth=2*scale, zorder=10, scalex=True, scaley=True, alpha=0.8)
        P = [(x[1]-x[0])/2 for x in X]                
        ax1.errorbar(T, Y, xerr = P, color=colors[frozenset(p)],marker=None, ls='none',alpha=1, zorder=13, linewidth=3*scale) 
        
        
        ax1.errorbar(T, Y, yerr = E, ecolor='black',marker=None, ls='none',alpha=0.2, linewidth=scale, zorder=15) 

    
    ax1.tick_params(axis='x', labelsize=fs) 
    ax1.tick_params(axis='y', labelsize=fs)
    ax1.xaxis.set_tick_params(width=scale)
    ax1.yaxis.set_tick_params(width=scale)
    ###ax1.grid(color='black', linestyle='-.', linewidth=1,alpha=0.5)
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(scale)
        
    #ax1.set_title(info['chr_id'].replace('chr', 'Chromosome '),fontsize=fs)
    ax1.set_title(f"{info['chr_id'].replace('chr', 'Chromosome '):s}, {info['depth']:.2f}x", fontsize=int(1.2 * fs), color='black', fontweight="bold")

    ax1.set_ylabel('Log-likelihood ratio', fontsize=fs,labelpad=2*scale)   
    ax1.set_xlabel('Chromosomal position', fontsize=fs,labelpad=2*scale)

    
    #Replace ticks along the x-axis 
    X_ticks = [i/10 for i in range(0,11,2)]
    X_labels = [('%g' % j) for j in X_ticks] 
    ax1.set_xticks(X_ticks)
    ax1.set_xticklabels(X_labels)
    
    #Y_ticks = [i for i in ax1.get_yticks()]
    #ax1.set_yticks(Y_ticks)
    #ax1.set_yticklabels(f'{j:g}' for j in Y_ticks) 
    
    mean_genomic_window_size = info['statistics']['window_size_mean']/chr_length(info['chr_id']) 
    ymin,ymax = ax1.get_ylim()
    ax1.errorbar( 0.9-mean_genomic_window_size, ymin + 0.08*(ymax-ymin),marker=None, ls='none', xerr=25*mean_genomic_window_size, linewidth=2*scale, color='k', capsize=4*scale)
    ax1.text(   0.9-mean_genomic_window_size, ymin + 0.05*(ymax-ymin), '25 GW',  horizontalalignment='center', verticalalignment='top',fontsize=2*fs//3)
    ax1.plot([0,1],[0,0],color='black', ls='dotted',alpha=0.5)
    ax1.set_ylim((ymin,ymax))
    ax1.set_xlim((0,1))    

    
    
    #handles, labels = ax1.get_legend_handles_labels()
    if len(H)>1:
        ax1.legend(handles=[i[0] for i in H.values()], title='', loc='upper right', ncol=len(H), fancybox=True,fontsize=int(0.8*fs))
    
    if save!='':
        print('Saving plot...')
        #ax1.set_title(save.rpartition('/')[-1].removesuffix('.png'))
        extension = 'svg'
        plt.tight_layout()
        plt.savefig('.'.join([save,extension]), format=extension, bbox_inches='tight')
        plt.close(fig)
        
    else:
       plt.tight_layout()
       plt.show()
       
def wrap_panel_plot_for_single_indv(identifier, **kwargs):
    """ Wraps the function panel_plot to show all the chromosomes from a single individual. """
    
    DATA = {}
    for i in [*range(1,23)]+['X']:
        llr_filename = kwargs.get('work_dir','.').rstrip('/') + '/' + f'{identifier:s}.chr{str(i):s}.LLR.p.bz2'
        likelihoods,info = load_likelihoods(llr_filename) 
        DATA[f"{info['chr_id'].replace('chr', 'Chromosome '):s}, {info['depth']:.2f}x"] = (likelihoods, info)
        show_info(llr_filename, info, kwargs.get('pairs', (('BPH','SPH'),)))
        kwargs['title'] = identifier
    panel_plot(DATA,**kwargs)
    
    return 0

def wrap_panel_plot_many_cases(filenames, **kwargs):
    """ Wraps the function panel_plot to show a panel with many cases. """
    
    DATA = {}
    for llr_filename in filenames:
        likelihoods,info = load_likelihoods(llr_filename)
        if llr_filename[-6:]=='.LLR.p':
            identifer = llr_filename[:-6].rsplit('/',1).pop()
        elif llr_filename[-10:]=='.LLR.p.bz2':
            identifer = llr_filename[:-10].rsplit('/',1).pop()
        elif llr_filename[-9:]=='.LLR.p.gz':
            identifer = llr_filename[:9].rsplit('/',1).pop()
        else:
            identifer = llr_filename.rsplit('/',1).pop()
        DATA[identifer]=(likelihoods, info)
        show_info(llr_filename, info, kwargs.get('pairs', (('BPH','SPH'),)))
    panel_plot(DATA, **kwargs)
    return 0

def wrap_single_plot(llr_filename, **kwargs):
    """ Wraps the function single_plot. """
    likelihoods,info = load_likelihoods(llr_filename)
    show_info(llr_filename, info, kwargs.get('pairs', (('BPH','SPH'),)))
    single_plot(likelihoods, info, **kwargs)
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plots log-likelihood ratios (LLR) vs. chromosomal position from a LLR file.')
    parser.add_argument('llr_filename', metavar='LLR_FILENAME', type=str, nargs='+',
                        help='One or more LLR files created by ANEUPLOIDY_TEST, containing likelihoods to observese reads under various aneuploidy landscapes .')
    parser.add_argument('-p', '--pairs', type=str, nargs='+', metavar='scenario_A,scenario_B', default=['BPH,SPH'],
                        help='Plots the LLR between scenario A and scenario B along the chromosome. The possible pairs are: BPH,disomy; disomy,SPH; SPH,monosomy; disomy,monosomy; BPH,SPH.'
                             'In addition, giving a list of pairs would plot the LLR of each pair in the same figure, e.g. \"BPH,SPH SPH,MONOSOMY\". The default value is BPH,SPH.')
    parser.add_argument('-b', '--bin-size', type=int, metavar='INT', default=4000000,
                        help='The bin size in which the chromosome is divided. The default value is 4,000,000 bp.')
    parser.add_argument('-z', '--z-score', type=int, metavar='INT', default=1.96,
                        help='The z-score value for the confidence intervals. The default value is 1.96, which corresponds to confidence level of 95%%.')

    kwargs = vars(parser.parse_args())
    kwargs['pairs'] = [j.split(',') for j in kwargs.get('pairs','')]
    
    if  len(kwargs['llr_filename'])==1:
        kwargs['llr_filename'] = kwargs['llr_filename'].pop()
        wrap_single_plot(**kwargs)
    else:
        kwargs['filenames'] = kwargs['llr_filename']
        del kwargs['llr_filename']
        wrap_panel_plot_many_cases(**kwargs)
    sys.exit(0)

else:
    print('The module PLOT_PANEL was imported.')

######## END OF FILE ########    
