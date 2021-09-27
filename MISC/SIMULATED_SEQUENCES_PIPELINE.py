#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

BUILD_SIMULATED_SEQUENCES

Daniel Ariad (daniel@ariad.org)
Dec 30, 2020

"""
import time, sys, random, os, operator, collections

from random import sample, choices, seed
from multiprocessing import Process

sys.path.append('../')
from MIX_HAPLOIDS import MixHaploids_wrapper
from EXTRACT_GENOTYPES import extract

leg_tuple = collections.namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
sam_tuple = collections.namedtuple('sam_tuple', ('sample_id', 'group1', 'group2', 'sex')) #Encodes the rows of the samples table
obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table


def read_ref(filename):
    with open(filename, 'r') as data_in:
        tab = tuple(str(line.replace('\n','')) for line in data_in)
    return tab

def runInParallel(*fns,**kwargs):
    proc = []
    for fn in fns:
        try:
            p = Process(target=fn,args=kwargs.get('args',tuple()))
            p.start()
            proc.append(p)
            time.sleep(5)
        except Exception as error:
            print('caution: a process failed!')
            print(error)
    for p in proc:
        try:
          p.join()
        except:
            None

def transitions(chr_id):
    """ Generates transitions between BPH to SPH regions for trisomy of meiosis II origin. """
    x = int(chr_id[3:]) if chr_id[3:].isnumeric() else chr_id[3:]
    if type(x) is int and 1<=x<=6:
        #BPH-SPH-BPH-SPH
        result = ('BPH',random.uniform(0,.25),random.uniform(.5,.75),random.uniform(.75,1))
    elif type(x) is int and 7<=x<=12:
        #BPH-SPH-BPH
        result = ('BPH',random.uniform(0,.333),random.uniform(.666,1))
    elif x=='X' or (type(x) is int and 13<=x<=22):
        #SPH-BPH
        result = ('SPH',random.uniform(.5,1))
    else:
        result = ('SPH',1)
    return (result,)

def aneuploidy_test_demo(obs_filename,chr_id,sp,model,min_reads,max_reads,output_dir,complex_admixture):
    from ANEUPLOIDY_TEST import aneuploidy_test
    args = dict(obs_filename = f'results_{sp:s}/ABC.obs.p',
                hap_filename = f'../../reference_panels/{sp:s}_panel.hg38/{chr_id:s}_{sp:s}_panel.hap.gz',
                leg_filename = f'../../reference_panels/{sp:s}_panel.hg38/{chr_id:s}_{sp:s}_panel.legend.gz',
                samp_filename = f'../../reference_panels/{sp:s}_panel.hg38/{sp:s}_panel.samples.gz',
                window_size = 0,
                subsamples = 100,
                offset = 0,
                min_reads = min_reads, #3,
                max_reads = max_reads, #8,
                min_HF = 0.05,
                minimal_score = 2,
                output_dir = output_dir, #f'results_{sp:s}/',
                output_filename = '',
                compress = 'bz2')
                #model = model)
    args['obs_filename'] = obs_filename #f'results_{sp:s}/' + obs_filename
    if complex_admixture:
        args['ancestral_proportion'] = (sp[:3],0.5)
    LLR_dict, info = aneuploidy_test(**args)
    return LLR_dict, info

def simulate_haploids(sample_id,sp,chr_id,genotypes,output_dir):
    """ Wraps the function 'extract'. """
    path = '../../reference_panels/'
    #path = f'../build_reference_panel/ref_panel.{sp:s}.hg38/'
    leg_filename = path + f'{sp:s}_panel.hg38/{chr_id:s}_{sp:s}_panel.legend.gz'
    hap_filename = path + f'{sp:s}_panel.hg38/{chr_id:s}_{sp:s}_panel.hap.gz'
    sam_filename = path + f'{sp:s}_panel.hg38/{sp:s}_panel.samples.gz'
    return extract(leg_filename,hap_filename,sam_filename,chr_id,sample_id,genotypes=genotypes,output_dir=output_dir)

def main(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,complex_admixture=False):
    work_dir = work_dir.rstrip('/') + '/' if len(work_dir)!=0 else ''
    #####################
    SPsorted = {('EUR_EAS'): 'EAS_EUR',
                ('EAS_EUR'): 'EAS_EUR',
                ('EUR_SAS'): 'SAS_EUR',
                ('SAS_EUR'): 'SAS_EUR',
                ('EAS_SAS'): 'EAS_SAS',
                ('SAS_EAS'): 'EAS_SAS',
                ('AFR_EUR'): 'AFR_EUR',
                ('EUR_AFR'): 'AFR_EUR',
                'EUR': 'EUR',
                'EAS': 'EAS',
                'SAS': 'SAS',
                'AMR': 'AMR',
                'AFR': 'AFR'}
    #####################
    seed(None, version=2)
    list_SP = sp.split('_')

    if len(list_SP)==1:
        INDIVIDUALS = read_ref(f"../../reference_panels/samples_per_panel/{sp:s}_panel.txt") #EAS_panel.txt')
        A = sample(INDIVIDUALS,k=3)
        B = choices(['A','B'],k=3)
    elif len(list_SP)==2 and not complex_admixture:
        B = choices(['A','B'],k=3)
        A = []
        for i,p in enumerate(random.sample(list_SP, len(list_SP)),start=1):
            INDIVIDUALS = read_ref(f"../../reference_panels/samples_per_panel/{p:s}_panel.txt") #EAS_panel.txt')
            A.extend(sample(INDIVIDUALS,k=i))
    elif len(list_SP)==2 and complex_admixture:
        B = choices(['A','B'],k=6)
        A = []
        for p in list_SP:
            INDIVIDUALS = read_ref(f"../../reference_panels/samples_per_panel/{p:s}_panel.txt") #EAS_panel.txt')
            A.extend(sample(INDIVIDUALS,k=3))
        A = operator.itemgetter(0,3,1,4,2,5)(A)
    else:
        print('error: unsupported sp value.')

    C = [i+j for i,j in zip(A,B)]
    print(C)
    #####################

    for a,b in zip(A,B): simulate_haploids(a, SPsorted[sp], chr_id, b, work_dir)
    sim_obs_tabs = [f'{work_dir:s}{c:s}.{chr_id:s}.hg38.obs.p' for c in C]
    filenames = MixHaploids_wrapper(*sim_obs_tabs, read_length=read_length, depth=depth, scenarios=('disomy','SPH'),
                                    output_dir=work_dir, complex_admixture=complex_admixture)
    print(filenames)

    for f in filenames: aneuploidy_test_demo(f,chr_id,SPsorted[sp],'MODELS/MODELS16.p',
                                             min_reads,max_reads,work_dir, complex_admixture)

    for c in C: os.remove(f'{work_dir:s}{c:s}.{chr_id:s}.hg38.obs.p')
    return 0


if __name__ == "__main__":
    random.seed(a=None, version=2)
    complex_admixture=False
    depth=0.1
    ###sp='EUR_AFR'
    #chr_id='chr21'
    read_length = 36
    min_reads,max_reads = 24,12

    for n in ([*range(16,0,-1),'X']):
        chr_id = 'chr' + str(n)
    #    runInParallel(*([main]*12),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir) )
    #main(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,complex_admixture)
        for sp in ('AFR_EUR','EAS_SAS','SAS_EUR','EAS_EUR'):
            work_dir = f"/mybox/F1-simulations/results_{sp:s}" #'../results' #'results_EAS'
            runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,complex_admixture) )
    print('DONE.')
    pass
else:
    print("The module BUILD_SIMULATED_SEQUENCES was imported.")
