#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 15:21:55 2018

@author: mishugeb
"""
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import nimfa 
import numpy as np
import pandas as pd
from sklearn import metrics
import time
import multiprocessing
from multiprocessing import current_process
from functools import partial
from numpy import linalg as LA
import sigProfilerPlotting as plot
from SigProfilerExtractor import SigProfilerPlotting_plot_png as plotE
import string 
import os
import scipy
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1"
import SigProfilerExtractor as cosmic
from scipy.stats import  ranksums
from SigProfilerExtractor import single_sample as ss
#from sklearn.cluster import KMeans
from SigProfilerExtractor import nmf_cpu
from sklearn.decomposition import NMF
from sklearn import mixture
from scipy.spatial.distance import cdist
from scipy.spatial.distance import correlation as cor
#import mkl
#mkl.set_num_threads(40)
#from numba import jit
#from numba import vectorize
import warnings as _warnings
_warnings.filterwarnings("ignore")

try:
    import torch
    from . import nmf_gpu
except ImportError:
    _warnings.warn('Cannot pytorch - GPU unavailable')

multiprocessing.set_start_method('spawn', force=True)
"""################################################################## Vivid Functions #############################"""

############################################################## FUNCTION ONE ##########################################

def format_integer(number, thousand_separator=','):
    def reverse(string):
        string = "".join(reversed(string))
        return string

    s = reverse(str(number))
    count = 0
    result = ''
    for char in s:
        count = count + 1
        if count % 3 == 0:
            if len(s) == count:
                result = char + result
            else:
                result = thousand_separator + char + result
        else:
            result = char + result
    return result

def make_letter_ids(idlenth = 10, mtype = None):
    
    listOfSignatures = []
    letters = list(string.ascii_uppercase)
    letters.extend([i+b for i in letters for b in letters])
    letters = letters[0:idlenth]
    
    for j,l in zip(range(idlenth),letters)  :
        listOfSignatures.append(mtype+l)
    listOfSignatures = np.array(listOfSignatures)
    return listOfSignatures

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))


def get_indeces(a, b):
    
    """ 
    Extracts the indices multiple items in a list.
    
    Parameters:
        a: list. where we want to get the index of the items.
        b; list. the items we want to get index of. 
    """

    indeces = []
    for i in b:
        try: 
            idx = a.index(i)
            indeces.append(idx)
        except: 
            next

    return indeces 

    """
    #example: 
    x = ['SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS8', 'SBS13', 'SBS40']
    y = ['SBS1',  'SBS5']
    get_indeces(x, y)
    #result
    >>> [1,3]
    """
    

def get_items_from_index(x,y):
    """ decipher the values of items in a list from their indices.
    """
    z = []
    for i in y:
        try:
            z.append(x[i])
        except:
            pass
    return z

    """
    #example: 
    x = [1,3,5,8]
    y = [1, 3]
    get_items_from_index(x, y)
    #result
    >>> [3,8]
    """


############################################################## FUNCTION ONE ##########################################    
def signature_plotting_text(value, text, Type):
    name = text + ": "
    name_list =[]
    total=np.sum(np.array(value))
    for i in value:
        
        if Type=="integer":  
            i = int(i)
            p=round(i/total*100,1)
            i = format(i, ',d')   
            tail = str(i)+'/'+str(p)+'%'
            name_list.append(name+tail)
        elif Type=="float":
            i = round(i,2)
            name_list.append(name + str(i))
    return(name_list)


############################################################## FUNCTION ONE ##########################################
def split_list(lst, splitlenth):
    newlst = []
    for i in range(0, len(lst), splitlenth):
        try: 
            newlst.append([lst[i], lst[i+1]])
        
        except: 
            newlst.append([lst[i]])
            
    return newlst


################################################################### FUNCTION ONE ###################################################################
# Calculates elementwise standard deviations from a list of matrices 
def mat_ave_std(matlst):

    matsum = np.zeros(matlst[0].shape)
    for i in matlst:
        matsum = matsum+i
    matavg = matsum/len(matlst)
    
    matvar=np.zeros(matlst[0].shape)
    for i in matlst:
        matvar = matvar+(i-matavg)**2
    matvar= matvar/len(matlst)
    matstd = matvar**(1/2)
    
    return matavg, matstd


"""############################################## Functions for modifications of Sample Matrices ##################"""   

def get_normalization_cutoff(data, manual_cutoff=50000):
    
    col_sums = np.array(np.sum(data, axis=0))

    # continue the loop if the differece the means is larger than the 2*2*STD of the larger cluster 
    while True:

        try:
            # doing Kmean clustering
            col_sums_for_cluster = col_sums.reshape(-1,1)
            
            #separate distributions using kmean
            #kmeans = KMeans(n_clusters=2, random_state=0).fit(col_sums_for_cluster)
            #labels = kmeans.labels_
            
            #separate distributions using mixture model
            clf = mixture.GaussianMixture(n_components=2, covariance_type='full')
            clf.fit(col_sums_for_cluster)
            labels = clf.predict(col_sums_for_cluster)
            
            unique, counts = np.unique(labels, return_counts=True)
            if len(unique)==1:
                break
            bigger_cluster = unique[np.argmax(counts)]
            smaller_cluster = unique[np.argmin(counts)]
            
            # estimating the magnitute of discripancy better the clusters
            bigger_cluster__dist = col_sums[labels==bigger_cluster]
            smaller_cluster__dist = col_sums[labels==smaller_cluster]
            bigger_cluster__dist_mean = np.mean(bigger_cluster__dist) 
            bigger_cluster__dist_std = np.std(bigger_cluster__dist) 
            smaller_cluster__dist_mean = np.mean(smaller_cluster__dist) 
            #smaller_cluster__dist_std = np.std(smaller_cluster__dist) 
            
             
            
            #print("bigger_cluster__dist_mean ", bigger_cluster__dist_mean)
            #print("bigger_cluster__dist_std ", bigger_cluster__dist_std)
            #print("smaller_cluster__dist_mean ", smaller_cluster__dist_mean)
            #print("smaller_cluster__dist_std ", smaller_cluster__dist_std)
            
            # continue the loop if the differece the means is larger than the 2*STD of the larger cluster 
            
            if abs(bigger_cluster__dist_mean-smaller_cluster__dist_mean)< 2*2*bigger_cluster__dist_std:
                break
                
            
            # col_sums will be equal to bigger_cluster__dist for the next iteration 
            col_sums = bigger_cluster__dist
        except:
            break
            
     
    mean = np.mean(col_sums)  
    std = np.std(col_sums)
    cutoff = (mean + 2*(std)).astype(int)
    
    if cutoff>manual_cutoff:
        cutoff = manual_cutoff
    
    return cutoff
        
def normalize_samples(bootstrapGenomes,totalMutations,norm="100X", normalization_cutoff=10000000):
        if norm=="gmm":
            bootstrapGenomes = np.array(bootstrapGenomes)
            indices = np.where(totalMutations>normalization_cutoff)[0]
            norm_genome = bootstrapGenomes[:,list(indices)]/totalMutations[list(indices)][:,np.newaxis].T*normalization_cutoff
            bootstrapGenomes[:,list(indices)] = norm_genome
            bootstrapGenomes = pd.DataFrame(bootstrapGenomes)
        elif norm == "100X":
            bootstrapGenomes = np.array(bootstrapGenomes)
            rows = bootstrapGenomes.shape[0]
            indices = np.where(totalMutations>(rows*100))[0]
            norm_genome = bootstrapGenomes[:,list(indices)]/totalMutations[list(indices)][:,np.newaxis].T*(rows*100)
            bootstrapGenomes = pd.DataFrame(bootstrapGenomes)
            
        elif norm == "log2":
            log2_of_tM = np.log2(totalMutations)
            bootstrapGenomes = bootstrapGenomes/totalMutations*log2_of_tM
        elif norm == "none":
            pass
        else:
            try:
                bootstrapGenomes = np.array(bootstrapGenomes)
                rows = bootstrapGenomes.shape[0]
                indices = np.where(totalMutations>int(norm))[0]
                norm_genome = bootstrapGenomes[:,list(indices)]/totalMutations[list(indices)][:,np.newaxis].T*(int(norm))
                bootstrapGenomes = pd.DataFrame(bootstrapGenomes)
            except:
                pass
        return bootstrapGenomes


    

def split_samples(samples, intervals, rescaled_items, colnames):

    colnames = np.array(colnames)
    total_mutations = samples.sum(axis =0)

    
    max_total = np.max(total_mutations)+1
    
    intervals = np.array(intervals )
    rescaled_items = np.array(rescaled_items)
    selected_indices, = np.where(intervals<=max_total)
    
    intervals=intervals[selected_indices]
    rescaled_items=rescaled_items[selected_indices]
   
    rescaled_items = list(rescaled_items)
   
    
    sample_list = []
    sample_total = []
    rescale_list = []
    rescale_values = []
    colnames_list = []
    ranges = list(intervals)+[max_total]
    #print(ranges)
    for i in range(len(ranges)-1):    
        sub_sample, = np.where((total_mutations<=ranges[i+1]) & (total_mutations>ranges[i]))
        if samples[:,sub_sample].shape[1]>0:
            sample_list.append(samples[:,sub_sample])
            colnames_list.append(list(colnames[sub_sample]))
            sample_total.append(np.sum(samples[:,sub_sample], axis=0))
            rescale_list.append(rescaled_items[i])
            rescale_values.append(ranges[i])
      
    return(sample_list, colnames_list, sample_total, rescale_list, rescale_values) 
    
   
def denormalize_samples(genomes, original_totals, normalization_value=30000):
    normalized_totals = np.sum(genomes, axis=0)
    original_totals = np.array(original_totals)
    results = genomes/normalized_totals*original_totals
    results = np.round(results,0)
    results = results.astype(int)
    return results 
##################################################################################################################    

"""###################################################### Fuctions for NON NEGATIVE MATRIX FACTORIZATION (NMF) #################"""

def matmul(a,b):
    return a*b


def matdiv(a,b):
    return a/b
    

    
    
    
def inhouse_nmf(v, w=0, h=0, k=2, iterations=200000,tol=None):
    
    v = np.float32(v)
    #min_v = np.min(v[v>0])
    n,m = v.shape[0], v.shape[1]
    EPS = np.float32(np.finfo(float).eps)
    #h[h <= min_v] = min_v
    #w[w <= min_v] = min_v
    initcost = 1000000000000 # a big number
    conv = 1
    for i in range(iterations):
        
        #updata rule for h
        x1 = np.repeat((np.sum(w,axis=0).T)[:,np.newaxis], m, axis=1)
        #print(np.linalg.norm(x1))
        dot1 = np.dot(w.T,(v/(np.dot(w,h))))
        h= (h*dot1)/x1
        #print(np.linalg.norm(h))
        
        #updata rule for w
        x2 = np.repeat((np.sum(h,axis=1).T)[np.newaxis,:], n, axis=0)
        #print(np.linalg.norm(x2))
        dot2 =np.dot((v/(np.dot(w,h))),h.T)
        w=(w*dot2)/x2
        #print(np.linalg.norm(w))
        
        #Adjust small values every ten steps to avoid undeflow
        if (i+1)%10==0:
            h[h <= EPS] = EPS
            w[w <= EPS] = EPS
        if (i+1)%1000==0:    
            est_v = np.dot(w, h)
            
            norm_v = np.linalg.norm(v)
            norm_est_v = np.linalg.norm(est_v)
            
            cost = norm_est_v-norm_v
            
            diff = abs(initcost-cost)
            initcost = cost
            
            
            if diff <= tol:
               #print(diff)
               #print("diff is smaller than tol")
               conv+=1
               if conv==3:
                  break
            elif diff>tol:
               conv =1
    
    
    return w, h
    

def nnmf_cpu(genomes, nfactors, init="nndsvd", excecution_parameters=None):
   
    
    genomes = torch.from_numpy(genomes).float()
    min_iterations=excecution_parameters["min_NMF_iterations"]
    max_iterations=excecution_parameters["max_NMF_iterations"]
    tolerance=excecution_parameters["NMF_tolerance"]
    test_conv=excecution_parameters["NMF_test_conv"]
    net = nmf_cpu.NMF(genomes,rank=nfactors, min_iterations=min_iterations, max_iterations=max_iterations, tolerance=tolerance,test_conv=test_conv, init_method=init,seed=None)
    net.fit()
    Ws = []
    Hs = []
    
    for H in net.H.detach().cpu().numpy():
        Hs.append(np.matrix(H))
    for W in net.W.detach().cpu().numpy():
        Ws.append(np.matrix(W))
    
    convergence = int(net.conv)
    
    W = Ws[0]
    H = Hs[0]
    #calculate L1, L2 and KL for the solution 
    
    est_genome = np.array(np.dot(W, H))
    genomes = np.array(genomes)
    similarities = calculate_similarities(genomes, est_genome, sample_names=False)[0].iloc[:,2:]
    similarities = np.array(np.mean(similarities, axis=0)).T
    similarities = np.append(similarities, convergence)
    return W, H, similarities

def nnmf_gpu(genomes, nfactors, init="nndsvd",excecution_parameters=None):
    p = current_process()
    identity = p._identity[0]
    #print(genomes.shape)
    gpu_id = identity % torch.cuda.device_count()
    genomes = torch.from_numpy(genomes).float().cuda(gpu_id)
    min_iterations=excecution_parameters["min_NMF_iterations"]
    max_iterations=excecution_parameters["max_NMF_iterations"]
    tolerance=excecution_parameters["NMF_tolerance"]
    test_conv=excecution_parameters["NMF_test_conv"]
    net = nmf_gpu.NMF(genomes,rank=nfactors,min_iterations=min_iterations,max_iterations=max_iterations, tolerance=tolerance,test_conv=test_conv, gpu_id=gpu_id, init_method=init,seed=None)
    net.fit()
    Ws = []
    Hs = []
    for H in net.H.detach().cpu().numpy():
        Hs.append(np.matrix(H))
    for W in net.W.detach().cpu().numpy():
        Ws.append(np.matrix(W))
    
    if len(Ws)==1:  
        convergence = int(net.conv)
    else:
        convergence = int(max_iterations)
        
    convergence=[convergence]*len(Ws)

    return Ws, Hs, convergence

def nnmf(genomes, nfactors, init="nndsvd"):
    
   
    genomes = np.array(genomes)
   
    #print(init)       
    #model_init = NMF(n_components=nfactors, max_iter=1, solver= "mu",init= "nndsvd", tol=0.0, beta_loss = 'kullback-leibler',  verbose=False)
    #w = model_init.fit_transform(genomes)
    #h = model_init.components_
    w,h=initialize_nmf(genomes, nfactors, init=init, eps=1e-6,random_state=None)
    W, H = inhouse_nmf(genomes, w=w, h=h, k=nfactors, iterations=200000, tol=0.000005)
   
    
    
    
    #calculate L1, L2 and KL for the solution 
    est_genome = np.array(np.dot(W, H))
    
    similarities = calculate_similarities(genomes, est_genome, sample_names=False)[0].iloc[:,2:]
    similarities = np.array(np.mean(similarities, axis=0)).T
    
    
    

    return W, H, similarities




def BootstrapCancerGenomes(genomes, seed=None):
    np.random.seed(seed) # Every time initiate a random seed so that all processors don't get the same seed
    
    """
    index = genomes.index
    cols = genomes.columns
    genomes = np.array(genomes)
    n = 100
    stderr = genomes*0.025
    std = stderr*math.sqrt(n)
    genomes = np.random.normal(genomes, std, size=None).astype(int)
    genomes[genomes<0] = 0

    dataframe = pd.DataFrame(genomes)
    dataframe.index = index
    dataframe.columns = cols 
    """
    #normalize the data 
    a = pd.DataFrame(genomes.sum(0))  #performing the sum of each column. variable a store the sum
    a = a.transpose()  #transfose the value so that if gets a proper dimention to normalize the data
    repmat = pd.concat([a]*genomes.shape[0])  #normalize the data to get the probility of each datapoint to perform the "random.multinomial" operation
    normGenomes = genomes/np.array(repmat)      #part of normalizing step
    
    
    #Get the normGenomes as Matlab/ alternation of the "mnrnd" fuction in Matlap
    #do the Boostraping 
    dataframe = pd.DataFrame() #declare an empty variable to populate with the Bootstrapped data
    for i in range(0,normGenomes.shape[1]):
        #print a.iloc[:,i]
        #print normGenomes.iloc[:,i]
        dataframe[i]=list(np.random.multinomial(a.iloc[:,i], normGenomes.iloc[:,i], 1)[0])
        
   
    return dataframe

# NMF version for the multiprocessing library
def pnmf(batch_seed_pair=[1,None], genomes=1, totalProcesses=1, resample=True, init="nndsvd", seeds=None, normalization_cutoff=10000000, norm="log2", gpu=False, excecution_parameters=None):
    tic = time.time()
    totalMutations = np.sum(genomes, axis =0)
    genomes = pd.DataFrame(genomes) #creating/loading a dataframe/matrix
    
    if gpu:
        batch_size=batch_seed_pair[0]
        seeds=batch_seed_pair[1]
        nmf_fn = nnmf_gpu
        results = []
        genome_list = []

        for b in range(batch_size):
            if resample == True:
                bootstrapGenomes= BootstrapCancerGenomes(genomes, seed=seeds)
            else: 
                bootstrapGenomes=genomes    
            
            bootstrapGenomes[bootstrapGenomes<0.0001]= 0.0001
            totalMutations = np.sum(bootstrapGenomes, axis=0)
            
              
            bootstrapGenomes=normalize_samples(bootstrapGenomes,totalMutations,norm=norm, normalization_cutoff=normalization_cutoff)
                    
            
            genome_list.append(bootstrapGenomes.values)
            
            #print(genomes.shape)
        #print(len(genome_list))      
        g = np.array(genome_list)
        W, H, Conv = nmf_fn(g, totalProcesses, init=init, excecution_parameters=excecution_parameters)
        for i in range(len(W)):
            
            _W = np.array(W[i])
            _H = np.array(H[i])
            total = _W.sum(axis=0)[np.newaxis]
            _W = _W/total
            _H = _H*total.T
            _H=denormalize_samples(_H, totalMutations) 
            _conv=Conv[i]
            results.append((_W, _H, _conv))
            print ("process " +str(totalProcesses)+" continues please wait... ")
            print ("execution time: {} seconds \n".format(round(time.time()-tic), 2))

        return results

    else:
        nmf_fn = nnmf_cpu
        seeds=batch_seed_pair[1]
       
        if resample == True:
            bootstrapGenomes= BootstrapCancerGenomes(genomes, seed=seeds)
        else:
            bootstrapGenomes=genomes
            
            #print(pool_constant)
            #print(bootstrapGenomes.iloc[0,:].T)
            #print("\n\n\n")
            
        
        bootstrapGenomes[bootstrapGenomes<0.0001]= 0.0001
        # normalize the samples to handle the hypermutators
       
        totalMutations = np.sum(bootstrapGenomes, axis=0)
        
        #print(normalization_cutoff)
        bootstrapGenomes=normalize_samples(bootstrapGenomes,totalMutations,norm=norm, normalization_cutoff=normalization_cutoff)
            
        bootstrapGenomes=np.array(bootstrapGenomes)
    
        W, H, kl = nmf_fn(bootstrapGenomes,totalProcesses, init=init, excecution_parameters=excecution_parameters)  #uses custom function nnmf
        
        
        #print ("initital W: ", W); print("\n");
        #print ("initial H: ", H); print("\n");
        W = np.array(W)
        H = np.array(H)
        total = W.sum(axis=0)[np.newaxis]
        #print ("total: ", total); print("\n");
        W = W/total
        H = H*total.T
        
        # denormalize H
        H = denormalize_samples(H, totalMutations) 
        print ("process " +str(totalProcesses)+" continues please wait... ")
        print ("execution time: {} seconds \n".format(round(time.time()-tic), 2))
        
        
        return W, H, kl



# =============================================================================
# def pnmf(seed=None, genomes=1, totalProcesses=1, resample=True, init="random", normalization_cutoff=10000000, gpu=False):
#     
#     tic = time.time()
#     
#     
#     totalMutations = np.sum(genomes, axis =0)
#     genomes = pd.DataFrame(genomes) #creating/loading a dataframe/matrix
#     
#     
#     if gpu:
#         nmf_fn = nnmf_gpu
#     else:
#         nmf_fn = nnmf
# 
#     if resample == True:
#         bootstrapGenomes= BootstrapCancerGenomes(genomes, seed=seed)
#         
#         #print(pool_constant)
#         #print(bootstrapGenomes.iloc[0,:].T)
#         #print("\n\n\n")
#         
#         
#         bootstrapGenomes[bootstrapGenomes<0.0001]= 0.0001
#         
#         # normalize the samples to handle the hypermutators
#         bootstrapGenomes = np.array(bootstrapGenomes)
#         #print(normalization_cutoff)
#         #bootstrapGenomes = normalize_samples(bootstrapGenomes, normalize=True, all_samples=False, number=normalization_cutoff)
#         #print(type(bootstrapGenomes))
#         
#         totalMutations = np.sum(bootstrapGenomes, axis=0)
#         log2_of_tM = np.log2(totalMutations)
#         bootstrapGenomes = bootstrapGenomes/totalMutations*log2_of_tM
#         W, H, kl = nmf_fn(bootstrapGenomes,totalProcesses, init=init)  #uses custom function nnmf
#         
#     else:
#         genomes = normalize_samples(genomes, normalize=False, all_samples=False, number=normalization_cutoff)
#         W, H, kl = nmf_fn(genomes,totalProcesses, init= init)  #uses custom function nnmf
#     #print ("initital W: ", W); print("\n");
#     #print ("initial H: ", H); print("\n");
#     W = np.array(W)
#     H = np.array(H)
#     total = W.sum(axis=0)[np.newaxis]
#     #print ("total: ", total); print("\n");
#     W = W/total
#     H = H*total.T
#     
#     # denormalize H
#     H = denormalize_samples(H, totalMutations, normalization_value=totalMutations) 
#     print ("process " +str(totalProcesses)+" continues please wait... ")
#     print ("execution time: {} seconds \n".format(round(time.time()-tic), 2))
#     
#     
#     return W, H, kl
# =============================================================================

def parallel_runs(excecution_parameters, genomes=1, totalProcesses=1, verbose = False):
    
    iterations = excecution_parameters["NMF_replicates"]
    seeds=excecution_parameters["seeds"]
    init=excecution_parameters["NMF_init"]
    normalization_cutoff=excecution_parameters["normalization_cutoff"]
    n_cpu=excecution_parameters["cpu"]
    resample=excecution_parameters["resample"]
    norm=excecution_parameters["matrix_normalization"]
    gpu=excecution_parameters["gpu"]
    batch_size=excecution_parameters["batch_size"]
    
    
    
    
    #random_seeds=True
    # set up the seeds generation same matrices for different number of signatures
    #if random_seeds==True:
        #seeds = np.random.randint(0, 10000000, size=iterations) # set the seeds ranging from 0 to 10000000 for resampling and same seeds are used in different number of signatures
    #else:
        #seeds = list(range(0,iterations))
    
    if verbose:
        print ("Process "+str(totalProcesses)+ " is in progress\n===================================>")
    if n_cpu==-1:
        pool = multiprocessing.Pool()
    else:
        pool = multiprocessing.Pool(processes=n_cpu)

    num_full_batches = iterations // batch_size
    last_batch_size = iterations % batch_size

    batches = [batch_size for _ in range(num_full_batches)]
    if last_batch_size != 0:
        batches.append(last_batch_size)
    
    batch_seed_pair = []
    for i,j in zip(batches,seeds):
        batch_seed_pair.append([i,j])
   
   
    if gpu==True:
        pool_nmf=partial(pnmf, genomes=genomes, totalProcesses=totalProcesses, resample=resample, seeds=seeds, init=init, normalization_cutoff=normalization_cutoff, norm=norm, gpu=gpu, excecution_parameters=excecution_parameters)
        result_list = pool.map(pool_nmf, batch_seed_pair) 
        pool.close()
        pool.join()
        #result_list_flattened = [s for sublist for sublist in result_list]
        flat_list = [item for sublist in result_list for item in sublist]

    else:
         pool_nmf=partial(pnmf, genomes=genomes, totalProcesses=totalProcesses, resample=resample, init=init, seeds=seeds,normalization_cutoff=normalization_cutoff, norm=norm, gpu=gpu, excecution_parameters=excecution_parameters)
         result_list = pool.map(pool_nmf, batch_seed_pair) 
         pool.close()
         pool.join()
         flat_list = result_list
    return flat_list#result_list
# =============================================================================
# def parallel_runs(genomes=1, totalProcesses=1, iterations=1,  n_cpu=-1, verbose = False, resample=True, seeds = None, init="random", normalization_cutoff=10000000, gpu=False):
#     if verbose:
#         print ("Process "+str(totalProcesses)+ " is in progress\n===================================>")
#     if n_cpu==-1:
#         pool = multiprocessing.Pool()
#     else:
#         pool = multiprocessing.Pool(processes=n_cpu)
#         
#     #print(seeds)
#     pool_nmf=partial(pnmf, genomes=genomes, totalProcesses=totalProcesses, resample=resample, init=init, normalization_cutoff=normalization_cutoff, gpu=gpu)
#     result_list = pool.map(pool_nmf, seeds) 
#     pool.close()
#     pool.join()
#     
#     return result_list
# =============================================================================
####################################################################################################################



"""
#############################################################################################################
#################################### Decipher Signatures ###################################################
#############################################################################################################
"""
def decipher_signatures(excecution_parameters, genomes=[0], i=1, totalIterations=1, cpu=-1, mut_context="96"):
    
    
        
    m = mut_context
    
    tic = time.time()
    # The initial values accumute the results for each number of 
    totalMutationTypes = genomes.shape[0];
    totalGenomes = genomes.shape[1];
    totalProcesses = i
    totalIterations=excecution_parameters["NMF_replicates"]
    gpu=excecution_parameters["gpu"]
    dist=excecution_parameters["dist"]
    norm = excecution_parameters["matrix_normalization"]
    normalization_cutoff=excecution_parameters["normalization_cutoff"]
    
    print ("Extracting signature {} for mutation type {}".format(i, m))  # m is for the mutation context
    
    if norm=="gmm":
        print("The matrix normalizig cutoff is {}\n\n".format(normalization_cutoff))
    else:
        print("The matrix normalizig cutoff is set for {}\n\n".format(norm))
    
    
    
    ##############################################################################################################################################################################         
    ############################################################# The parallel processing takes place here #######################################################################  
    ############################################################################################################################################################################## 
    if gpu==True:
        results = []
        flat_list = parallel_runs(excecution_parameters, genomes=genomes, totalProcesses=totalProcesses,  verbose = False)
        
        for items in range(len(flat_list)):
            
            W = flat_list[items][0]
            H = flat_list[items][1]
            conv=flat_list[items][2]
            #calculate L1, L2 and KL for the solution 
            est_genome = np.array(np.dot(W, H))
            
            similarities = calculate_similarities(genomes, est_genome, sample_names=False)[0].iloc[:,2:]
            similarities = np.array(np.mean(similarities, axis=0)).T
            similarities = np.append(similarities, conv)
            
            results.append([W,H,similarities])
    else:
        results = parallel_runs(excecution_parameters, genomes=genomes, totalProcesses=totalProcesses, verbose = False)
    #print(results[0][2])     
    toc = time.time()
    print ("Time taken to collect {} iterations for {} signatures is {} seconds".format(totalIterations , i, round(toc-tic, 2)))
    ##############################################################################################################################################################################       
    ######################################################### The parallel processing ends here ##################################################################################      
    ##############################################################################################################################################################################        
    
    
    ################### Achieve the best clustering by shuffling results list using a few iterations ##########        
    Wall = np.zeros((totalMutationTypes, totalProcesses * totalIterations));
    #print (Wall.shape)
    Hall = np.zeros((totalProcesses * totalIterations, totalGenomes));
    converge_information = np.zeros((totalIterations, 6))
    
    finalgenomeErrors = np.zeros((totalMutationTypes, totalGenomes, totalIterations));
    finalgenomesReconstructed = np.zeros((totalMutationTypes, totalGenomes, totalIterations))
    
    processCount=0
    for j in range(len(results)):
        W = results[j][0]
        H = results[j][1]
        converge_information[j,:] = results[j][2][:]
        finalgenomeErrors[:, :, j] = genomes -  np.dot(W,H);
        finalgenomesReconstructed[:, :, j] = np.dot(W,H);
        Wall[ :, processCount : (processCount + totalProcesses) ] = W;
        Hall[ processCount : (processCount + totalProcesses), : ] = H;
        processCount = processCount + totalProcesses;
    
    
    processes=i #renamed the i as "processes"    
    processAvg, exposureAvg, processSTE,  exposureSTE, avgSilhouetteCoefficients, clusterSilhouetteCoefficients = cluster_converge_outerloop(Wall, Hall, processes, dist=dist, gpu=gpu)
    reconstruction_error = round(LA.norm(genomes-np.dot(processAvg, exposureAvg), 'fro')/LA.norm(genomes, 'fro'), 2)   
    

    return  processAvg, exposureAvg, processSTE, exposureSTE, avgSilhouetteCoefficients, np.round(clusterSilhouetteCoefficients,3), finalgenomeErrors, finalgenomesReconstructed, Wall, Hall, converge_information, reconstruction_error, processes



"""################################################################### FUNCTIONS TO CALCULATE DISTANCES BETWEEN VECTORS ###################################################################"""
################################################################### FUNCTION ONE ###################################################################
#function to calculate the cosine similarity
def cos_sim(a, b):
      
    
    """Takes 2 vectors a, b and returns the cosine similarity according 
    to the definition of the dot product
    
    Dependencies: 
    *Requires numpy library. 
    *Does not require any custom function (constructed by me)
    
    Required by:
    * pairwise_cluster_raw
    	"""
    if np.sum(a)==0 or np.sum(b) == 0:
        return 0.0      
    dot_product = np.dot(a, b)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    return dot_product / (norm_a * norm_b)

def cor_sim(a, b):
      
    
    """Takes 2 vectors a, b and returns the corrilation similarity according 
    to the definition of the dot product
    
    Dependencies: 
    *Requires numpy library. 
    *Does not require any custom function (constructed by me)
    
    Required by:
    * pairwise_cluster_raw
    	"""
    if np.sum(a)==0 or np.sum(b) == 0:
        return 0.0      
    corr =1-cor(a, b)
    return corr


    
################################################################### FUNCTION ONE ###################################################################
#function to calculate multiple similarities/distances
def calculate_similarities(genomes, est_genomes, sample_names=False):
    from numpy import inf
    
    if  sample_names is False:
        sample_names = ["None"]*genomes.shape[1]
        
    cosine_similarity_list = []
    kl_divergence_list = []
    l1_norm_list = []
    l2_norm_list = []
    total_mutations_list = []
    relative_l1_list = []
    relative_l2_list = []
    for i in range(genomes.shape[1]):
        p_i = genomes[:,i]
        q_i = est_genomes[:, i]
        #maxqp_i = np.max(np.concatenate((p_i[:,np.newaxis], q_i[:,np.newaxis]), axis=1), axis=1)
        #p_i = p_i/np.sum(p_i)*100
        #q_i = q_i/np.sum(q_i)*100
        
        cosine_similarity_list.append(round(cos_sim(p_i,q_i ),3))
        kl_divergence_list.append(round(scipy.stats.entropy(p_i,q_i),4))
        l1_norm_list.append(round(np.linalg.norm(p_i-q_i , ord=1),3))
        relative_l1_list.append(round((l1_norm_list[-1]/np.linalg.norm(p_i, ord=1))*100,3))
        l2_norm_list.append(round(np.linalg.norm(p_i-q_i , ord=2),3))
        relative_l2_list.append(round((l2_norm_list[-1]/np.linalg.norm(p_i, ord=2))*100,3))
        total_mutations_list.append(np.sum(p_i))
    kl_divergence_list = np.array(kl_divergence_list)
    kl_divergence_list[kl_divergence_list == inf] =1000
    similarities_dataframe = pd.DataFrame({"Sample Names": sample_names, \
                                           "Total Mutations":total_mutations_list, \
                                           "Cosine Similarity": cosine_similarity_list, \
                                           "L1 Norm": l1_norm_list, \
                                           "L1_Norm_%":relative_l1_list, \
                                           "L2 Norm": l2_norm_list, \
                                           "L2_Norm_%": relative_l2_list, \
                                           "KL Divergence": kl_divergence_list})
    similarities_dataframe = similarities_dataframe.set_index("Sample Names")
    return [similarities_dataframe, cosine_similarity_list]




"""################################################################### CLUSTERING FUNCTIONS ###################################################################"""
################################################################### FUNCTION ONE ###################################################################
# function to calculate the centroids

################################################################### FUNCTION  ###################################################################
def pairwise_cluster_raw(mat1=([0]), mat2=([0]), mat1T=([0]), mat2T=([0]), dist="cosine", gpu=False):  # the matrices (mat1 and mat2) are used to calculate the clusters and the lsts will be used to store the members of clusters
    
    """ Takes a pair of matrices mat1 and mat2 as arguments. Both of the matrices should have the 
    equal shapes. The function makes a partition based clustering (the number of clusters is equal 
    to the number of colums of the matrices, and not more column is assigned into a cluster from 
    a single matrix). It return the list of clusters  as "lstCluster" and the list of clusters 
    based on their indeces of the vectors in their original matrices. Please run the 
    test of function/example code provided below the for better understanding. 
    
     Dependencies:
        *cos_sim
        
    Required by:
        *pairwise_cluster_init 
        *pairwise_cluster_elong
    
    
    """

    
   
        
    if dist=="cosine":
        con_mat = 1-cdist(mat1.T, mat2.T, "cosine")
    elif dist=="correlation":
        con_mat = 1-cdist(mat1.T, mat2.T, "correlation")
        
    
    maximums = np.argmax(con_mat, axis = 1) #the indices of maximums
    uniques = np.unique(maximums) #unique indices having the maximums

        
    lstCluster=[]
    lstClusterT=[]
    idxPair= []  
    
    # check if more than one centroids have their maximums in the same indices
    if len(maximums)==1000000000000:#!=len(uniques):
        diffRank = np.zeros([len(con_mat), 2])    
        for i in range(len(con_mat)):
            centroid = con_mat[i]
            value_order= np.argsort(centroid)
            first_highest_value_idx = value_order[-1]
            second_highest_value_idx = value_order[-2]
            diffRank[i,0] =  i
            diffRank[i, 1] = centroid[first_highest_value_idx]-centroid[second_highest_value_idx]
        diffRank = diffRank[diffRank[:,1].argsort(kind='mergesort')]
        diffRank = diffRank.astype(int)

  
    
        for a in range(len(diffRank)-1,-1, -1):
            i = diffRank[a,0]
            j = con_mat[i,:].argmax()
            con_mat[i,:]=0
            con_mat[:,j]=0
            lstCluster.append([mat1[:,i], mat2[:,j]])
            idxPair.append([i,j])   #for debugging
            lstClusterT.append([mat1T[i,:], mat2T[j, :]])
    
    
    else:
        for a in range(0,mat1.shape[1]):
        
            i,j=np.unravel_index(con_mat.argmax(), con_mat.shape) 
            lstCluster.append([mat1[:,i], mat2[:,j]])
            idxPair.append([i,j])   #for debugging
            lstClusterT.append([mat1T[i,:], mat2T[j, :]])
            con_mat[i,:]=0
            con_mat[:,j]=0
        
    return lstCluster, idxPair, lstClusterT
         





################################################################### FUNCTION  ###################################################################
def reclustering(tempWall=0, tempHall=0, processAvg=0, exposureAvg=0, dist="cosine",gpu=False):
    # exposureAvg is not important here. It can be any matrix with the same size of a single exposure matrix
    iterations = int(tempWall.shape[1]/processAvg.shape[1])
    processes =  processAvg.shape[1]
    idxIter = list(range(0, tempWall.shape[1], processes))
   
    processes3D = np.zeros([processAvg.shape[1], processAvg.shape[0], iterations])
    exposure3D = np.zeros([exposureAvg.shape[0], iterations, exposureAvg.shape[1]])
    
    for  iteration_number in range(len(idxIter)):
        
        
        #print(i)
        statidx = idxIter[iteration_number]
        loopidx = list(range(statidx, statidx+processes))
        lstCluster, idxPair, lstClusterT = pairwise_cluster_raw(mat1=processAvg, mat2=tempWall[:, loopidx], mat1T=exposureAvg, mat2T=tempHall[loopidx,:],dist=dist, gpu=gpu)
        
        for cluster_items in idxPair:
            cluster_number = cluster_items[0]
            query_idx = cluster_items[1]
            processes3D[cluster_number,:,iteration_number]=tempWall[:,statidx+query_idx]
            exposure3D[cluster_number, iteration_number, :] = tempHall[statidx+query_idx,:]
            
    
    
    count = 0
    labels=[]
    clusters = pd.DataFrame()
    for cluster_id in range(processes3D.shape[0]):
        cluster_vectors = pd.DataFrame(processes3D[cluster_id,:,:])
        clusters = clusters.append(cluster_vectors.T)
        for k in range(cluster_vectors.shape[1]):
            labels.append(count)
        count= count+1
    
    
    
    try:
        if dist=="cosine":
            SilhouetteCoefficients = metrics.silhouette_samples(clusters, labels, metric='cosine')
        if dist=="correlation":
            SilhouetteCoefficients = metrics.silhouette_samples(clusters, labels, metric='correlation')
            
    
        
    
    except:
        SilhouetteCoefficients = np.ones((len(labels),1))
        
        
        
        
    avgSilhouetteCoefficients = np.mean(SilhouetteCoefficients)
    
    #clusterSilhouetteCoefficients 
    splitByCluster = np.array_split(SilhouetteCoefficients, processes3D.shape[0])
    clusterSilhouetteCoefficients = np.array([])
    for i in splitByCluster:
        
        clusterSilhouetteCoefficients=np.append(clusterSilhouetteCoefficients, np.mean(i))
        
    processAvg = np.mean(processes3D, axis=2).T
    processSTE = scipy.stats.sem(processes3D, axis=2, ddof=1).T
    exposureAvg = np.mean(exposure3D, axis=1) 
    exposureSTE = scipy.stats.sem(exposure3D, axis=1, ddof=1)
    
        
    return  processAvg, exposureAvg, processSTE,  exposureSTE, avgSilhouetteCoefficients, clusterSilhouetteCoefficients




def cluster_converge_innerloop(Wall, Hall, totalprocess, iteration=1, dist="cosine", gpu=False):
    
    processAvg = np.random.rand(Wall.shape[0],totalprocess)
    exposureAvg = np.random.rand(totalprocess, Hall.shape[1])
    
    result = 0
    convergence_count = 0
    while True:
        processAvg, exposureAvg, processSTE,  exposureSTE, avgSilhouetteCoefficients, clusterSilhouetteCoefficients = reclustering(Wall, Hall, processAvg, exposureAvg, dist=dist, gpu=gpu)
        
        if result == avgSilhouetteCoefficients:
            break
        elif convergence_count == 10:
            break
        else:
            result = avgSilhouetteCoefficients
            convergence_count = convergence_count + 1
        
    return processAvg, exposureAvg, processSTE,  exposureSTE, avgSilhouetteCoefficients, clusterSilhouetteCoefficients



def parallel_clustering(Wall, Hall, totalProcesses, iterations=50,  n_cpu=-1, dist= "cosine", gpu=False):
    
    if n_cpu==-1:
        pool = multiprocessing.Pool()
    else:
        pool = multiprocessing.Pool(processes=n_cpu)
        
    pool_nmf=partial(cluster_converge_innerloop, Wall, Hall, totalProcesses, dist=dist, gpu=gpu)
    result_list = pool.map(pool_nmf, range(iterations)) 
    pool.close()
    pool.join()
    return result_list

# To select the best clustering converge of the cluster_converge_innerloop
def cluster_converge_outerloop(Wall, Hall, totalprocess, dist="cosine", gpu=False):
    avgSilhouetteCoefficients = -1  # intial avgSilhouetteCoefficients 
    
    #do the parallel clustering 
    result_list = parallel_clustering(Wall, Hall, totalprocess, iterations=50,  n_cpu=-1,  dist=dist, gpu=gpu)
    
    for i in range(50):  # using 10 iterations to get the best clustering 
        
        temp_processAvg, temp_exposureAvg, temp_processSTE,  temp_exposureSTE, temp_avgSilhouetteCoefficients, temp_clusterSilhouetteCoefficients = result_list[i][0], result_list[i][1], result_list[i][2], result_list[i][3], result_list[i][4], result_list[i][5]
        
        if avgSilhouetteCoefficients < temp_avgSilhouetteCoefficients:
              processAvg, exposureAvg, processSTE,  exposureSTE, avgSilhouetteCoefficients, clusterSilhouetteCoefficients =   temp_processAvg, temp_exposureAvg, temp_processSTE,  temp_exposureSTE, temp_avgSilhouetteCoefficients, temp_clusterSilhouetteCoefficients
        
      
    return  processAvg, exposureAvg, processSTE,  exposureSTE, avgSilhouetteCoefficients, clusterSilhouetteCoefficients




################################### Dicompose the new signatures into global signatures   #########################
def signature_decomposition(signatures, mtype, directory, genome_build="GRCh37", mutation_context=None):
    
    # open the log file for signature decomposition 
    lognote = open(directory+"/decomposition_logfile.txt", "w") 
    lognote.write("############################ Signature Decomposition Details ################################\n\n\n")
    lognote.write("Context Type: {}\n".format(mtype))
    lognote.write("Genome Build: {}\n".format(genome_build))
    
    paths = cosmic.__path__[0]
    
    if signatures.shape[0]==96:
        
        if genome_build=="GRCh37":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="GRCh37")
            sigDatabase,_,_,_ = read_csv(sigDatabase)
        elif genome_build=="GRCh38":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="GRCh38")
            sigDatabase,_,_,_ = read_csv(sigDatabase) 
        elif genome_build=="mm9":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="mm9")
            sigDatabase,_,_,_ = read_csv(sigDatabase)  
        elif genome_build=="mm10":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="mm10")
            sigDatabase,_,_,_ = read_csv(sigDatabase) 
            
        signames = sigDatabase.columns 
        
        
    elif signatures.shape[0]==78:
        if genome_build=="GRCh37":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="GRCh37").iloc[:,1:]
        elif genome_build=="GRCh38":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="GRCh38").iloc[:,1:]
        elif genome_build=="mm9":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="mm9").iloc[:,1:]
        elif genome_build=="mm10":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="mm10").iloc[:,1:]
            
        signames = sigDatabase.columns
        
    elif signatures.shape[0]==83:
        sigDatabase = pd.read_csv(paths+"/data/sigProfiler_ID_signatures.csv", sep=",").iloc[:,1:]
        signames = sigDatabase.columns
    else:
        sigDatabase = np.random.rand(signatures.shape[0],2)
        signames = list(range(signatures.shape[0]))
        
    sigDatabases = sigDatabase
    letters = list(string.ascii_uppercase)
    letters.extend([i+b for i in letters for b in letters])
    letters = letters[0:signatures.shape[1]]
    
    
    # replace the probability data of the process matrix with the number of mutation
    for i in range(signatures.shape[1]):
        signatures[:, i] =  signatures[:, i]*5000      #(np.sum(exposureAvg[i, :]))
    
    sigDatabase = np.array(sigDatabase)
    allsignatures = np.array([])
    newsig = list() # will create the letter based id of newsignatures
    newsigmatrixidx = list() # will create the original id of newsignature to help to record the original matrix
    fh = open(directory+"/comparison_with_global_ID_signatures.csv", "w")
    fh.write("De novo extracted, Global NMF Signatures, L1 Error %, L2 Error %, Cosine Similarity\n")
    fh.close()
    dictionary = {}
    
    # get the names of denovo signatures
    denovo_signature_names = make_letter_ids(signatures.shape[1], mtype=mutation_context)
    #lognote.write("\n********** Starting Signature Decomposition **********\n\n")
    activity_percentages=[]
    for i, j in zip(range(signatures.shape[1]), denovo_signature_names):
        
        # Only for context SBS96
        if signatures.shape[0]==96:
            lognote = open(directory+"/decomposition_logfile.txt", "a")  
            lognote.write("\n\n\n\n######################## Decomposing "+j+" ########################\n"  )
            lognote.close()
            if genome_build=="mm9" or genome_build=="mm10":
                check_rule_negatives = [1,16]
                check_rule_penalty=1.50
            else:
                check_rule_negatives = []
                check_rule_penalty=1.0
            
            #exposures, _, similarity = ss.add_signatures(sigDatabase, signatures[:,i][:,np.newaxis], presentSignatures=[], solver = "nnls", metric = "l2",check_rule_negatives=check_rule_negatives, check_rule_penalty=check_rule_penalty)
            #print("Exposure after adding", exposures)
            #exposures, _, similarity = ss.remove_all_single_signatures(sigDatabase, exposures, signatures[:,i], metric="l2", solver = "nnls", cutoff=0.01, background_sigs= [], verbose=False)
            #print(exposures)
            _, exposures,L2dist,similarity, cosine_similarity_with_four_signatures = ss.add_remove_signatures(sigDatabase, 
                                                                                                         signatures[:,i], 
                                                                                                         metric="l2", 
                                                                                                         solver="nnls", 
                                                                                                         background_sigs = [0,4], 
                                                                                                         permanent_sigs = [0,4], 
                                                                                                         candidate_sigs="all", 
                                                                                                         allsigids = signames, 
                                                                                                         penalty = 0.05, 
                                                                                                         check_rule_negatives = check_rule_negatives, 
                                                                                                         checkrule_penalty = check_rule_penalty, 
                                                                                                         directory = directory+"/decomposition_logfile.txt", 
                                                                                                         verbose=False)
            #print(exposures)
            #print("######################################################################")
            #ss.remove_all_single_signatures(sigDatabase, exposures, signatures[:,i], metric="cosine", solver = "nnls", cutoff=0.05, background_sigs= [0,4], verbose=True)
            #print("Expousre after remove", exposures)
            #print("\n\n\n\n\n\n\n\n")
        # for other contexts     
        else:
            _, exposures,L2dist,similarity, cosine_similarity_with_four_signatures = ss.add_remove_signatures(sigDatabase, 
                                                                                                         signatures[:,i], 
                                                                                                         metric="l2", 
                                                                                                         solver="nnls", 
                                                                                                         background_sigs = [], 
                                                                                                         candidate_sigs="all", 
                                                                                                         penalty = 0.05, 
                                                                                                         check_rule_negatives = [], 
                                                                                                         checkrule_penalty = [], 
                                                                                                         directory = directory+"/decomposition_logfile.txt", 
                                                                                                         verbose=False)
        #print(signames[np.nonzero(exposures)], similarity)
        
        # calculate the L1 Error %
        L1dist = np.linalg.norm(signatures[:,i]-np.dot(sigDatabase,exposures) , ord=1)/np.linalg.norm(signatures[:,i], ord=1)
        
        #print(exposures[np.nonzero(exposures)]/np.sum(exposures[np.nonzero(exposures)])*100)
        exposure_percentages = exposures[np.nonzero(exposures)]/np.sum(exposures[np.nonzero(exposures)])*100
        listofinformation = list("0"*len(np.nonzero(exposures)[0])*3)
        
        count =0
        decomposed_signatures = []
        contribution_percentages = []
        for j in np.nonzero(exposures)[0]:
            listofinformation[count*3] = signames[j]
            listofinformation[count*3+1] = round(exposure_percentages[count],2)
            contribution_percentages.append(round(exposure_percentages[count],2))
            listofinformation[count*3+2]="%"
            decomposed_signatures.append(signames[j])
            count+=1
        ListToTumple = tuple([mtype, letters[i]]+listofinformation+[L1dist]+[L2dist]+[similarity])
        activity_percentages.append(contribution_percentages)
        
        strings ="Signature %s-%s,"+" Signature %s (%0.2f%s) &"*(len(np.nonzero(exposures)[0])-1)+" Signature %s (%0.2f%s), %0.2f,  %0.2f, %0.2f\n" 
        #print(strings%(ListToTumple))
        ##print(np.nonzero(exposures)[0])
        ##print(similarity)
        ##print("\n")
        #print(strings%(ListToTumple))
        new_signature_thresh_hold = 0.8
        if  similarity>new_signature_thresh_hold and cosine_similarity_with_four_signatures > new_signature_thresh_hold: ########### minimum signtatures and cosine similarity needs to be fitted to become a unique signature 
            allsignatures = np.append(allsignatures, np.nonzero(exposures))
            fh = open(directory+"/comparison_with_global_ID_signatures.csv", "a")
            fh.write(strings%(ListToTumple))
            fh.close()
            
            dictionary.update({"{}".format(mutation_context+letters[i]):decomposed_signatures}) 
            
        else:
            newsig.append(mutation_context+letters[i])
            newsigmatrixidx.append(i)
            fh = open(directory+"/comparison_with_global_ID_signatures.csv", "a")
            fh.write("Signature {}-{}, Signature {}-{}, {}, {}, {}\n".format(mtype, letters[i], mtype, letters[i], 0, 0, 1))
            fh.close()
            dictionary.update({"{}".format(mutation_context+letters[i]):["{}".format(mutation_context+letters[i])]}) 
            #dictionary.update({letters[i]:"Signature {}-{}, Signature {}-{}, {}\n".format(mtype, letters[i], mtype, letters[i], 1 )}) 
       
    
    different_signatures = np.unique(allsignatures)
    different_signatures=different_signatures.astype(int)
    if mtype == "96":
        different_signatures = list(set().union(different_signatures, [0,4]))
        different_signatures.sort()
    
    
        
    # add connected signatures   
    different_signatures = ss.add_connected_sigs(different_signatures, list(signames))
    
    #get the name of the signatures
    try:
        detected_signatures = signames[different_signatures]
        globalsigmats= sigDatabases.loc[:,list(detected_signatures)]
    except:
        detected_signatures=[None]
        globalsigmats=None
    
    newsigsmats=signatures[:,newsigmatrixidx]
    
    #for k, v in dictionary.items():
        #print('{}: {}'.format(k, v))
        
    #only for SBS96
    if mtype == "96" or mtype =="1536":        
        background_sigs = get_indeces(list(detected_signatures), ['SBS1', 'SBS5'])
    
    #for other contexts
    else:
        background_sigs = []
        
    # close the lognote
    lognote.close()
    
    #return values
    return {"globalsigids": list(detected_signatures), "newsigids": newsig, "globalsigs":globalsigmats, "newsigs":newsigsmats/5000, "dictionary": dictionary, 
            "background_sigs": background_sigs, "activity_percentages": activity_percentages} 










################################################### Generation of probabilities for each processes given to A mutation type ############################################
def probabilities(W, H, index, allsigids, allcolnames):  
    
    # setting up the indices 
    rows = index
    cols = allcolnames
    sigs = allsigids
    
    W = np.array(W)
    H= np.array(H)
    # rebuild the original matrix from the estimated W and H 
    genomes = np.dot(W,H)
    
    
    result = 0
    for i in range(H.shape[1]): #here H.shape is the number of sample
        
        M = genomes[:,i][np.newaxis]
        #print (M.shape)
        
        probs = W*H[:,i]/M.T
        #print ("Sum of Probabilities for Sample {} is : ".format(i+1))
        #print(probs)
        #print ("Sum of Probabilities for Sample {} is: ".format(i+1))
        #print(probs.sum(axis=1)[np.newaxis].T) 
        
        
        probs = pd.DataFrame(probs)
        probs.columns = sigs
        col1 = [cols[i]]*len(rows)
        probs.insert(loc=0, column='Sample Names', value=col1)
        probs.insert(loc=1, column='MutationTypes', value = rows)
        #print (probs)
        #print ("\n")
        if i!=0:
            result = pd.concat([result, probs], axis=0)
        else:
            result = probs
    
        
    return result
#########################################################################################################################################################################        





"""        
##########################################################################################################################################################################
#################################################################### Data loading and Result Exporting  ###################################################################################
##########################################################################################################################################################################
"""

########################################################### Functions to load the MATLAB OBJECTS ###################

def extract_arrays(data, field, index=True):
    accumulator = list()
    if index==False:
        for i in data[field]:
            accumulator.append(i)
        return accumulator 
    else:        
        for i in data[field]:
            accumulator.append(i[0][0])
        return accumulator 

def extract_input(data):         
    originalGenome=data['originalGenomes']        
    cancerType = extract_arrays(data, 'cancerType', False)
    sampleNames = extract_arrays(data, 'sampleNames', True)
    subTypes = extract_arrays(data, 'subtypes', True)
    types = extract_arrays(data, 'types', True)
    allFields = [cancerType, originalGenome, sampleNames, subTypes, types]
    return allFields
########################################################################################################################

##################################### Get Input From CSV Files ##############################################

def read_csv(filename, folder = False):
    
  
    if folder==False:
        if type(filename) == str:
            genomes = pd.read_csv(filename, sep=",").iloc[:, :]   
        else:
            genomes = filename
    else:
    
        count = 1
        for file in os.listdir(filename):
            #print(count) 
            df = pd.read_csv(filename+"/"+file)
            if count==1:
                    genomes = df
            else:
                    genomes = pd.merge(genomes, df, on=["Mutation type", "Trinucleotide"])
            count += 1


    

    
    
    mtypes = [str(genomes.shape[0])]
    
    if mtypes == ['96']:
        # create the list to sort the mutation types
        orderlist = ['A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T', 'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T',
                     'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T',
                     'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T', 'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T',
                     'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T',
                     'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T',
                     'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T',
                     'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T',
                     'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T']
        #Contruct the indeces of the matrix
        #setting index and columns names of processAvg and exposureAvg
        index1 = genomes.iloc[:,1]
        index2 = genomes.iloc[:,0]
        index = []
        for i, j in zip(index1, index2):
            index.append(i[0]+"["+j+"]"+i[2])
    
    
        index = np.array(pd.Series(index))
        genomes["index"] = index
    
    
        #sort the mutationa types
        genomes["index"]= pd.Categorical(genomes["index"], orderlist)
        genomes = genomes.sort_values("index")
        
    
        genomes = genomes.iloc[:,2:genomes.shape[1]]
    
        
    
        # set the index 
        genomes = genomes.set_index("index")
    
        # prepare the index and colnames variables
        index = np.array(orderlist)
        colnames = np.array(pd.Series(genomes.columns.tolist()))
    
    else:
        
        index = np.array(genomes.iloc[:,0])
        genomes = genomes.iloc[:,1:]
        genomes = genomes.loc[:, (genomes != 0).any(axis=0)]
        colnames = genomes.columns
       
        
        
   
    
    
    return(genomes, index, colnames, mtypes)


###################################################################################### Export Results ###########################################
def export_information(loopResults, mutation_context, output, index, colnames, wall=False):
    
  
    
    
    # get the number of processes
    i = loopResults[-1]
    #get Wall and Hall
    Wall = loopResults[-3]
    Hall = loopResults[-2]
    
   
    
    # get the mutational contexts    
    #print ("The mutaion type is", mutation_type)    
    m = mutation_context
    if not (m=="DINUC"or m=="INDEL"):
        mutation_type = "SBS"+m
            
    else:
        if m == "DINUC":
            mutation_type = "DBS78"
        elif m== "INDEL":
            mutation_type = "ID83"
    # Create the neccessary directories
    subdirectory = output+"/All_solutions/"+mutation_type+"_"+str(i)+"_Signatures"
    if not os.path.exists(subdirectory):
        os.makedirs(subdirectory)
    
    
    
    #Export the loopResults as pickle objects
    
    #resultname = "signature"+str(i)
    
    # =============================================================================
    #         f = open(output+"/pickle_objects/"+resultname, 'wb')
    #         
    #         pickle.dump(loopResults, f)
    #         f.close()
    # =============================================================================
            
       
    #preparing the column and row indeces for the Average processes and exposures:  
    listOfSignatures = []
    letters = list(string.ascii_uppercase)
    letters.extend([i+b for i in letters for b in letters])
    letters = letters[0:i]
    
    for j,l in zip(range(i),letters)  :
        listOfSignatures.append(mutation_type+l)
    listOfSignatures = np.array(listOfSignatures)
    
    #print("print listOfSignares ok", listOfSignatures)
        
    
    #Extract the genomes, processAVG, processStabityAvg
    genome= loopResults[0]
    #print ("genomes are ok", genome)
    processAvg= (loopResults[1])
    exposureAvg= (loopResults[2])
    process_stabililities = np.array(loopResults[6])
    minProcessStability= round(np.min(process_stabililities), 2) 
    meanProcessStability = round(np.mean(process_stabililities), 2)
    #print ("processStabityAvg is ok", processStabityAvg)
    
    # Calculating and listing the reconstruction error, process stability and signares to make a csv file at the end
    reconstruction_error = round(LA.norm(genome-np.dot(processAvg, exposureAvg), 'fro')/LA.norm(genome, 'fro'), 4)
    
    #print ("reconstruction_error is ok", reconstruction_error)
    #print (' Initial reconstruction error is {} and the process stability is {} for {} signatures\n\n'.format(reconstruction_error, round(processStabityAvg,4), i))
    # Preparing the results to export as textfiles for each signature
    
    #First exporting the Average of the processes
    processAvg= pd.DataFrame(processAvg)
    processes = processAvg.set_index(index)
    processes.columns = listOfSignatures
    processes = processes.rename_axis("MutationType", axis="columns")
    #print(processes)
    #print("process are ok", processes)
    processes.to_csv(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Signatures"+".txt", "\t", index_label=[processes.columns.name]) 
    
    #Second exporting the Average of the exposures
    exposureAvg = pd.DataFrame(exposureAvg.astype(int))
    exposures = exposureAvg.set_index(listOfSignatures)
    exposures.columns = colnames
    #plot exposures
    """for p in range(exposures.shape[0]):
      plt.bar(colnames, exposures.iloc[p], bottom = np.sum(exposures.iloc[:p], axis = 0), label = listOfSignatures[p])
        
    plt.legend(loc=(1.01,0.0))
    plt.title("Signature Activities on Samples")
    plt.xlabel("Samples")
    plt.ylabel("Mutation Count")
    plt.xticks(colnames, rotation='vertical')
    plt.tight_layout()
    plt.savefig(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Activities_Plot.pdf", dpi=300)  
    
    plt.close()"""
    exposures = exposures.T
    exposures = exposures.rename_axis("Samples", axis="columns")
    #print("exposures are ok", exposures)
    exposures.to_csv(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Activities.txt", "\t", index_label=[exposures.columns.name]) 
    
    
    # get the standard errors of the processes
    processSTE = loopResults[3]      
    #export the processStd file
    processSTE = pd.DataFrame(processSTE)
    processSTE = processSTE.set_index(index)
    processSTE.columns = listOfSignatures
    processSTE = processSTE.rename_axis("MutationType", axis="columns")
    processSTE.to_csv(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Signatures_SEM_Error"+".txt", "\t", float_format='%.2E', index_label=[processes.columns.name]) 
    
    # get the standard errors of the exposures   
    exposureSTE = loopResults[4] 
    #export the exposureStd file
    exposureSTE = pd.DataFrame(exposureSTE)
    exposureSTE = exposureSTE.set_index(listOfSignatures)
    exposureSTE.columns = colnames
    exposureSTE = exposureSTE.T
    exposureSTE = exposureSTE.rename_axis("Samples", axis="columns")
    exposureSTE.to_csv(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Activities_SEM_Error.txt", "\t", float_format='%.2E', index_label=[exposures.columns.name]) 
    
    all_similarities = loopResults[8].copy() 
    #all_similarities = all_similarities.replace(None, 1000)
    all_similarities['L1_Norm_%'] = all_similarities['L1_Norm_%'].astype(str) + '%'
    all_similarities['L2_Norm_%'] = all_similarities['L2_Norm_%'].astype(str) + '%'
    all_similarities.to_csv(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Samples_stats.txt", sep="\t")
    
    signature_stats =  loopResults[9]               
    signature_stats = signature_stats.set_index(listOfSignatures)
    signature_stats = signature_stats.rename_axis("Signatures", axis="columns")
    signature_stats.to_csv(subdirectory+"/"+mutation_type+"_S"+str(i)+"_"+"Signatures_stats.txt", "\t", index_label="Signatures") 
    
    #export convergence information 
    converge_information = loopResults[13]
    converge_information = pd.DataFrame(np.around(converge_information, decimals=3))
    conv_index = list(range(1,len(converge_information)+1)) 
    colmetrices = ['L1', 'L1 %', 'L2', 'L2 %', 'KL Divergence', "Convergence Iterations"]
    converge_information.index = conv_index 
    converge_information.columns = colmetrices
    converge_information.to_csv(subdirectory+"/"+mutation_type+"_S"+str(i)+"_"+"NMF_Convergence_Information.txt", "\t", index_label="NMF_Replicate")
    
    # export Wall and Hall if "wall" argument is true
    if wall==True:
        
        np.savetxt(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Wall.txt", Wall, delimiter="\t")
        np.savetxt(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Hall.txt", Hall, delimiter="\t")
    
    fh = open(output+"/All_solutions_stat.csv", "a") 
    print ('The reconstruction error is {}, average process stability is {} and \nthe minimum process stability is {} for {} signatures\n\n'.format(reconstruction_error,  round(meanProcessStability,2), round(minProcessStability,2), i))
    fh.write('{}, {}, {}, {}\n'.format(i , round(minProcessStability, 2), round(reconstruction_error, 4), round(meanProcessStability,2)))
    fh.close()
    
    
    
    ########################################### PLOT THE SIGNATURES ################################################
    #prepare the texts lists:
    stability_list = signature_plotting_text(loopResults[6], "Stability", "float")
    total_mutation_list = signature_plotting_text(loopResults[7], "Sig. Mutations", "integer")
    
    
    
    if m=="DINUC" or m=="78":        
        plot.plotDBS(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Signatures"+".txt", subdirectory+"/Signature_plot" , "S"+str(i), "78", True, custom_text_upper=stability_list, custom_text_middle=total_mutation_list)
    elif m=="INDEL" or m=="83":
        plot.plotID(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Signatures"+".txt", subdirectory+"/Signature_plot" , "S"+str(i), "83", True, custom_text_upper=stability_list, custom_text_middle=total_mutation_list)
    else:
        plot.plotSBS(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Signatures"+".txt", subdirectory+"/Signature_plot", "S"+str(i), m, True, custom_text_upper=stability_list, custom_text_middle=total_mutation_list)
     
        
# =============================================================================
#     processAvg = pd.read_csv(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Signatures_"+".txt", sep="\t", index_col=0)
#     exposureAvg = pd.read_csv(subdirectory+"/"+mutation_type+"_S"+str(i)+"_Sig_activities.txt", sep="\t", index_col=0)
#     probability = probabilities(processAvg, exposureAvg)
#     probability=probability.set_index("Sample")
#     probability.to_csv(subdirectory+"/mutation_probabilities.txt", "\t") 
# =============================================================================
    

#############################################################################################################
######################################## MAKE THE FINAL FOLDER ##############################################
#############################################################################################################
def make_final_solution(processAvg, allgenomes, allsigids, layer_directory, m, index, allcolnames, process_std_error = "none", signature_stabilities = " ", \
                        signature_total_mutations= " ", signature_stats = "none",  remove_sigs=False, attribution= 0, denovo_exposureAvg  = "none", penalty=0.05, \
                        background_sigs=0, genome_build="GRCh37",  verbose=False):
    
    # Get the type of solution from the last part of the layer_directory name
    solution_type = layer_directory.split("/")[-1]
    
    # Create the lognote file
    lognote = open(layer_directory+"/Signature_assaignment_logfile.txt", "w")
    lognote.write("************************ Stepwise Description of Signature Assaignment to Samples ************************")
    lognote.close()
    
    
    #Get the type of Signatures
    if m == 83:
        signature_type = "INDEL"
    elif m==78:
        signature_type = "DINUC"
    else:
        signature_type = "SBS"+str(m)
    
    
    
    allgenomes = np.array(allgenomes)
    if (m=="96" or m=="1536" or m=="288") and (genome_build=="mm9" or genome_build=="mm10"):
        check_rule_negatives = [1,16]
        check_rule_penalty=1.50
    else:
        check_rule_negatives = []
        check_rule_penalty=1.0
    exposureAvg = np.zeros([processAvg.shape[1], allgenomes.shape[1]] )  
    if remove_sigs==True:
        denovo_exposureAvg = denovo_exposureAvg.T
        
        #print("\n")
        for r in range(allgenomes.shape[1]):
            if verbose==True:
                print("\n\n\n\n\n                                        ################ Sample "+str(r+1)+ " #################")
                     
            # Record information to lognote
            lognote = open(layer_directory+"/Signature_assaignment_logfile.txt", "a")
            lognote.write("\n\n\n\n\n                    ################ Sample "+str(r+1)+ " #################\n") 
            
            sample_exposure = np.array(denovo_exposureAvg.iloc[:,r])
            
            init_sig_idx = np.nonzero(sample_exposure)[0]
            init_sigs = denovo_exposureAvg.index[init_sig_idx]
            
            
            init_decomposed_sigs = []
            for de_novo_sig in init_sigs:
                
                init_decomposed_sigs = union(init_decomposed_sigs, list(attribution[de_novo_sig]))
                
            #print(init_decomposed_sigs) 
            init_decomposed_sigs_idx = get_indeces(allsigids, init_decomposed_sigs)
            init_decomposed_sigs_idx.sort()
            init_decomposed_sigs_idx = list(set().union(init_decomposed_sigs_idx, background_sigs))
            #print(init_decomposed_sigs_idx)
            
            # get the indices of the background sigs in the initial signatures
            background_sig_idx = get_indeces(init_decomposed_sigs_idx, background_sigs)
            
            
            fit_signatures = processAvg[:,init_decomposed_sigs_idx]
            #fit signatures
            newExposure, newSimilarity = ss.fit_signatures(fit_signatures, allgenomes[:,r])
            
            #create the exposureAvg vector
            #print(init_decomposed_sigs_idx)
            #print(newExposure)
            for nonzero_idx, nozero_exp in zip(init_decomposed_sigs_idx, newExposure):
                exposureAvg[nonzero_idx, r] = nozero_exp
            if verbose==True:
                print("############################# Initial Composition #################################### ") 
                print(pd.DataFrame(exposureAvg[:, r],  index=allsigids).T)   
                print("L2%: ", newSimilarity)  
                
            lognote.write("############################# Initial Composition ####################################\n")
            exposures = pd.DataFrame(exposureAvg[:, r],  index=allsigids).T
            lognote.write("{}\n".format(exposures.iloc[:,exposures.to_numpy().nonzero()[1]])) 
            lognote.write("L2 Error %: {}\nCosine Similarity: {}\n".format(round(newSimilarity,2), round(cos_sim(allgenomes[:,r], np.dot(processAvg, exposureAvg[:, r] )),2)))
            #remove signatures 
            exposureAvg[:,r],L2dist,cosine_sim = ss.remove_all_single_signatures(processAvg, exposureAvg[:, r], allgenomes[:,r], metric="l2", \
                       solver = "nnls", cutoff=0.05, background_sigs= [], verbose=False)
            if verbose==True:
                print("############################## Composition After Initial Remove ############################### ")
                print(pd.DataFrame(exposureAvg[:, r],  index=allsigids).T)  
                print("L2%: ", L2dist)
            lognote.write("############################## Composition After Initial Remove ###############################\n")
            exposures = pd.DataFrame(exposureAvg[:, r],  index=allsigids).T
            lognote.write("{}\n".format(exposures.iloc[:,exposures.to_numpy().nonzero()[1]])) 
            lognote.write("L2 Error %: {}\nCosine Similarity: {}\n".format(round(L2dist,2), round(cosine_sim,2)))
            lognote.write("\n############################## Performing Add-Remove Step ##############################\n")
            #Close the Lognote file
            lognote.close()
            
            init_add_sig_idx = list(set().union(list(np.nonzero(exposureAvg[:, r])[0]), background_sigs))
            #print(init_add_sig_idx)
            
            #get the background_sig_idx for the add_remove function only for the decomposed solution:
            if background_sigs != 0:  # in the decomposed solution only 
                background_sig_idx = get_indeces(allsigids, ["SBS1", "SBS5"])
                #print(background_sig_idx)
            
             
            
            
            # if the there is no other signatures to be added on top the existing signatures
            try:
                
                _, exposureAvg[:, r],L2dist,similarity, cosine_similarity_with_four_signatures = ss.add_remove_signatures(processAvg, 
                                                                                                  allgenomes[:,r], 
                                                                                                  metric="l2", 
                                                                                                  solver="nnls", 
                                                                                                  background_sigs = init_add_sig_idx, 
                                                                                                  permanent_sigs = background_sig_idx, 
                                                                                                  candidate_sigs="all", 
                                                                                                  allsigids = allsigids, 
                                                                                                  penalty = 0.05, 
                                                                                                  check_rule_negatives = check_rule_negatives, 
                                                                                                  checkrule_penalty = check_rule_penalty, 
                                                                                                  directory = layer_directory+"/Signature_assaignment_logfile.txt", 
                                                                                                  verbose=False)
                if verbose==True:
                    print("####################################### Composition After Add-Remove #######################################\n") 
                    print(exposureAvg[:, r])
                    print("L2%: ", L2dist)
                # Recond the information in the log file
                lognote = open(layer_directory+"/Signature_assaignment_logfile.txt", "a")
                lognote.write("####################################### Composition After Add-Remove #######################################\n")
                exposures = pd.DataFrame(exposureAvg[:, r],  index=allsigids).T
                lognote.write("{}\n".format(exposures.iloc[:,exposures.to_numpy().nonzero()[1]])) 
                lognote.write("L2 Error %: {}\nCosine Similarity: {}\n".format(round(L2dist,2), round(similarity,2)))
                lognote.close()
            except:
                pass
            
            
            """
            # add signatures
            exposureAvg[:, r], _, similarity = ss.add_signatures(processAvg, allgenomes[:,r][:,np.newaxis], presentSignatures=copy.deepcopy(init_add_sig_idx),cutoff=penalty, metric="l2", solver = "nnls",check_rule_negatives=check_rule_negatives, check_rule_penalty=check_rule_penalty)
            if verbose==True:
                print("############################################################# After adding :")
                print(exposureAvg[:, r])
            #print("\n")
            #remove signatures 
            exposureAvg[:,r],_,_ = ss.remove_all_single_signatures(processAvg, exposureAvg[:, r], allgenomes[:,r], metric="l2", \
                       solver = "nnls", cutoff=0.01, background_sigs= background_sig_idx, verbose=False)
            if verbose==True:
                print("############################################################# After Remove :")
                print(exposureAvg[:, r])
            
            init_add_sig_idx = list(set().union(list(np.nonzero(exposureAvg[:, r])[0]), background_sigs))
            #print(init_add_sig_idx)
            
            # add signatures
            exposureAvg[:, r], _, similarity = ss.add_signatures(processAvg, allgenomes[:,r][:,np.newaxis], presentSignatures=copy.deepcopy(init_add_sig_idx),cutoff=penalty, metric="l2", solver = "nnls", check_rule_negatives=check_rule_negatives, check_rule_penalty=check_rule_penalty)
            
            if verbose==True:
                print("############################################################# After adding :") 
                print(exposureAvg[:, r])
            """      
               
    else:   
        
        # when refilt de_novo_signatures 
        if denovo_exposureAvg=="none":
            for g in range(allgenomes.shape[1]):
                
                exposures, _, similarity = ss.add_signatures(processAvg, 
                                                             allgenomes[:,g][:,np.newaxis], 
                                                             presentSignatures=[],
                                                             cutoff=penalty,
                                                             check_rule_negatives=[], 
                                                             check_rule_penalty=1.0)
                exposureAvg[:,g] = exposures
                
        # when use the exposures from the initial NMF
        else:
            
            exposureAvg=denovo_exposureAvg
            
        
    
    processAvg= pd.DataFrame(processAvg.astype(float))
    processes = processAvg.set_index(index)
    processes.columns = allsigids
    processes = processes.rename_axis("MutationsType", axis="columns")
    processes.to_csv(layer_directory+"/"+solution_type+"_"+"Signatures_"+signature_type+".txt", "\t", float_format='%.8f',index_label=[processes.columns.name]) 
    
     
    exposureAvg = pd.DataFrame(exposureAvg.astype(int))
    
    allsigids = np.array(allsigids)
    exposures = exposureAvg.set_index(allsigids)
    exposures.columns = allcolnames
    
    #plot exposures
    for p in range(exposures.shape[0]):
      plt.bar(allcolnames, exposures.iloc[p], bottom = np.sum(exposures.iloc[:p], axis = 0), label = allsigids[p])
        
    plt.legend(loc=(1.01,0.0))
    plt.title("Signature Activities on Samples")
    plt.xlabel("Samples")
    plt.ylabel("Mutation Count")
    plt.xticks(allcolnames, rotation='vertical')
    try:
        plt.tight_layout()
    except: 
        pass
    plt.savefig(layer_directory+"/"+solution_type+"_"+"Activities_Plot_"+signature_type+".pdf", dpi=300)  
    
    plt.close()
    
    exposures = exposures.T
    exposures = exposures.rename_axis("Samples", axis="columns")
    exposures.to_csv(layer_directory+"/"+solution_type+"_"+"Activities_"+signature_type+".txt", "\t", index_label=[exposures.columns.name]) 
    
    
    
    # Calcutlate the similarity matrices
    est_genomes = np.dot(processAvg, exposureAvg)
    all_similarities, cosine_similarities = calculate_similarities(allgenomes, est_genomes, allcolnames)
    all_similarities.iloc[:,[3,5]] = all_similarities.iloc[:,[3,5]].astype(str) + '%'
    all_similarities.to_csv(layer_directory+"/"+solution_type+"_Samples_stats_"+signature_type+".txt", sep="\t")
    
    
    if type(process_std_error) != str:
        process_std_error= pd.DataFrame(process_std_error)
        processSTE = process_std_error.set_index(index)
        processSTE.columns = allsigids
        processSTE = processSTE.rename_axis("MutationType", axis="columns")
        processSTE.to_csv(layer_directory+"/"+solution_type+"_"+"Signatures_SEM_Error_"+signature_type+".txt", "\t", float_format='%.2E', index_label=[processes.columns.name]) 
    
    if type(signature_stats)!=str:
        signature_stats = signature_stats.set_index(allsigids)
        signature_stats = signature_stats.rename_axis("Signatures", axis="columns")
        signature_stats.to_csv(layer_directory+"/"+solution_type+"_"+"Signatures_stats_"+signature_type+".txt", "\t", index_label=[exposures.columns.name]) 
        signature_total_mutations = np.sum(exposureAvg, axis =1).astype(int)
        signature_total_mutations = signature_plotting_text(signature_total_mutations, "Sig. Mutations", "integer")
    else: #when it works with the decomposed solution
        signature_total_mutations = np.sum(exposureAvg, axis =1).astype(int)
        signature_total_mutations = signature_plotting_text(signature_total_mutations, "Sig. Mutations", "integer")
        if m == "1536" or m=="288": # collapse the 1536 to 96
            m = "96"  
       
    ########################################### PLOT THE SIGNATURES ################################################
    
    if m=="DINUC" or m=="78":
        plot.plotDBS(layer_directory+"/"+solution_type+"_"+"Signatures_"+signature_type+".txt", layer_directory+"/Signature_plot" , solution_type, "78", True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations )
        #plotE.plotDBS(layer_directory+"/"+solution_type+"_"+"Signatures_"+signature_type+".txt", layer_directory+"/Signature_plot" , solution_type, "78", True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations )
        
    elif m=="INDEL" or m=="83":
        plot.plotID(layer_directory+"/"+solution_type+"_"+"Signatures_"+signature_type+".txt", layer_directory+"/Signature_plot" , solution_type, "94", True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations )
        #plotE.plotID(layer_directory+"/"+solution_type+"_"+"Signatures_"+signature_type+".txt", layer_directory+"/Signature_plot" , solution_type, "94", True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations )
    else:
        plot.plotSBS(layer_directory+"/"+solution_type+"_"+"Signatures_"+signature_type+".txt", layer_directory+"/Signature_plot", solution_type, m, True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations )
        #plotE.plotSBS(layer_directory+"/"+solution_type+"_"+"Signatures_"+signature_type+".txt", layer_directory+"/Signature_plot", solution_type, m, True, custom_text_upper= signature_stabilities, custom_text_middle = signature_total_mutations, path=layer_directory )
        
        
        
    #processAvg = pd.read_csv(layer_directory+"/"+solution_type+"_"+"Signatures.txt", sep="\t", index_col=0)
    #exposureAvg = pd.read_csv(layer_directory+"/"+solution_type+"_"+"Activities.txt", sep="\t", index_col=0)
    
    probability = probabilities(processAvg, exposureAvg, index, allsigids, allcolnames)
    probability=probability.set_index("Sample Names" )
    probability.to_csv(layer_directory+"/Mutation_Probabilities.txt", "\t") 
    
    """
    try:
        clusters = dendrogram(exposureAvg, 0.05, layer_directory)
        clusters.to_csv(layer_directory+"/Cluster_of_Samples.txt", "\t") 
    except:
        pass
    """
    
    return exposures
"""
#############################################################################################################
######################################### PLOTTING FUNCTIONS ##############################################
#############################################################################################################
"""

######################################## Clustering ####################################################
def dendrogram(data, threshold, layer_directory):
    colnames = data.columns
    data = np.array(data)

    Z = hierarchy.linkage(data.T, 'single',  'cosine')
    plt.figure(figsize=(15, 9))
    dn = hierarchy.dendrogram(Z, labels = colnames, color_threshold=threshold)
    plt.title("Clustering of Samples Based on Mutational Signatures" )
    plt.ylabel("Cosine Distance")
    plt.xlabel("Sample IDs")
    #plt.ylim((0,1))
    plt.savefig(layer_directory+'/dendrogram.pdf',figsize=(10, 8), dpi=300)
    # which datapoints goes to which cluster
    # The indices of the datapoints will be displayed as the ids 
    Y = hierarchy.fcluster(Z, threshold, criterion='distance', R=None, monocrit=None)
    dataframe = pd.DataFrame({"Cluster":Y, "Sample Names":list(colnames)})
    dataframe = dataframe.set_index("Sample Names")
    #print(dataframe)
    dictionary = {"clusters":Y, "informations":dn}
    
    return dataframe 


######################################## Plot the reconstruction error vs stabilities and select the optimum number of signature ####################################################   
def stabVsRError(csvfile, output, title, all_similarities_list, mtype= ""):
    
    
    data = pd.read_csv(csvfile, sep=",")
    
    
    #exracting and preparing other similiry matrices from the all_similarities_list
    mean_l1 = []; maximum_l1 = []; mean_l2 = []; maximum_l2 = []; mean_kl = []; maximum_kl = []; wilcoxontest=[]; #all_mean_l2=[];
    #median_l1= Median L1 , maximum _l1 = Maximum L1, median_l2= Median L2 , maximum _l2 = Maximum L2, median_kl= Median KL , maximum _kl = Maximum KL, wilcoxontest = Wilcoxontest significance (True or False); all_med_l2=[] = list of all Median L2
    
    
    #statistical test we need to set a previous L2 (pre_l2) which have a same length as the sample number
    pre_l2_dist = np.zeros([all_similarities_list[0].shape[0]])
    pre_mean = np.inf  
    
    for values in range(len(all_similarities_list)): # loop through the opposite direction
      all_similarities = all_similarities_list[values].iloc[:,[3,5,6]]
      avg_stability=data["avgStability"][values]
      min_stability=data["Stability"][values]
      thresh_hold=avg_stability+min_stability
      
      #record the statistical test between the l2_of the current and previous signatures first
      cur_l2_dist = all_similarities["L2_Norm_%"]
      cur_mean = all_similarities["L2_Norm_%"].mean()
     
      wiltest = ranksums(np.array(cur_l2_dist), np.array(pre_l2_dist))[1]
      
      
      if (wiltest<0.05) and (pre_mean-cur_mean>0):
          wilcoxontest.append("True")
          #all_mean_l2.append(cur_mean)
          #pre_l2_dist = cur_l2_dist
          #pre_mean = cur_mean
          if thresh_hold>=1.0:
              pre_l2_dist = cur_l2_dist
              pre_mean = cur_mean
      else:
          wilcoxontest.append("False")
          #all_mean_l2.append(pre_mean)
          #pre_l2_dist = cur_l2_dist
      
      
      mean_l1.append(round(all_similarities["L1_Norm_%"].mean(),2))
      maximum_l1.append(round(all_similarities["L1_Norm_%"].max(),2))
      mean_l2.append(round(all_similarities["L2_Norm_%"].mean(),2))
      maximum_l2.append(round(all_similarities["L2_Norm_%"].max(),2))
      mean_kl.append(round(all_similarities["KL Divergence"].mean(), 4))
      maximum_kl.append(round(all_similarities["KL Divergence"].max(), 4))
      
     
    data.iloc[:,2] = np.round(data.iloc[:,2]*100, 2)
    data = data.assign(**{'Mean Sample L1%': mean_l1, 'Maximum Sample L1%': maximum_l1, 'Mean Sample L2%': mean_l2, 'Maximum Sample L2%': maximum_l2, 'Significant Decrease of L2':wilcoxontest, 'Mean Sample KL': mean_kl, 'Maximum Sample KL': maximum_kl})  
    data=data.rename(columns = {'Stability': 'Minimum Stability'})
    data = data.set_index("Total Signatures")
     
    
    #get the solution
    probable_solutions = data.copy()
    avg_stability = data.iloc[:,2]
    data = data.drop(columns="avgStability")
    
    probable_solutions["Significant Decrease of L2"] = probable_solutions["Significant Decrease of L2"].astype(str)
    probable_solutions = probable_solutions[probable_solutions["Significant Decrease of L2"]=="True"]
    #probable_solutions = probable_solutions[probable_solutions["avgStability"]>=0.8]
    if mtype=="DBS78" or mtype=="ID83":
        probable_solutions = probable_solutions[probable_solutions["Minimum Stability"]>=0.80]
    else:    
        probable_solutions = probable_solutions[probable_solutions["Minimum Stability"]+probable_solutions["avgStability"]>=1.0]
    #probable_solutions = probable_solutions[probable_solutions["avgStability"]>=0.8]
    #print(probable_solutions.columns)
    #print(probable_solutions)
    #probable_solutions = probable_solutions[probable_solutions["Minimum Stability"]+probable_solutions["avgStability"]>=1.25]
    
    
    try: 
        alternative_solution = int(probable_solutions.index[-1])
    except:
        alternative_solution = int(data.index[0])
        print("There is no solution over the thresh-hold minimum stability. We are selecting the minimum number of signature which could be wrong.")
#---------------------------- alternative solution end -------------------------------------------------#
    
   
    # Create some mock data
    
    t = np.array(data.index)
    
    #data1 = np.array(data.iloc[:,2])  #reconstruction error
    data1 = np.array(mean_l2)/100
    data2 = np.array(avg_stability)  #process stability
    

    
    
    shadow_alternative_start = alternative_solution-0.2
    shadow_alternative_end=alternative_solution+0.2
    
    
    fig, ax1 = plt.subplots(num=None, figsize=(10, 6), dpi=300, facecolor='w', edgecolor='k')
    
    color = 'tab:red'
    ax1.set_xlabel('Total Signatures')
    ax1.set_ylabel('Mean L2 %', color=color)
    ax1.set_title(title)
    lns1 = ax1.plot(t, data1, marker='o', linestyle=":", color=color, label = 'Mean L2 %')
    
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.xaxis.set_ticks(np.arange(min(t), max(t)+1, 1))
    #ax1.axvspan(shadow_start, shadow_end, alpha=0.20, color='#ADD8E6')
    ax1.axvspan(shadow_alternative_start,  shadow_alternative_end, alpha=0.20, color='#696969')         
    # manipulate the y-axis values into percentage 
    vals = ax1.get_yticks()
    ax1.set_xticklabels(np.arange(min(t), max(t)+1, 1),list(), rotation=30)
    
    ax1.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    
    #ax1.legend(loc=0)
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:blue'
    ax2.set_ylabel('Avg Stability', color=color)  # we already handled the x-label with ax1
    lns2 = ax2.plot(t, data2, marker='s', linestyle="-.", color=color, label = 'Avg Stability')
    ax2.tick_params(axis='y', labelcolor=color)
    #ax2.legend(loc=1)
    
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    #plt.show()
    
    # added these three lines
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=0)
    
    plt.savefig(output+'/'+mtype+'_selection_plot.pdf')    
    
    plt.close()
    
    #put * in the selected solution
    index = list(data.index.astype(str))
    solution = get_indeces(index, [str(alternative_solution)])[0]
    index[solution] = index[solution]+"*" 
    data.index = index
    
    # add % signs
    
    data.iloc[:,1:6] = data.iloc[:,1:6].astype(str) + '%'
    return alternative_solution, data
######################################## Plot Samples ####################################################
def plot_csv_sbs_samples(filename, output_path, project, mtype="96", percentage=False, custom_text_upper=" " ):

    
    
    data, index, colnames, _ = read_csv(filename, folder = False)
    data.index.names= ["MutationType"]
    data.to_csv("new_file.text", sep="\t")

    plot.plotSBS("new_file.text", output_path, project, mtype, False, custom_text_upper=" ")
    
    os.remove("new_file.text")
    
#plot_csv_sbs_samples("CNS-Oligo.96.csv" , "/Users/mishugeb/Desktop/new_plot", "project", mtype="96", percentage=False, custom_text_upper=" " )
    
    
    
########## Functions From Scikit Learn #################
from sklearn.utils.extmath import randomized_svd, squared_norm
from sklearn.utils.validation import check_non_negative
from sklearn.utils import check_random_state
from math import sqrt


def norm(x):
    """Dot product-based Euclidean norm implementation
    See: http://fseoane.net/blog/2011/computing-the-vector-norm/
    Parameters
    ----------
    x : array-like
        Vector for which to compute the norm
    """
    return sqrt(squared_norm(x))

def initialize_nmf(X, n_components, init=None, eps=1e-6,
                    random_state=None):
    """Algorithms for NMF initialization.
    Computes an initial guess for the non-negative
    rank k matrix approximation for X: X = WH
    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        The data matrix to be decomposed.
    n_components : integer
        The number of components desired in the approximation.
    init :  None | 'random' | 'nndsvd' | 'nndsvda' | 'nndsvdar'
        Method used to initialize the procedure.
        Default: None.
        Valid options:
        - None: 'nndsvd' if n_components <= min(n_samples, n_features),
            otherwise 'random'.
        - 'random': non-negative random matrices, scaled with:
            sqrt(X.mean() / n_components)
        - 'nndsvd': Nonnegative Double Singular Value Decomposition (NNDSVD)
            initialization (better for sparseness)
        - 'nndsvda': NNDSVD with zeros filled with the average of X
            (better when sparsity is not desired)
        - 'nndsvdar': NNDSVD with zeros filled with small random values
            (generally faster, less accurate alternative to NNDSVDa
            for when sparsity is not desired)
        - 'custom': use custom matrices W and H
    eps : float
        Truncate all values less then this in output to zero.
    random_state : int, RandomState instance or None, optional, default: None
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`. Used when ``random`` == 'nndsvdar' or 'random'.
    Returns
    -------
    W : array-like, shape (n_samples, n_components)
        Initial guesses for solving X ~= WH
    H : array-like, shape (n_components, n_features)
        Initial guesses for solving X ~= WH
    References
    ----------
    C. Boutsidis, E. Gallopoulos: SVD based initialization: A head start for
    nonnegative matrix factorization - Pattern Recognition, 2008
    http://tinyurl.com/nndsvd
    """
    check_non_negative(X, "NMF initialization")
    n_samples, n_features = X.shape

    if (init is not None and init != 'random'
            and n_components > min(n_samples, n_features)):
        raise ValueError("init = '{}' can only be used when "
                         "n_components <= min(n_samples, n_features)"
                         .format(init))

    if init is None:
        if n_components <= min(n_samples, n_features):
            init = 'nndsvd'
        else:
            init = 'random'

    # Random initialization
    if init == 'random':
        avg = np.sqrt(X.mean() / n_components)
        rng = check_random_state(random_state)
        H = avg * rng.randn(n_components, n_features)
        W = avg * rng.randn(n_samples, n_components)
        # we do not write np.abs(H, out=H) to stay compatible with
        # numpy 1.5 and earlier where the 'out' keyword is not
        # supported as a kwarg on ufuncs
        np.abs(H, H)
        np.abs(W, W)
        return W, H

    # NNDSVD initialization
    U, S, V = randomized_svd(X, n_components, random_state=random_state)
    W, H = np.zeros(U.shape), np.zeros(V.shape)

    # The leading singular triplet is non-negative
    # so it can be used as is for initialization.
    W[:, 0] = np.sqrt(S[0]) * np.abs(U[:, 0])
    H[0, :] = np.sqrt(S[0]) * np.abs(V[0, :])

    for j in range(1, n_components):
        x, y = U[:, j], V[j, :]

        # extract positive and negative parts of column vectors
        x_p, y_p = np.maximum(x, 0), np.maximum(y, 0)
        x_n, y_n = np.abs(np.minimum(x, 0)), np.abs(np.minimum(y, 0))

        # and their norms
        x_p_nrm, y_p_nrm = norm(x_p), norm(y_p)
        x_n_nrm, y_n_nrm = norm(x_n), norm(y_n)

        m_p, m_n = x_p_nrm * y_p_nrm, x_n_nrm * y_n_nrm

        # choose update
        if m_p > m_n:
            u = x_p / x_p_nrm
            v = y_p / y_p_nrm
            sigma = m_p
        else:
            u = x_n / x_n_nrm
            v = y_n / y_n_nrm
            sigma = m_n

        lbd = np.sqrt(S[j] * sigma)
        W[:, j] = lbd * u
        H[j, :] = lbd * v

    W[W < eps] = 0
    H[H < eps] = 0

    if init == "nndsvd":
        pass
    elif init == "nndsvda":
        avg = np.mean(X)
        W[W == 0] = avg
        H[H == 0] = avg
    elif init == "nndsvdar":
        rng = check_random_state(random_state)
        #avg = X.mean()
        avg = np.mean(X)
        W[W == 0] = abs(avg * rng.randn(len(W[W == 0])) / 100)
        H[H == 0] = abs(avg * rng.randn(len(H[H == 0])) / 100)
    elif init == "alexandrov-lab-custom":
         min_X = np.min(X[X>0])
         H[H <= min_X] = min_X
         W[W <= min_X] = min_X
    else:
        raise ValueError(
            'Invalid init parameter: got %r instead of one of %r' %
            (init, (None, 'random', 'nndsvd', 'nndsvda', 'nndsvdar')))

    return W, H    
