#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 10:55:51 2019

@author: mishugeb
"""
import numpy as np
from numpy import linalg as LA
from scipy.optimize import nnls
from scipy.optimize import minimize
from sigproextractor import subroutines as sub
import copy 
"""
#############################################################################################################
#################################### Functions For Single Sample Algorithms #############################
#############################################################################################################
"""
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



def parameterized_objective2_custom(x, signatures, samples):
    rec = np.dot(signatures, x)
    try:
        y = LA.norm(samples-rec[:,np.newaxis])
    except:
        y = LA.norm(samples-rec)
    return y


def constraints1(x, samples):
    sumOfSamples=np.sum(samples, axis=0)
    #print(sumOfSamples)
    #print(x)
    result = sumOfSamples-(np.sum(x))
    #print (result)
    return result


def create_bounds(idxOfZeros, samples, numOfSignatures):
    total = np.sum(samples)
    b = (0.0, float(total))
    lst =[b]*numOfSignatures
    
    for i in idxOfZeros:
        lst[i] = (0.0, 0.0) 
    
    return lst 

# required for removing signatures with background
def get_changed_background_sig_idx(exposures, background_sigs):
    
    background_sigs_values = sub.get_items_from_index(exposures,background_sigs)
    temp_exposures = exposures[:]
    temp_exposures[:] = (value for value in temp_exposures if value != 0)
    
    # remove the background signatures with zero values
    background_sigs[:] = (value for value in background_sigs_values if value != 0)
    
    # get the new indices of the background signatures 
    background_sigs = sub.get_indeces(temp_exposures, background_sigs_values)
    
    return background_sigs


# Fit signatures 
def fit_signatures(W, genome, metric="l2"):
    genome = np.array(genome)
    maxmutation = round(np.sum(genome))
     
    
    
    #initialize the guess
    x0 = np.random.rand(W.shape[1], 1)*maxmutation
    x0= x0/np.sum(x0)*maxmutation
    
    #the optimization step
    #set the bounds and constraints
    bnds = create_bounds([], genome, W.shape[1]) 
    cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genome]}
    
    
    #sol = scipy.optimize.minimize(parameterized_objective2_custom, x0, args=(W, genome),  bounds=bnds, constraints =cons1, tol=1e-30)
    #solution = sol.x
    
    
    ### using NNLS algorithm 
    reg = nnls(W,genome)
    weights = reg[0]
    normalised_weights = weights/sum(weights)
    solution = normalised_weights*sum(genome)
    
    #print(W1)
    #convert the newExposure vector into list type structure
    newExposure = list(solution)
    
    
    # get the maximum value of the new Exposure
    maxcoef = max(newExposure)
    idxmaxcoef = newExposure.index(maxcoef)
    
    newExposure = np.round(newExposure)
    
    # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
    if np.sum(newExposure)!=maxmutation:
        newExposure[idxmaxcoef] = round(newExposure[idxmaxcoef])+maxmutation-sum(newExposure)
     
    # compute the estimated genome
    est_genome = np.dot(W, newExposure)
    
    if metric=="cosine":
        newSimilarity = cos_sim(genome, est_genome)
    else:
        newSimilarity = np.linalg.norm(genome-est_genome , ord=2)/np.linalg.norm(genome, ord=2)
    
    return (newExposure, newSimilarity)

# to do the calculation, we need: W, samples(one column), H for the specific sample, tolerence
# Fit signatures with multiprocessing
def fit_signatures_pool(total_genome, W, index):
    genome = total_genome[:,index]
    genome = np.array(genome)
    maxmutation = round(np.sum(genome))
     
    
    
    #initialize the guess
    x0 = np.random.rand(W.shape[1], 1)*maxmutation
    x0= x0/np.sum(x0)*maxmutation
    
    #the optimization step
    #set the bounds and constraints
    bnds = create_bounds([], genome, W.shape[1]) 
    cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genome]}
    
    
    #sol = scipy.optimize.minimize(parameterized_objective2_custom, x0, args=(W, genome),  bounds=bnds, constraints =cons1, tol=1e-30)
    #solution = sol.x
    
    
    ### using NNLS algorithm 
    reg = nnls(W,genome)
    weights = reg[0]
    normalised_weights = weights/sum(weights)
    solution = normalised_weights*sum(genome)
    
    #print(W1)
    #convert the newExposure vector into list type structure
    newExposure = list(solution)
    
    
    # get the maximum value of the new Exposure
    maxcoef = max(newExposure)
    idxmaxcoef = newExposure.index(maxcoef)
    
    newExposure = np.round(newExposure)
    
    # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
    if np.sum(newExposure)!=maxmutation:
        newExposure[idxmaxcoef] = round(newExposure[idxmaxcoef])+maxmutation-sum(newExposure)
     
    # compute the estimated genome
    est_genome = np.dot(W, newExposure)
    
    newSimilarity = cos_sim(genome, est_genome)
    
    return (newExposure, newSimilarity)


def add_signatures(W, genome, cutoff=0.05, presentSignatures=[], toBeAdded="all", metric="l2", solver="nnls", check_rule_negatives=[], check_rule_penalty=1.0, verbose = False):
     # This function takes an array of signature and a single genome as input, returns a dictionray of cosine similarity, exposures and presence 
     # of signatures according to the indices of the original signature array
    
    if toBeAdded == "all":
        notToBeAdded=[]
    else:               
        #print(len(list(range(W.shape[1]))))
        #print(len(toBeAdded))
        notToBeAdded = list(set(list(range(W.shape[1])))-set(toBeAdded)) # get the indices of the signatures to be excluded
        
        #print(len(notToBeAdded))
    originalSimilarity = 100 # it can be also written as oldsimilarity. wanted to put a big number
    maxmutation = round(np.sum(genome))
    init_listed_idx = presentSignatures 
    
    init_nonlisted_idx = list(set(list(range(W.shape[1])))-set(init_listed_idx)-set(notToBeAdded))
    
    finalRecord = [["similarity place-holder" ], ["newExposure place-holder"], ["signatures place-holder"]] #for recording the cosine difference, similarity, the new exposure and the index of the best signauture
    
    # get the initial original similarity if some signatures are already given to be present
    if len(init_listed_idx)!=0:
        
        loop_liststed_idx=init_listed_idx
        loop_liststed_idx.sort()
        
        W1 = W[:,loop_liststed_idx]
       
        #print (W1.shape)
        #initialize the guess
        x0 = np.random.rand(W1.shape[1], 1)*maxmutation
        x0= x0/np.sum(x0)*maxmutation
        
        #set the bounds and constraints
        bnds = create_bounds([], genome, W1.shape[1]) 
        cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genome]} 
        
        #the optimization step
        if solver == "slsqp":
            sol = minimize(parameterized_objective2_custom, x0, args=(W1, genome),  bounds=bnds, constraints =cons1, tol=1e-30)
            newExposure = list(sol.x)
        
        if solver == "nnls":
        ### using NNLS algorithm        
            reg = nnls(W1, genome[:,0])
            weights = reg[0]
            normalised_weights = weights/sum(weights)
            solution = normalised_weights*sum(genome)
            newExposure = list(solution)
        #print(W1)
        #convert the newExposure vector into list type structure
        
        #newExposure = list(solution) #for nnls
                           
                          
        # get the maximum value of the new Exposure
        maxcoef = max(newExposure)
        idxmaxcoef = newExposure.index(maxcoef)
        
        newExposure = np.round(newExposure)
        
        # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
        if np.sum(newExposure)!=maxmutation:
            newExposure[idxmaxcoef] = round(newExposure[idxmaxcoef])+maxmutation-sum(newExposure)
         
        # compute the estimated genome
        est_genome = np.dot(W1, newExposure)
        if metric=="cosine":
            originalSimilarity = 1-cos_sim(genome[:,0], est_genome )
        elif metric=="l2":    
            originalSimilarity = np.linalg.norm(genome[:,0]-est_genome , ord=2)/np.linalg.norm(genome[:,0], ord=2)
        
            
        finalRecord = [originalSimilarity, newExposure, init_listed_idx]
        
        
    while True:
        
        bestDifference = -100  #### wanted to put a big initila number
        bestSimilarity = 100 #### actually bestdistance
        loopRecord = [["newExposure place-holder"], ["signatures place-holder"], ["best loop signature place-holder"]]
        
        for sig in init_nonlisted_idx:
            
            
            
            if len(init_listed_idx)!=0:
                loop_liststed_idx=init_listed_idx+[sig]
                loop_liststed_idx.sort()
                #print(loop_liststed_idx)
                W1 = W[:,loop_liststed_idx]
                #print (W1.shape)
                #initialize the guess
                x0 = np.random.rand(W1.shape[1], 1)*maxmutation
                x0= x0/np.sum(x0)*maxmutation
                
                #set the bounds and constraints
                bnds = create_bounds([], genome, W1.shape[1]) 
                cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genome]} 
            # for the first time addintion  
            else:
                W1 = W[:,sig][:,np.newaxis]
                #print (W1.shape)        
                #initialize the guess
                x0 = np.ones((1,1))*maxmutation    
            
                #set the bounds and constraints
                bnds = create_bounds([], genome, 1) 
                cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genome]} 
            
            if solver == "slsqp":
                #the optimization step
                sol = minimize(parameterized_objective2_custom, x0, args=(W1, genome),  bounds=bnds, constraints =cons1, tol=1e-30)
                newExposure = list(sol.x)
            if solver == "nnls":
                reg = nnls(W1, genome[:,0])
                weights = reg[0]
                normalised_weights = weights/sum(weights)
                solution = normalised_weights*sum(genome)
                newExposure = list(solution) #for nnls
            #print(W1)
            #convert the newExposure vector into list type structure
            
            
            
            # get the maximum value of the new Exposure
            maxcoef = max(newExposure)
            idxmaxcoef = newExposure.index(maxcoef)
            
            newExposure = np.round(newExposure)
            
            # We may need to tweak the maximum value of the new exposure to keep the total number of mutation equal to the original mutations in a genome
            if np.sum(newExposure)!=maxmutation:
                newExposure[idxmaxcoef] = round(newExposure[idxmaxcoef])+maxmutation-sum(newExposure)
             
            # compute the estimated genome
            est_genome = np.dot(W1, newExposure)
            
            if metric=='cosine':
                newSimilarity = 1-cos_sim(genome[:,0], est_genome)
            elif metric=='l2':
                newSimilarity = np.linalg.norm(genome[:,0]-est_genome , ord=2)/np.linalg.norm(genome[:,0], ord=2)
            
            if sig in check_rule_negatives:
                
                newSimilarity = newSimilarity*check_rule_penalty
                
            difference = originalSimilarity -  newSimilarity 
            
            # record the best values so far that creates the best difference
           
            if difference>bestDifference:
                bestDifference = difference
                bestSimilarity = newSimilarity
                loopRecord = [newExposure, W1, sig]  #recording the cosine difference, the new exposure and the index of the best signauture
                
        
         
         
        if originalSimilarity - bestSimilarity>cutoff: 
            
            
            originalSimilarity = bestSimilarity
            init_listed_idx.append(loopRecord[2])  #add the qualified signature in the signature list
            init_nonlisted_idx.remove(loopRecord[2]) #remover the qualified signature from the to be added list
            init_listed_idx.sort()
            #print(originalSimilarity)
            finalRecord = [originalSimilarity, loopRecord[0], init_listed_idx, loopRecord[1], genome]
            if verbose==True:
                print("Current signaure indices: ", init_listed_idx)
                if metric == "cosine":
                    print("Cosine similarity: ", 1-originalSimilarity)
                elif metric == "l2":
                    print("L2 norm: ", originalSimilarity)
                print("\n")
            if len(init_nonlisted_idx)!= 0:
                
                continue
            else:
                break
        else:
            
            break
     
    #print("Good so far")
    #print(finalRecord)
    dictExposure= {"similarity":finalRecord[0], "exposures":finalRecord[1], "signatures": finalRecord[2]}  
    addExposure = np.zeros([W.shape[1]])
    addExposure[dictExposure["signatures"]]=dictExposure["exposures"]
    #print(cos_sim(genome[:,0], np.dot(W,addExposure)))
    if metric == "cosine":
        return  addExposure, cos_sim(genome[:,0], np.dot(W,addExposure)), cos_sim(genome[:,0], np.dot(W,addExposure)) #first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity
    elif metric == "l2":
        return  addExposure, np.linalg.norm(genome[:,0]-np.dot(W, addExposure) , ord=2)/np.linalg.norm(genome[:,0], ord=2), cos_sim(genome[:,0], np.dot(W,addExposure)) #first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity







def remove_all_single_signatures(W, H, genomes, metric="ls", solver = "nnls", cutoff=0.05, background_sigs = [], verbose=False):
    # make the empty list of the successfull combinations
    successList = [0,[],0] 
    background_sig = copy.deepcopy(background_sigs)
    if metric == "cosine":
        # get the cos_similarity with sample for the oringinal W and H[:,i]
        originalSimilarity= 1-cos_sim(genomes, np.dot(W, H))
        if verbose==True:
            print("originalSimilarity", 1-originalSimilarity)
    elif metric == "l2":
       
        originalSimilarity = np.linalg.norm(genomes-np.dot(W, H) , ord=2)/np.linalg.norm(genomes, ord=2)
        if verbose==True:
            print("originalSimilarity", originalSimilarity)
    # make the original exposures of specific sample round
    oldExposures = np.round(H)
    
    # set the flag for the while loop
    if len(oldExposures[np.nonzero(oldExposures)])>1:
        Flag = True
    else: 
        Flag = False
        if metric=="cosine":
            return oldExposures, 1-originalSimilarity, cos_sim(genomes, np.dot(W,H))  #first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity
        else:
            return oldExposures, originalSimilarity, cos_sim(genomes, np.dot(W,H))    #first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity
    # The while loop starts here
    while Flag: 
        
        # get the list of the indices those already have zero values
        if len(successList[1]) == 0:
            initialZerosIdx = list(np.where(oldExposures==0)[0]) 
            #get the indices to be selected
            selectableIdx = list(np.where(oldExposures>0)[0]) 
        elif len(successList[1]) > 1: 
            initialZerosIdx = list(np.where(successList[1]==0)[0]) 
            #get the indices to be selected
            selectableIdx = list(np.where(successList[1]>0)[0])
        else:
            print("iteration is completed")
            #break
        
        
        # get the total mutations for the given sample
        maxmutation = round(np.sum(genomes))
        
        # new signature matrix omiting the column for zero
        #Winit = np.delete(W, initialZerosIdx, 1)
        Winit = W[:,selectableIdx]
        
        # set the initial cos_similarity or other similarity distance
        record  = [100, [], 1]   #
        # get the number of current nonzeros
        l= Winit.shape[1]
        
        background_sig = get_changed_background_sig_idx(list(oldExposures), background_sig)
        
        for i in range(l):
            
            if i in background_sig:
                continue
            
            loopSelection = list(range(l))
            del loopSelection[i]
            #print (loopSelection)
            W1 = Winit[:,loopSelection]
           
            
            #initialize the guess
            x0 = np.random.rand(l-1, 1)*maxmutation
            x0= x0/np.sum(x0)*maxmutation
            
            #set the bounds and constraints
            bnds = create_bounds([], genomes, W1.shape[1]) 
            cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genomes]} 
            
            #the optimization step
            if solver == "slsqp":
                sol = minimize(parameterized_objective2_custom, x0, args=(W1, genomes),  bounds=bnds, constraints =cons1, tol=1e-15)
                
                #print (sol.success)
                #print (sol.x)
                
                #convert the newExposure vector into list type structure
                newExposure = list(sol.x)
            if solver == "nnls":
                ### using NNLS algorithm 
                reg = nnls(W1,genomes)
                weights = reg[0]
                normalised_weights = weights/sum(weights)
                solution = normalised_weights*sum(genomes)                
                newExposure = list(solution)
                
            #insert the loopZeros in its actual position 
            newExposure.insert(i, 0)
            
            #insert zeros in the required position the newExposure matrix
            initialZerosIdx.sort()
            for zeros in initialZerosIdx:
                newExposure.insert(zeros, 0)
            
            # get the maximum value the new Exposure
            maxcoef = max(newExposure)
            idxmaxcoef = newExposure.index(maxcoef)
            
            newExposure = np.round(newExposure)
            
            
            if np.sum(newExposure)!=maxmutation:
                newExposure[idxmaxcoef] = round(newExposure[idxmaxcoef])+maxmutation-sum(newExposure)
                
            
            newSample = np.dot(W, newExposure)
            if verbose==True:
                print(newExposure)
            
            if metric == "cosine":
                newSimilarity = 1-cos_sim(genomes, newSample) 
                if verbose==True:
                    print("newSimilarity",1-newSimilarity) 
            elif metric == "l2":
                newSimilarity = np.linalg.norm(genomes-newSample, ord=2)/np.linalg.norm(genomes, ord=2)
                if verbose==True:
                    print("newSimilarity",newSimilarity) 
            difference =  newSimilarity -originalSimilarity 
            
            if verbose==True:
                print("difference", difference) 
            if difference<record[0]:
                record = [difference, newExposure, newSimilarity]
                
                #print("difference", difference)
            if verbose==True:
                print("--------------------------------------")
        if verbose==True:    
            print("\n") 
            print("Selected Exposure")
            print(record[1])
            print("New Similarity", record[2])
            print("****************************")
            #print ("This loop's selection is {}".format(record))
        
        if record[0]>cutoff:   
            Flag=False
        elif len(record[1][np.nonzero(record[1])])==1:
            successList = record 
            Flag=False
        else:
            successList = record
            background_sig = get_changed_background_sig_idx(list(record[1]), background_sig)
        #print("The loop selection is {}".format(successList))
        
        #print (Flag)
        #print ("\n\n")
    
    #print ("The final selection is {}".format(successList))
    
    if len(successList[1])==0:
        successList = [0.0, oldExposures, originalSimilarity]
    
    if verbose==True:
        print("\n") 
        print("Final Exposure")
        print(successList[1])
        if metric=="cosine":
            print("Final Similarity: ", 1-successList[2])
        if metric=="l2":
            print("Final Similarity: ", successList[2])
    H = successList[1]
    return H, successList[2], cos_sim(genomes, np.dot(W,H)) #first value is the exprosure, second value is the similarity (based on the similarity matic), third value is the cosine similarity
    


def remove_all_single_signatures_pool(indices, W, exposures, totoalgenomes):
    i = indices
    H = exposures[:,i]
    genomes= totoalgenomes[:,i]
    # make the empty list of the successfull combinations
    successList = [0,[],0] 
    # get the cos_similarity with sample for the oringinal W and H[:,i]
    originalSimilarity= cos_sim(genomes, np.dot(W, H))
    # make the original exposures of specific sample round
    oldExposures = np.round(H)
    
    # set the flag for the while loop
    if len(oldExposures[np.nonzero(oldExposures)])>1:
        Flag = True
    else: 
        Flag = False
        return oldExposures
    # The while loop starts here
    while Flag: 
        
        # get the list of the indices those already have zero values
        if len(successList[1]) == 0:
            initialZerosIdx = list(np.where(oldExposures==0)[0]) 
            #get the indices to be selected
            selectableIdx = list(np.where(oldExposures>0)[0]) 
        elif len(successList[1]) > 1: 
            initialZerosIdx = list(np.where(successList[1]==0)[0]) 
            #get the indices to be selected
            selectableIdx = list(np.where(successList[1]>0)[0])
        else:
            print("iteration is completed")
            #break
        
        
        # get the total mutations for the given sample
        maxmutation = round(np.sum(genomes))
        
        # new signature matrix omiting the column for zero
        #Winit = np.delete(W, initialZerosIdx, 1)
        Winit = W[:,selectableIdx]
        
        # set the initial cos_similarity
        record  = [0.11, []] 
        # get the number of current nonzeros
        l= Winit.shape[1]
        
        for i in range(l):
            #print(i)
            loopSelection = list(range(l))
            del loopSelection[i]
            #print (loopSelection)
            W1 = Winit[:,loopSelection]
           
            
            #initialize the guess
            x0 = np.random.rand(l-1, 1)*maxmutation
            x0= x0/np.sum(x0)*maxmutation
            
            #set the bounds and constraints
            bnds = create_bounds([], genomes, W1.shape[1]) 
            cons1 ={'type': 'eq', 'fun': constraints1, 'args':[genomes]} 
            
            #the optimization step
            sol = minimize(parameterized_objective2_custom, x0, args=(W1, genomes),  bounds=bnds, constraints =cons1, tol=1e-15)
            
            #print (sol.success)
            #print (sol.x)
            
            #convert the newExposure vector into list type structure
            newExposure = list(sol.x)
            
            #insert the loopZeros in its actual position 
            newExposure.insert(i, 0)
            
            #insert zeros in the required position the newExposure matrix
            initialZerosIdx.sort()
            for zeros in initialZerosIdx:
                newExposure.insert(zeros, 0)
            
            # get the maximum value the new Exposure
            maxcoef = max(newExposure)
            idxmaxcoef = newExposure.index(maxcoef)
            
            newExposure = np.round(newExposure)
            
            
            if np.sum(newExposure)!=maxmutation:
                newExposure[idxmaxcoef] = round(newExposure[idxmaxcoef])+maxmutation-sum(newExposure)
                
            
            newSample = np.dot(W, newExposure)
            newSimilarity = cos_sim(genomes, newSample) 
             
            difference = originalSimilarity - newSimilarity
            #print(originalSimilarity)
            #print(newSample)
            #print(newExposure)
            #print(newSimilarity)
            
            #print(difference)
            #print (newExposure)
            #print (np.round(H))
            #print ("\n\n")
             
            if difference<record[0]:
                record = [difference, newExposure, newSimilarity]
            
            
        #print ("This loop's selection is {}".format(record))
        
        if record[0]>0.01:   
            Flag=False
        elif len(record[1][np.nonzero(record[1])])==1:
            successList = record 
            Flag=False
        else:
            successList = record
        #print("The loop selection is {}".format(successList))
        
        #print (Flag)
        #print ("\n\n")
    
    #print ("The final selection is {}".format(successList))
    
    if len(successList[1])==0:
        successList = [0.0, oldExposures, originalSimilarity]
    #print ("one sample completed")
    return successList[1]




def add_remove_signatures(W,sample, metric="l2", solver="nnls", background_sigs = [], permanent_sigs = [], candidate_sigs="all", penalty = 0.05, check_rule_negatives =[],checkrule_penalty = 1.00, verbose=False):
    
    
    always_background = copy.deepcopy(permanent_sigs)
    M = sample   
    if candidate_sigs == "all":
        candidate_sigs = list(range(W.shape[1]))
    #first add each signatures
    original_distance = 100
    layer = 0
    
    # check the cosine_similarity with 4 signatures (highest)
    cosine_similarity_with_four_signatures = 1.0 # a value that will allow the signature not be reported as novel signature
    while True:
        if verbose:
            print("\n\n\n\n!!!!!!!!!!!!!!!!!!STARTING LAYER: ", layer)
        layer=layer+1
        layer_original_distance = 100
        sigsToBeAdded = list(set(candidate_sigs)-set(background_sigs))
        
        for i in sigsToBeAdded:
            loop_sig = [i] 
            
            background_sigs = list(set(background_sigs).union(set(always_background)))
            #print(background_sigs)
            add_exposures, add_distance, _ = add_signatures(W, M[:,np.newaxis], presentSignatures=copy.deepcopy(background_sigs), toBeAdded=loop_sig, 
                                                  metric="l2", verbose=False, check_rule_negatives=check_rule_negatives, check_rule_penalty=1.00, cutoff=penalty) 
            if verbose:
                print("\n\n\n################## Add Index {} ########################".format(i)) 
                print(np.nonzero(add_exposures)[0])
                print(sigsToBeAdded)
                print(add_distance)
                
                
            remove_exposures, remove_distance, _ = remove_all_single_signatures(W, add_exposures, M, background_sigs = always_background, metric="l2", verbose = False, cutoff=0.01)
            if verbose:
                print("\n################## Remove ########################")
                print(np.nonzero(remove_exposures)[0])
                print(remove_distance)
            # check if there is any change between the add and remove signatures and assign the distance and exposure accoridingly    
            if (np.nonzero(add_exposures)[0].all()==np.nonzero(remove_exposures)[0].all()) and np.nonzero(add_exposures)[0].shape == np.nonzero(remove_exposures)[0].shape:
                distance = add_distance
                exposures = add_exposures
            else:
                distance = remove_distance
                exposures = remove_exposures
                
            if distance < layer_original_distance:
                selected_signatures = np.nonzero(exposures)[0]
                layer_original_distance = distance 
                activities = exposures
        if verbose:
            print("\n\n#################### After Add-Remove #################################")
            print(selected_signatures)
            print(layer_original_distance)
        if layer_original_distance < original_distance:
            original_distance = layer_original_distance
            background_sigs = list(selected_signatures)
            finalactivities = activities
            if len(background_sigs)==4:
                cosine_similarity_with_four_signatures = cos_sim(M, np.dot(W,finalactivities))
                #print(cosine_similarity_with_four_signatures)
        else:
            cosine_similarity = cos_sim(M, np.dot(W,finalactivities))
            break
        
    if verbose:
        print("\n########################## Final ###########################")        
        print(background_sigs)
        print(original_distance)
        print(finalactivities)
    
    #newExposure, newSimilarity = fit_signatures(W[:,list(background_sigs)], M, metric="l2")
    #print(newExposure, newSimilarity)
    return (background_sigs, finalactivities, original_distance, cosine_similarity, cosine_similarity_with_four_signatures)