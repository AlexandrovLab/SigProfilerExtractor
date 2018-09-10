#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 11:04:39 2018

@author: mishugeb
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 13:39:29 2018

@author: mishugeb
"""

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


def nnmf(genomes, nfactors):
    genomes = np.array(genomes)
    nmf = nimfa.Nmf(genomes, max_iter=100000, rank=nfactors, update='divergence', objective='div', test_conv= 30)
    nmf_fit = nmf()
    W = nmf_fit.basis()
    H = nmf_fit.coef()
    return W, H

def BootstrapCancerGenomes(genomes):
    np.random.seed() # Every time initiate a random seed so that all processors don't get the same seed
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




def nmf(pool_constant=1, genomes=1, totalProcesses=1):
    
    genomes = pd.DataFrame(genomes) #creating/loading a dataframe/matrix
    bootstrapGenomes= BootstrapCancerGenomes(genomes)
    bootstrapGenomes[bootstrapGenomes<0.0001]= 0.0001
    W, H = nnmf(bootstrapGenomes,totalProcesses)  #uses custom function nnmf
    #print ("initital W: ", W); print("\n");
    #print ("initial H: ", H); print("\n");
    W = np.array(W)
    H = np.array(H)
    total = W.sum(axis=0)[np.newaxis]
    #print ("total: ", total); print("\n");
    W = W/total
    H = H*total.T
    print ("process " +str(totalProcesses)+" continues please wait")
    return W, H


def parallel_runs(genomes=1, totalProcesses=1, iterations=1,  n_cpu=-1, verbose = False):
    if verbose:
        print ("Precess "+str(totalProcesses)+ " is in progess")
    if n_cpu==-1:
        pool = multiprocessing.Pool()
    else:
        pool = multiprocessing.Pool(processes=n_cpu)
        
  
    pool_nmf=partial(nmf, genomes=genomes, totalProcesses=totalProcesses)
    result_list = pool.map(pool_nmf, range(iterations)) 
    pool.close()
    pool.join()
    return result_list




def extract(genomes=1, totalProcesses=1, totalIterations=1, verbose=False, n_cpu=-1):
    

    totalMutationTypes = genomes.shape[0];
    totalGenomes = genomes.shape[1];
    Wall = np.zeros((totalMutationTypes, totalProcesses * totalIterations));
    Hall = np.zeros((totalProcesses * totalIterations, totalGenomes));
    genomeErrors = np.zeros((totalMutationTypes, totalGenomes, totalIterations));
    genomesReconstructed = np.zeros((totalMutationTypes, totalGenomes, totalIterations))
    
    
    results=parallel_runs(genomes=genomes, totalProcesses=totalProcesses, iterations=totalIterations, n_cpu=n_cpu, verbose=verbose)   
    
    processCount=0
    for i in range(len(results)):
        W = results[i][0]
        H = results[i][1]
        genomeErrors[:, :, i] = genomes -  np.dot(W,H);
        genomesReconstructed[:, :, i] = np.dot(W,H);
        Wall[ :, processCount : (processCount + totalProcesses) ] = W;
        Hall[ processCount : (processCount + totalProcesses), : ] = H;
        processCount = processCount + totalProcesses;
    return genomeErrors, genomesReconstructed, Wall, Hall


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
	dot_product = np.dot(a, b)
	norm_a = np.linalg.norm(a)
	norm_b = np.linalg.norm(b)
	return dot_product / (norm_a * norm_b)





################################################################### FUNCTION ONE ###################################################################
# function to calculate the centroids

def compute_pairwise_centroid(vec1, vec2):
    """Takes 2 vectors vec1, vec2 and returns the centroids of 
    vec1 and vec2 as a vector named "centroids". Here, the "centroid"
    vector is a list type variable. For better understanding please the run 
    code and check the variables in the "test of function/example section below 
    the function"
    
    Dependencies:
        *Does not require any custom function (constructed by me)
        
    Required by:
        *pairwise_cluster_init 
        *pairwise_cluster_elong
        
	"""
    
    centroids= []  
    for i in range(0, len(vec1)):
        centroid= (vec1[i]+vec2[i])/2
        centroids.append(centroid)
        
    return centroids


################################################################### FUNCTION  ###################################################################
def pairwise_cluster_raw(mat1=([0]), mat2=([0]), mat1T=([0]), mat2T=([0])):  # the matrices (mat1 and mat2) are used to calculate the clusters and the lsts will be used to store the members of clusters
    
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
    con_mat = np.zeros((mat1.shape[1],mat2.shape[1]))
    for i in range(0, mat1.shape[1]):
        for j in range(0,mat2.shape[1]):
            con_mat[i, j] = cos_sim(mat1[:,i], mat2[:,j])  #used custom function
    #debug
    #print (con_mat); print("\n")
        
    lstCluster=[]
    lstClusterT=[]
    idxPair= []  
    
    
    for a in range(0,mat1.shape[1]):
    
        i,j=np.unravel_index(con_mat.argmax(), con_mat.shape) 
        lstCluster.append([mat1[:,i], mat2[:,j]])
        idxPair.append([i,j])   #for debugging
        lstClusterT.append([mat1T[i,:], mat2T[j, :]])
        con_mat[i,:]=0
        con_mat[:,j]=0
        

        
    
    return lstCluster, idxPair, lstClusterT
         





################################################################### FUNCTION  ###################################################################
def pairwise_cluster_init(matPair, matPairT):
    
    """Takes a list of a pair of matrices, "matPair". It run the "pairwise_cluster_raw" inside
    it and returns a list that contains the list of cluster and the list of the corresponding 
    centroids of the each clusters. 
    
    Dependencies:
        *compute_pairwise_centroid
        *pairwise_cluster_raw
        
    Required by:
        *cluster_init
        
    """
    
    mat1 = matPair[0]; mat2 = matPair[1]
    mat1T = matPairT[0]; mat2T = matPairT[1]
    
    lstCluster,idxPairs, lstClusterT= pairwise_cluster_raw(mat1, mat2, mat1T, mat2T)   #used custom function "

    
    lstCen=[]
    for i in lstCluster: 
                
        lstCen.append(compute_pairwise_centroid(i[0],i[1]))
    return [lstCluster, np.array(lstCen), lstClusterT]





################################################################### FUNCTION  ###################################################################    
def cluster_init(listOfArrays, listOfArraysT):
    
    """Takes an array of list of matrix pairs "listOfArrays" (see the example below the function). 
    It returns a list of list of clusters-centroids pairs which are used for the next iteration 
    where "cluster_elong" is used . 
    
    Dependencies:
        *pairwise_cluster_init
    
    Required by: 
        *find_clusters
    """
    
    fulllist=[] 
    for i,j in zip(listOfArrays, listOfArraysT):
        #i an j contains have the similar structure as matPair and matPairT, respectively. see "pairwise_cluster_init".
        
        if len(i)>1:
             newlst=pairwise_cluster_init(i, j)  #custom function 
             
            
        else: 
            oddMat=[]#store the information the odd member/last matrix
            oddcen=[]
            oddMatT =[]
            for k in range(0,i[0].shape[1]):
                oddMat.append([i[0][:,k]])
                oddcen.append(i[0][:,k])
                oddMatT.append([j[0][k,:]])
            newlst = [oddMat, list(oddcen), oddMatT] 
            #print ("The odd member cannot be clustered in this round")
        
        fulllist.append(newlst)
    return fulllist


  
  
################################################################### FUNCTION  ###################################################################
    
def pairwise_cluster_elong(clust1, clust2, cen1, cen2, clust1T, clust2T):
    
    """ Takes a pair of clusters clust1, clust2 and the curresponding 
    centroids cen1, cen2. It merges the list of clusters to make a bigger 
    clusters with a double member in each cluster and computes the conbined 
    centroids for the correstponding clusters. It retrunts a list which contains 
    the list of the clusters "clusterNew" and the corresponding list of their
    centroids "cenNew"
    
     Dependencies:
        *compute_pairwise_centroid
        *pairwise_cluster_raw
        
    Required by:
        *cluster_elong
    """
    
    mat1=np.zeros([len(cen1[0]), len(cen1)])  # get the matrix size equal to datapoints or mutaion types / number of processes or ranks
    mat2=np.zeros([len(cen2[0]), len(cen2)])  # same as above for the matrix 2
    
    # the length of the list  in the first iteration will have one 
    

    
    
    for i in range(0, mat1.shape[1]):
        
        # the length of the list  in the first iteration will have one. We need to put use list1 
        
        mat1[:,i]= cen1[i] # populate the matrix with the calculated centroids
        mat2[:,i]= cen2[i]  # same as above for the matrix 2 
    
        
    lstCluster,idxPairs, lstClusterT = pairwise_cluster_raw(mat1, mat2, mat1.T, mat2.T)   #used custom function "pairwise_cluster_raw
    
    #the list of new cluster
    clusterNew = []
    cenNew = []
    clusterNewT = []
    for k in idxPairs:
        x=(k[0])
        y=(k[1])
        #print (x, y)
        newList = clust1[x]+clust2[y]
        newListT = clust1T[x] + list(clust2T[y])
        newClust = compute_pairwise_centroid(cen1[x], cen2[y])
        clusterNew.append(newList)
        cenNew.append(newClust)
        clusterNewT.append(newListT)
        
    return [clusterNew, cenNew, clusterNewT]
        
    



################################################################### FUNCTION  ###################################################################
def cluster_elong(listOfArrays):
    
    """Takes a list of list of clusters-centroids pairs "listOfArrays" and returns
    a list of list of clusters-centroids pairs (members in each cluster are double in 
    count compared to the previous clusters). The output is used for the next iteration 
    until there is only one cluster-centroid function in the list. 
    
    
    Dependencies:
        *pairwise_cluster_elong
    
    Required by: 
        *find_clusters
    """
    
    fulllist=[]
    
    
    for i in listOfArrays:
        
      
        

        
        if len(i)>1:
             clust1, clust2, cen1, cen2, clust1T, clust2T = i[0][0], i[1][0], i[0][1], i[1][1], i[0][2], i[1][2]
             newlst=pairwise_cluster_elong(clust1, clust2,  cen1, cen2, clust1T, clust2T)  #custom function 
             
             
            
        else: 
            newlst = [i[0][0],i[0][1], i[0][2]]
            
    
            
        
        fulllist.append(newlst)
        
    return fulllist


################################################################### FUNCTION  ###################################################################

#####################################################
#Contvert the list of clusters to list of dataframes
#No description given since the function is quite strait forward
#This function is required to execute the "find_clusters" function. 


def to_dataframes(clusters):
    listOfClusters=[]
    for i in clusters:
        listOfClusters.append(pd.DataFrame(i).T)
    return listOfClusters







################################################################### FUNCTION  ###################################################################
    
def find_clusters_v1(matList, matListT):
    
    """Takes a list of matrices "listOfMatrices" as argument. All of the matrices should have the 
    equal shapes. The function makes a partition based clustering (the number of clusters is equal 
    to the number of colums of the matrices, and not more column is assigned into a cluster from 
    a single matrix). It return the list of clusters  as "clusters". Please run the 
    test of function/example code provided below the for better understanding.  
    """
    
    W= split_list(matList, 2)
    H= split_list(matListT, 2)
    
    iteration="first"
    clusterSizeDifference=1  #Flag
    while clusterSizeDifference>0:
        #print ("processing")  ############################# for debugging 
        if  iteration=="first":
            listOfMatrices = cluster_init(W, H)
            listOfMatrices= split_list(listOfMatrices, 2)
            iteration = "next"
        else: 
            before = len(listOfMatrices[0][0][0][0])
            listOfMatrices = cluster_elong(listOfMatrices)
            listOfMatrices = split_list(listOfMatrices, 2)
            after = len(listOfMatrices[0][0][0][0])
            clusterSizeDifference = after-before
    clusters =  listOfMatrices[0][0][0]
    clustersT = listOfMatrices[0][0][2]
    
    
    
    listOfClusters = to_dataframes(clusters)
    listOfClustersT = to_dataframes(clustersT)
    
    




# =============================================================================
#    Calculation of Silhoutte coefficient    
# =============================================================================
    #from sklearn import metrics
    
    
    count = 0
    labels=[]
    clusters =[]
    for i in listOfClusters:
        
        clusters.append(i.T)
        for k in i:
            labels.append(count)
        count= count+1
    
    
    array = pd.concat(clusters, axis=0, ignore_index=True)
    
    try:
        SilhouetteCoefficients = metrics.silhouette_samples(array, labels, metric='cosine')
    
        
    
    except:
        SilhouetteCoefficients = np.ones((len(labels),1))
        
        
        
        
    avgSilhouetteCoefficients = np.mean(SilhouetteCoefficients)
    
    #clusterSilhouetteCoefficients 
    splitByCluster = np.array_split(SilhouetteCoefficients, len( listOfClusters))
    clusterSilhouetteCoefficients = np.array([])
    for i in splitByCluster:
        clusterSilhouetteCoefficients=np.append(clusterSilhouetteCoefficients, np.mean(i))
        
    
    
    return listOfClusters, listOfClustersT, avgSilhouetteCoefficients, clusterSilhouetteCoefficients
        


def analysis_signatures(genomes=[0], startprocesses =2, endprocesses=4, totalIterations= 10, verbose=False, n_cpu=-1 ):
    
    results = []
    
    for i in range(startprocesses , endprocesses+1, 1):
          
        processes=i
        
        print ("\nnumber of process is: ", processes )
        
        genomeErrors, genomesReconstructed, Wall, Hall = extract(genomes=genomes, totalProcesses=i , totalIterations=totalIterations, verbose=verbose, n_cpu=n_cpu)
       
       
            
        
        
        W= np.array_split(Wall, totalIterations, axis=1)
        H= np.array_split(Hall, totalIterations, axis=0)
        
              
        
        
        processclust, exposerclust, avgSilhouetteCoefficients, clusterSilhouetteCoefficients= find_clusters_v1(W, H)
        
       
        
        meanGenomeErrors = np.mean(genomeErrors, axis=2)
        meanGenomeReconstructed = np.mean(genomesReconstructed)    
        
        # computing the avg and std of the processes and exposures:
        
        processAvg = np.zeros((genomes.shape[0], processes))
        exposureAvg = np.zeros((processes, genomes.shape[1]))
        processStd = np.zeros((genomes.shape[0], processes))
        exposureStd = np.zeros((processes, genomes.shape[1]))
        
        for i in range(0, processes):
            processAvg[:,i]=np.mean(processclust[i], axis=1)
            processStd[:,i] = np.std(processclust[i], axis=1)
            exposureAvg[i,:] = np.mean(exposerclust[i], axis=1)
            exposureStd[i,:] = np.std(exposerclust[i], axis=1)
            
        results.append([genomes, processAvg, exposureAvg, processStd, exposureStd, avgSilhouetteCoefficients, clusterSilhouetteCoefficients, genomeErrors, genomesReconstructed, Wall, Hall, processes]) 
        
    return results


def cluster_similarity(mat1=([0]), mat2=([0])):  # the matrices (mat1 and mat2) are used to calculate the clusters and the lsts will be used to store the members of clusters
    
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
    con_mat = np.zeros((mat1.shape[1],mat2.shape[1]))
    for i in range(0, mat1.shape[1]):
        for j in range(0,mat2.shape[1]):
            con_mat[i, j] = cos_sim(mat1[:,i], mat2[:,j])  #used custom function
    #debug
    #print (con_mat); print("\n")
    similarity_mat   = pd.DataFrame(con_mat.copy())  
    lstCluster=[]
    idxPair= []  
    similarity = []
    
    for a in range(0,mat1.shape[1]):
        
        i,j=np.unravel_index(con_mat.argmax(), con_mat.shape) 
        similarity.append(con_mat.max())
        lstCluster.append([mat1[:,i], mat2[:,j]])
        idxPair.append([i,j])   #for debugging
        
        con_mat[i,:]=0
        con_mat[:,j]=0
        
     
        
    
    return  similarity_mat, idxPair, similarity, lstCluster





    
def main():        
    parser = argparse.ArgumentParser(description='Extraction of mutational signatures from Cancer genomes')
    
    parser.add_argument('project', type =str,
                        help= 'Name of the project file')
    
    parser.add_argument('refgen', type =str,
                        help= 'Name of the reference genome')
    
    parser.add_argument('minprocesses', type=int, 
                        help= 'The minimum number of processes to be extracted')
    
    parser.add_argument('maxprocesses', type=int, 
                        help= 'The maximum number of processes to be extracted')
    
    parser.add_argument('iterations', type=int, 
                        help= 'The number of iterations to be executed')
    
    
    parser.add_argument('--n_cpu', type=int, 
                        help= 'The number of cores to be used in the excecution.')
    
    
  
    
    parser.add_argument("--exome", help="Optional parameter instructs script to create the catalogues using only the exome regions. Whole genome context by default", action='store_true')
    parser.add_argument("--indel", help="Optional parameter instructs script to create the catalogue for limited INDELs", action='store_true')
    parser.add_argument("--extended_indel", help="Optional parameter instructs script to create the catalogue for extended INDELs", action='store_true')
    
    parser.add_argument('--mtypes', type =str,
                        help= 'The types of mutation. User should pass the inteded mutation types among to be analyzed\
                        separeted by coma "," with no space. The sigporfiler engine will analyze the specific mutation\
                        types those are passed to this argument. The valid mutation type are  96, 1536, 192, 3072 and  DINUC.\
                        For example, if the user wants analyze mutation type 96, 192 and DINUC, that person should pass\
                        "--mtypes 96,192,DINUC" as in the argument. If the argument is not used, all the mutations will\
                        be analyzed')
    
    args=parser.parse_args()
    
     
    
    
    data = datadump.sigProfilerMatrixGeneratorFunc (args.project, args.refgen, exome=False, indel=False, indel_extended=False)

   
    
    if args.n_cpu:
        n_cpu = args.n_cpu
    else:
        n_cpu = -1
    
     if args.exome:
        exome = True

    if args.extended_indel:
        indel = True

    if args.indel:
        indel = True
        limited_indel = True
    
    mlist = []
    for m in data:
        mlist.append(m)
    
    if mtypes:
        mtypes = args.mtypes.split(",")
        if any(x not in mlist for x in mtypes):
            print ("Please pass valid mutation types seperated by comma with no space. Also please use the uppercase characters")
            return None
    else:
        mtypes = mlist
   
    
    for m in mtypes:
        genomes = data[m]
        
        results = analysis_signatures(genomes=genomes, startprocesses = args.minprocesses, endprocesses=args.maxprocesses, totalIterations= args.iterations, n_cpu = n_cpu, verbose=True )
    
    
        f = open('../output/results_'+m, 'wb')
        pickle.dump(results, f)
        f.close()
    


if __name__ == "__main__":
    

    import argparse
    import scipy.io
    import numpy as np
    import pandas as pd
    import scipy.sparse as spr
    import nimfa
    import time
    from sklearn import metrics
    import pickle
    import multiprocessing
    from functools import partial
    import sigProfilerMatrixGeneratorFunc as datadump 
    main()