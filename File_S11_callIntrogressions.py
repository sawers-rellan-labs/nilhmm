# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 13:05:58 2024

@author: jholland
"""

#Function to implement HMM to call introgressions in NILs
import numpy as np
from hmmlearn import hmm

######################
#make a function to fit one NIL genotype sequence from one chromosome at a time
#return a numpy array with 1 dimenions representing most likely states of one line
def call_intros_one_chrom(i, hmm_model, geno_mat):
    nil_i  = np.identity(4)[np.nan_to_num(geno_mat[i,:]).astype(int)]
    preds = hmm_model.predict(nil_i).reshape(1,-1)
    return(preds)

######################

################
#Make a function to build HMM model, call introgressions,
#and summarize percent introgression het and homozygous calls
#using specific values of nir and scer
################

def call_intros(geno, marker_dict, nir, germ, gert, p, mr, r, f_1, f_2, return_calls = True):

    """
    for the HMM we need starting probabilities, a vector of hidden states, and a vector of observed states 
    a transition probability matrix and an emmision probability matrix

    ARGUMENTS:
        geno = numpy array of numeric genotype calls, individuals in rows, markers in columns, 0,1,2 are minor allele counts, 3 is missing data
        marker_dict = dictionary of lists of markers as components of a dictionary
        nir = non-informative rate
        germ = SNP calling error rate on true introgression homozygotes
        gert = SNP calling error rate on true introgression heterozygotes
        p = proportion of homozygous SNP call errors resulting in a het call
        mr = missing call rate
        r = expected recombination rate between adjacent markers
        f_1 = expected frequency of heterozyzotes
        f_2 = expected frequency of homozygous introgressions
        return_calls = T/F flag whether to return Viterbi best calls

    """
    #starting probabilities:

    #f(B73 homoz) = 1 - others
    f_0 = 1 - f_1 - f_2

    states = [0, 1, 2] #"B73", "het",  "donor"

    #prob. transition from 0 to 0 = no recombination
    p00 = 1 - r
    
    #prob. transition from 0 to 1 | given that first gene is 0 class = prob(recomb)*(prob het at 2nd locus)
    # = prob(recomb)*((freq het)/(freq het + freq homozygous alt allele at 2nd locus))
    p01 = r*(f_1/(f_1 + f_2))
    
    #prob. transition from 0 to 2 = prob(recomb)*(prob homoz alt at 2nd locus)
    # = prob(recomb)*(1 - ((freq het)/(freq het + freq homozygous alt allele at 2nd locus))
    p02 = r*(f_2/(f_1 + f_2))
    
    #prob transition from 1 to 0 = prob(recom) * f_0/(f0+f2)
    #if there was a recomb between double hets, we get recurrent homozygote with proportion of f_0 out of total of homozygous classes
    p10 = r*f_0/(f_0 + f_2)
    
    #prob. transition from 1 to 1 = no recombination
    p11 = 1 - r
    
    #if there was a recomb between double hets, we get donor homozygote with proportion of f_2 out of total of homozygous classes
    p12 = r*f_2/(f_0 + f_2)
    
    #prob transition from 2 to 0 = freq recombination * freq recurrent homoz among recurrent homoz + het classes 
    p20 = r*f_0/(f_0 + f_1)

    #prob transition from 2 to 1 = freq recombination * freq het among recurrent homoz + het classes 
    p21 = r*f_1/(f_0 + f_1)
    
    #prob. transition from 2 to 2 = no recombination
    p22 = 1 - r
    
    tmat = np.array([[p00, p01, p02],
                    [p10, p11, p12],
                    [p20, p21, p22]]) #ignoring double recombinations in these small regions
    
    #transition matrix rows must sum to 1 across columns
    # print("transmission matrix:")
    # print(tmat)
    # print("transmission matrix sum across columns:")
    # print(np.sum(tmat, axis = 1))

    #Emission matrix is probability of observing a value of 0,1,2, or 3 (NA) given each true state
    #rows are true states, columns are observed states, so it's a 3 x 4 matrix in this case
    emimat = np.array([[(1-germ)*(1-mr), p*germ*(1-mr), (1-p)*germ*(1-mr), mr],
                       [(((1-nir)*0.5*gert) + nir*(1-germ))*(1-mr), (((1-nir)*(1-gert)) + (nir*germ*p))*(1-mr),(((1-nir)*0.5*gert) + nir*germ*(1-p))*(1-mr), mr],
                       [((1-nir)*germ*(1-p) + (nir*(1-germ)))*(1-mr), germ*p*(1-mr), ((1-nir)*(1-germ) + (nir*germ*(1-p)))*(1-mr), mr]])
                       
    # print("emission matrix:")
    # print(emimat)
    # print("emission matrix sum across columns:")
    # print(np.sum(emimat, axis = 1))
    
    #set up the HMM
    model = hmm.MultinomialHMM(
        n_components=len(states),
        n_trials=1,
        init_params='')
    
    model.n_features = len(states) + 1 #we have an additional NA 'feature'

    model.startprob_ = np.array([f_0, f_1, f_2]) #note that this is correct for NILs, but not parents
    model.transmat_ = tmat
    model.emissionprob_ = emimat
    
    results = {} #make a dictionary, each element will be results from one chromosome
    for chr in [1,2,3,4,5,6,7,8,9,10]:
        results[chr] = [] #make empty list for current chromosome
        geno_current_chr = geno[:,marker_dict[chr]] #subset genos for current chrom
    
        #print("start processing markers on chrom " + str(chr))
        #start = time.process_time()
        for i in range(geno.shape[0]):
            results[chr].append(call_intros_one_chrom(i, hmm_model = model, geno_mat = geno_current_chr))
        
        #print("finish processing markers on chrom " + str(chr))
        #print("time elapsed")
        #print(time.process_time() - start)
    
    #now have to pack results back into a big array
    resultsByChrom = {}
    
    for chr in [1,2,3,4,5,6,7,8,9,10]:
        resultsByChrom[chr] = np.concatenate(results[chr], axis = 0)
        
    #Now we have a dict of arrays, one per chromosome. Last step is to merge these into one array
    NIL_calls = np.concatenate(list(resultsByChrom.values()), axis = 1)
    
    
    if return_calls:
        return(NIL_calls)


####################################
