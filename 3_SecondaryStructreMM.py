# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 19:49:59 2022

@author: danie
"""
import csv
import numpy as np
import collections
from sys import exit
import sys
import json
from tqdm import tqdm
import pandas as pd
from PIL import Image
import copy
import pickle
import math
from os import listdir
from os.path import isfile, join
import os
import random

def main():
    # alignments = np.load('nw_Alignments.npy')
    nwMat = np.load('nwMattemp.npy')
    
    
    with open("alignments.txt", "rb") as fp:   # Unpickling
        alignment = pickle.load(fp)
        
    size = math.sqrt(np.size(nwMat))
    size = copy.deepcopy(int(size))
    currentDir = os.getcwd()
    
    fileNames = [f for f in listdir(currentDir) if isfile(join(currentDir, f))]
    
    
    # easier to make a dictionary where tuples map to aligment entries, than to play
    # with numpy 
    
    
    alignments = {}
    for i in range(int(size)):
        for j in range(int(size)):
            alignments[(i,j)] = alignment[j+size*i]
    
    PDBnames = []
    
    
    
    file1 = open('PDBnames.txt', 'r')
    for line in file1:
        PDBnames.append(line.strip())
    PDBnames.reverse()
    file1.close()
    
    HMMs = []
    counter = 0
    for i in range(len(PDBnames)):
        for j in range(len(PDBnames)):
            counter += 1
            HMMs.append([PDBnames[i], PDBnames[j], secondaryStructureReader(counter,PDBnames[0], PDBnames[i], PDBnames[j])])
    
    print(HMMs)
    print(len(HMMs))
    print(HMMs[0])
    exit()
    


def secondaryStructureReader(counter:int,proteinName:str, name1:str, name2:str):
    
    
    a = ['ala','arg','asn','asp','asx','cys','glu','gln','glx','gly','his','ile',
         'leu','lys','met','phe','pro','ser','thr','trp','tyr','val']
    b = 'ARNDBCEQZGHILKMFPSTWYV'
    
    
    proteinKeyForward = {}
    proteinKeyReverse = {}
    for q in range(len(a)):
        proteinKeyForward[a[q]] = b[q]
        proteinKeyReverse[b[q]] = a[q]
     
     
    secondary_1 = []
    secondary_2 = []
    inFile = open('%s.txt'%name1)
    for line in inFile:
        if line[0:3] == "LOC":
            secondary_1.append(line[0:30])
    inFile.close()
    
    inFile = open('%s.txt'%name2)
    for line in inFile:
        if line[0:3] == "LOC":
            secondary_2.append(line[0:24])
    inFile.close()
    secondary1 = {}
    secondary2 = {}
    secondary1
    
    for i in range(len(secondary_2)):
        secondary2[i] = secondary_2[i].split()[1].replace('\'','').replace('\"','')
    
    for i in range(len(secondary_1)):
        secondary1[i] = secondary_1[i].split()[1].replace('\'','').replace('\"','')
    
    second1 = []
    second2 = []
    for item in secondary_1:
        item.replace('\'','').replace('\"','')
        item1 = item.split()
        second1.append(item1[1:3])
    for item in secondary_2:
        item.replace('\'','').replace('\"','')
        item2 = item.split()
        second2.append(item2[1:3])
    
    
    
    #Hidden Markov Model Starts Here
    states1 = list((set([x.split()[1].replace('\'','').replace('\"','') for x in secondary_1])))
    states2 = list((set([x.split()[1].replace('\'','').replace('\"','') for x in secondary_2])))
    start_probability = returnStartProb(list(set(states1+states2)),secondary1,secondary2)
    # print(start_probability)
    # exit()
    transmission_probability = returnTransmissionProbability(list(set(states1+states2)),secondary1,secondary2)
    emission_probability = returnEmissionProbability(a, second1,second2)
    
    HMM1 = ''
    HMM2 = ''
    
    #random friendly start
    prev = 0
    rstart = {}
    for key,value in start_probability.items():
        rstart.update({key:round(value+prev,5)})
        prev = value+prev
    
    # random friendly transmission
    rtms1 = {}
    prev = 0
    for prob in transmission_probability:
        for key,value in prob.items():
            prev = 0
            for k,v in value.items():
                try:
                    rtms1[key]
                    rtms1[key].update({k:round(v+prev)})
                    prev += v
                except Exception as e:
                    # print(e)
                    rtms1[key] = {}
                    try:
                        rtms1[key].update({k:round(v+prev,5)})
                        prev+= v
                    except Exception as e2:
                        print(e2)
        break
        print('somehow it passed the break')
    rtms2 = {}
    for i in range(len(transmission_probability)):
        if i == 0:
            continue
        for key,value in transmission_probability[1].items():
            prev = 0
            for k,v in value.items():
                try:
                    rtms2[key]
                    rtms2[key].update({k:round(v+prev,5)})
                    prev+=v
                except Exception as e:
                    # print(e)
                    rtms2[key] = {}
                    try:
                        rtms2[key].update({k:round(v+prev,5)})
                        prev+=v
                    except Exception as e2:
                        print(e2)
        
    #random friendly emission
    
    e1 = emission_probability[0]
    re1 = {}
    e2 = emission_probability[1]
    re2 = {}
    prev = 0
    for key,value in e1.items():
        prev = 0
        for k,v in value.items():
            try:
                re1[key].update({k:round(v/sum(value.values()),5)})
                prev += v/sum(value.values())
            except:
                re1[key] = {}
                re1[key].update({k:round(v/(sum(value.values())),5)})
                prev+=v/sum(value.values())
    for key,value in e2.items():
        prev = 0
        for k,v in value.items():
            try:
                re2[key].update({k:round(v/sum(value.values()),5)})
                prev += v/sum(value.values())
            except:
                re2[key] = {}
                re2[key].update({k:round(v/(sum(value.values())),5)})
                prev+=v/sum(value.values())
    
    
    probability = []
    revSt = {}
    for k,v in start_probability.items():
        revSt.update({v:k})
    # print(transmission_probability[0])
    t1 = transmission_probability[0]
    t2 = transmission_probability[1]
    
    s1 = []
    s2 = []
    for sublist in second1:
        s1.append(sublist[1])
    for sublist in second2:
        s2.append(sublist[1])
    
    
    
    probability.append([[s1[0],max(start_probability.values()), revSt[max(start_probability.values())]]])
    # print(probability)
    # exit()
    
    for step in range(1,len(secondary1)):
        try:
            grab = grabber(re1,probability[step-1][0][2],s1[step], t1)
            probability.append(grab)
        except Exception as e:
            print(e)
            print(t1)
            print(probability[step-1][0][2])
            print(step)
            print(name1)
            print(name2)
            exit(2)
        
    
    # print(probability)
    # print(len(probability),'\t',len(s1))
    # print(counter)
    # print(probability)
    # exit()
    HMM1 = []
    for q in probability:
        HMM1.append(q[0])
    
    
    probability = []
    probability.append([[s2[0],max(start_probability.values()), revSt[max(start_probability.values())]]])
    # print(probability)
    # exit()
    
    for step in range(1,len(secondary2)):
        try:
            grab = grabber(re2,probability[step-1][0][2],s2[step], t2)
            probability.append(grab)
        except Exception as e:
            print(e)
            print(t2)
            print(probability[step-1][0][2])
            print(step)
            
            print(name2)
            exit(2)
        
    
    # print(probability)
    # print(len(probability),'\t',len(s1))
    # print(counter)
    # print(probability)
    # exit()
    HMM2 = []
    for q in probability:
        HMM2.append(q[0])
    
    
    # print('got to here')
    return(HMM1,HMM2)
    


def grabber(re1:set,target:str, amino, t1):
    
    interim = []
    
    targets = {}
    candidates = {}
    # print(t1)
    # for k,v in t1.items():
    #     targets.update(t1[target])
    
    candidates = {}
    for k,v in re1.items():
        multiplier = t1[k]
        
        for key,value in v.items():
            
            for ke,val in multiplier.items():
                
                if key == amino:
                    r = random.uniform(-0.01,0.01)
                    
                    candidates.update({ke:multiplier[k]*(abs(value+r))})
                    
                    
                    
            inv_cand = {a: b for b, a in candidates.items()}
            
    interim.append([key,max(inv_cand.keys())*abs(value +r),inv_cand[max(inv_cand.keys())]])
    # print(interim)
    # exit()
                
    
    return(interim)


def returnEmissionProbability(proteins:list, second1:list, second2:list):
    for item in second1:
        item[0] = item[0].replace('\'','').replace('\"','')
        if '\'' in item[0]:
            print('big problem chief')
            print(item[0])
            exit()
    for item in second2:
        item[0] = item[0].replace('\'','').replace('\"','')
        if '\'' in item[0]:
            print('big problem chief')
            print(item[0])
            exit()
    
    emission1 = {}
    indexer = 0
    for pair in second1:
        try:
            emission1[pair[0]]
            if pair[0] in emission1.keys():
                if pair[1] in emission1[pair[0]].keys():
                    emission1[pair[0]].update({pair[1]:emission1[pair[0]][pair[1]]+1})
                else:
                    emission1[pair[0]].update({pair[1]:+1})
        except Exception as e:
            emission1[pair[0]] = {pair[1]:1}
        indexer += 1
    
    emission2 = {}
    indexer = 0
    for pair in second2:
        try:
            emission2[pair[0]]
            if pair[0] in emission2.keys():
                if pair[1] in emission2[pair[0]].keys():
                    emission2[pair[0]].update({pair[1]:emission2[pair[0]][pair[1]]+1})
                else:
                    emission2[pair[0]].update({pair[1]:+1})
        except Exception as e:
            emission2[pair[0]] = {pair[1]:1}
        indexer += 1
    
    
    
    
    unpackedE1 = []
    unpackedE2 = []
    for key, value in emission1.items():
        unpackedE1.append([key,value])
    for key,value in emission2.items():
        unpackedE2.append([key,value])
    
    unpacked = {}
    temp = {}
    
    emission = {}
    
    for pair in second2:
        try:
            emission[pair[0]]
            if pair[0] in emission.keys():
                if pair[1] in emission[pair[0]].keys():
                    emission[pair[0]].update({pair[1]:emission[pair[0]][pair[1]]+1})
                else:
                    emission[pair[0]].update({pair[1]:+1})
        except Exception as e:
            emission[pair[0]] = {pair[1]:1}
        indexer += 1
    
    
    indexer = 0
    for pair in second1:
        try:
            emission[pair[0]]
            if pair[0] in emission.keys():
                if pair[1] in emission[pair[0]].keys():
                    emission[pair[0]].update({pair[1]:emission[pair[0]][pair[1]]+1})
                else:
                    emission[pair[0]].update({pair[1]:+1})
        except Exception as e:
            emission[pair[0]] = {pair[1]:1}
        indexer += 1
    
    sum1 = 0
    for key,value in emission.items():
        sum1 = 0
        for k,v in value.items():
            sum1+= v
        for k,v in value.items():
            value.update({k:v/sum1})
    # print(emission)
    # print(emission1)
    return(emission1,emission2,emission)
    # exit()
    
    
    # for k1,v1 in emission1.items():
    #     for k2,v2 in emission2.items():
    #         if k1 == k2:
    #             for key1, value1 in v1.items():
    #                 for key2, value2 in v2.items():
    #                     try:
    #                         print(emission1[k2])
    #                         break
    #                     except Exception as e:
    #                         print(e)
    #                         exit()
    #                     print(emission1)
    #                     exit()
    
    
    # for item1 in unpackedE1:
    #     for item2 in unpackedE2:
    #         if item1[0] == item2[0]:
    #             print(item1[1])
    #             # exit() 
    #             try:
    #                 unpacked[item1[0]]
    #                 # if item1[1] 
    #             except Exception as e:
    #                 print(e)
    #                 unpacked[item1[0]] = item1[1]
    #                 print(unpacked)
    #             exit()
                    
    # print(emission1)
    # print(emission2)
    
    # for key,value in emission1:
    #     pass
    
    # exit()
    # return()
    # pass

def secondaryStructure(nwScore:int, alignment:str, structure:str):
    pass
    
def tertiaryStructure(nwScore:int, alignment:str):
    pass


def returnTransmissionProbability(inputTypes:list,secondary1,secondary2):
    types1 = set(secondary1.keys())
    types2 = set(secondary2.keys())
    values1 = list(secondary1.values())
    values2 = list(secondary2.values())
    setValues1 = list(set(values1))
    setValues2 = list(set(values2))
    #creating a nonapplicable value
    for i in range(3):
        values1.insert(0,'null')
        values1.append('null')
    
    transmissionProb1 = {}
    index = 3
    # print(values1)
    # print('\n')
    for value in values1[3:-3]:
        for i in range(3):
            try:
                if values1[index+i] == 'null':
                    continue
            except:
                print(values1)
                print(i)
                print(index)
                exit()
            offset = 1
            if values1[index] == values1[index+i]:
                offset = 3
            else:
                offset = 0.5
            try:
                transmissionProb1[values1[index]]
                if values1[index] in transmissionProb1[values1[index]].keys():
                    # print('UPPER\t%s in %s'%(values1[index],transmissionProb1[values1[index]].keys()))
                    if values1[index] == values1[index+i]:
                        offset = 2
                    transmissionProb1[values1[index]].update({values1[index+i]:transmissionProb1[values1[index]][values1[index+i]]+ (1/(offset+i))})
                else:
                    transmissionProb1[values1[index]].update({values1[index+i]:1/(offset+i)})
            except Exception as e:
                # print(e) 
                # print('initializing new: ', values1[index+i])
                try:
                    transmissionProb1[values1[index]].update({values1[index+i]:1/(offset+i)})
                except Exception as e2:
                    # print(e2,'\tnot in  transmissionProb1[values1[index]]')
                    # print('printing keys from ^^^\t:\t',transmissionProb1.items())
                    transmissionProb1[values1[index]] = {values1[index+i]:1/(offset+i)}
                # print(transmissionProb1)
        
        
        
        for i in range(3):
            if values1[index-i] == 'null':
                continue
            offset = 1
            if values1[index] == values1[index-i]:
                offset = 3
            else:
                offset = 0.5
            try:                    
                
                transmissionProb1[values1[index]]
                if values1[index] in transmissionProb1[values1[index]].keys():
                    # print('LOWER\t%s in %s'%(values1[index],transmissionProb1[values1[index]].keys()))
                    transmissionProb1[values1[index]].update({values1[index-i]:transmissionProb1[values1[index]][values1[index-i]]+ (1/(offset+i))})
                else:
                    transmissionProb1[values1[index]].update({values1[index-i]:1/(offset+i)})
            except Exception as e:
                # print(e) 
                # print('initializing new: ', values1[index-i])
                try:
                    transmissionProb1[values1[index]].update({values1[index-i]:1/(offset+i)})
                except Exception as e2:
                    # print(e2,'\tnot in  transmissionProb1[values1[index]]')
                    # print('printing keys from ^^^\t:\t',transmissionProb1.items())                    
                    transmissionProb1[values1[index]] = {values1[index-i]:1/(offset+i)}
                # print(transmissionProb1)

        index += 1
        
    for i in range(3):
        values2.insert(0,'null')
        values2.append('null')
    
    transmissionProb2 = {}
    index = 3
    # print(values1)
    # print('\n')
    for value in values2[3:-3]:
        for i in range(3):
            if values2[index+i] == 'null':
                continue
            offset = 1
            if values2[index] == values2[index+i]:
                offset = 3
            else:
                offset = 0.5
            try:
                transmissionProb2[values2[index]]
                if values2[index] in transmissionProb2[values2[index]].keys():
                    # print('UPPER\t%s in %s'%(values1[index],transmissionProb1[values1[index]].keys()))
                    if values2[index] == values2[index+i]:
                        offset = 2
                    transmissionProb2[values2[index]].update({values2[index+i]:transmissionProb2[values2[index]][values2[index+i]]+ (1/(offset+i))})
                else:
                    transmissionProb2[values2[index]].update({values2[index+i]:1/(offset+i)})
            except Exception as e:
                # print(e) 
                # print('initializing new: ', values1[index+i])
                try:
                    transmissionProb2[values2[index]].update({values2[index+i]:1/(offset+i)})
                except Exception as e2:
                    # print(e2,'\tnot in  transmissionProb1[values1[index]]')
                    # print('printing keys from ^^^\t:\t',transmissionProb1.items())
                    transmissionProb2[values2[index]] = {values2[index+i]:1/(offset+i)}
                # print(transmissionProb1)
        
        
        
        for i in range(3):
            if values2[index-i] == 'null':
                continue
            offset = 1
            if values2[index] == values2[index-i]:
                offset = 3
            else:
                offset = 0.5
            try:                    
                
                transmissionProb2[values2[index]]
                if values2[index] in transmissionProb2[values2[index]].keys():
                    # print('LOWER\t%s in %s'%(values1[index],transmissionProb1[values1[index]].keys()))
                    transmissionProb2[values2[index]].update({values2[index-i]:transmissionProb2[values2[index]][values2[index-i]]+ (1/(offset+i))})
                else:
                    transmissionProb2[values2[index]].update({values2[index-i]:1/(offset+i)})
            except Exception as e:
                # print(e) 
                # print('initializing new: ', values1[index-i])
                try:
                    transmissionProb2[values2[index]].update({values2[index-i]:1/(offset+i)})
                except Exception as e2:
                    # print(e2,'\tnot in  transmissionProb1[values1[index]]')
                    # print('printing keys from ^^^\t:\t',transmissionProb1.items())                    
                    transmissionProb2[values2[index]] = {values2[index-i]:1/(offset+i)}
                # print(transmissionProb1)

        index += 1
        
        
        # if index>=5:
        #     break
    # transmissionProbA = copy.deepcopy(transmissionProb1)
    sum1 = 0
    for key,value in transmissionProb1.items():
        sum1 = 0
        for k,v in value.items():
            sum1+= v
        for k,v in value.items():
            value.update({k:v/sum1})
    sum1 = 0
    for key,value in transmissionProb2.items():
        sum1 = 0
        for k,v in value.items():
            sum1+= v
        for k,v in value.items():
            value.update({k:v/sum1})
    
    return(transmissionProb1,transmissionProb2)
    # exit()
    # sums = 0
    # for v in transmissionProb1.values():
    #     for value in v.values():
    #         sums+=int(value)
    # print(sums)
    # exit()
    
    #the transmission values will be the frequency with which a type of secondary
    #structure continues or changes into a different one.
    
    
    
    # for value in setValues[]:
    #     pass
    
    
    
    # return(startProbs)


def returnStartProb(inputTypes:list,secondary1,secondary2):
    
    alltypes = list(set(['Strand', '310Helix', 'TurnIV', 'TurnII', 'AlphaHelix', 'GammaInv', 'TurnI','Strand',
     'TurnVIb', 'TurnVIII', 'TurnIV', 'TurnII', 'AlphaHelix', 'TurnI',
     'Strand', 'TurnVIb', 'TurnVIII', 'TurnIV', 'TurnII', 'AlphaHelix', 'TurnI',
    'Strand', 'TurnVIII', 'Disulfide', 'TurnIV', 'TurnII', 'TurnI',
    'Strand', 'GammaClassic', 'TurnVIII', '310Helix', 'TurnIV', 'TurnII', 'AlphaHelix',
    'GammaInv', 'TurnI','Strand', '310Helix', 'TurnIV', 'TurnII', 'TurnI']))
    
    types = []
    
    revSecondary1 = {}
    values = list(secondary1.values())+list(secondary2.values())
    values.sort()
    # print(len(values)) 
    # exit()
    setValues = list(set(values))
    setValues.sort()
    size = 0
    tempsize = 0
    
    startProbs = {}
    
    for value in setValues:
        try:
            tempsize = len(values[0:values.index(setValues[setValues.index(value)+1])]) - len(values[0:values.index(value)])
            startProbs[value] = tempsize/(values.index(value)**2.5 +120)
        except:
            tempsize = len(values) - len(values[0:values.index(value)])
            startProbs[value] = tempsize/(values.index(value)**2.5 +120)
        # print(startProbs)
        size +=tempsize
    
    # print('Im here')
    # print(startProbs)
    
    sum1 = 0
    for key,value in startProbs.items():
        sum1 += value
    for key,value in startProbs.items():
        startProbs.update({key:value/(sum1)})
    
    
    return(startProbs)
    
    # revSecondary2 = {}
    # for key,value in secondary1.items():
    #     revSecondary2[value] = key
    
    
    # for secondaryType in alltypes:
    #     if secondaryType in secondary1.values() or secondaryType in secondary2.values():
    #         types.append(secondaryType)
    # # print(types)
    # for type in types:
    #     pass
    
    # exit()
    # pass


if __name__=='__main__':
    main()


