# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 20:55:03 2021

@author: danie
"""
import collections
from sys import exit
import sys
import json
from tqdm import tqdm
import pandas as pd
from PIL import Image
import copy
import pickle

def progress(count, total, status=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()
def difference(a,b):
    return sum(a[i] != b[i] for i in range(len(a)))


def LD(s, t):
        if s == "":
            return len(t)
        if t == "":
            return len(s)
        if s[-1] == t[-1]:
            cost = 0
        else:
                cost = 1
       
        res = min([LD(s[:-1], t)+1,
               LD(s, t[:-1])+1, 
               LD(s[:-1], t[:-1]) + cost])
    
        return res


def asdfasdf(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += 1

def find_all(string,substring):
    return(list(asdfasdf(string,substring)))

def treeSearch(seq1,seq2,letter,k,kernType):
    # x = list(find_all(seq1[:-len(letter)],letter))  
    x = list(find_all(seq1[:-1],letter))  
    # x = list(find_all(seq1,letter))    
    # y = list(find_all(seq2[:-len(letter)],letter))    
    y = list(find_all(seq2[:-1],letter))    
    # y = list(find_all(seq2,letter))
    if x != [] or y != []:
        memx = []
        memy = []
        for i in x:
            memx.append(seq1[i+len(letter)])
        for i in y:
            memy.append(seq2[i+len(letter)])
        memtemp = list(set(memx))+list(set(memy))
        mem = [item for item, count in collections.Counter(memtemp).items() if count > 1]
        if mem != [] and len(mem) == 1:
            return(treeSearch(seq1,seq2,letter+mem[0],k,kernType))
            # treeSearch(seq1,seq2,letter+mem[0],k,kernType)
        elif mem !=[] and len(mem) > 1:
            for l in mem:
                return(treeSearch(seq1,seq2,letter+l,k,kernType))
        else:
            if len(letter) == k:
                return(letter)
            else:
                pass
    elif x != [] or y != []:
        # if len(letter) == 3:
        #   return(letter)  
        return []
            


def spectrumkernel(seq1:str,seq2:str,k:int):
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    frames = []
    for l in alphabet:
        n = 0
        kernType = 'spec'
        frames.append(treeSearch(seq1,seq2,l,k,kernType))

    return(frames)
        
            
def count(lCount,rCount):
    # if lCount * rCount != 0:
        return(lCount + rCount)

def realmismatchKernel(seq1:str, seq2:str, k:int, m):
    kmers1 = []
    kmers2 = []
    for i in range(0,len(seq1)-k+1):
        kmers1.append(seq1[i:i+k])
    for i in range(0,len(seq2)-k+1):
        kmers2.append(seq2[i:i+k])
    kmers1.sort()
    kmers2.sort()
    sharedKmers = sorted(set(kmers1+kmers2))
    merged = {}
    lMerged = {}
    rMerged = {}
    mergedDict = {}
    lCount = {}
    rCount = {}
    
    for i in range(len(sharedKmers)):
        if sharedKmers[i] in kmers1:
            try:
                lCount[i] = lCount[i]+1
            except:
                lCount[i] = 1
            kmers1.pop(kmers1.index(sharedKmers[i]))
        else:
            lCount[i] = 0
        if sharedKmers[i] in kmers2:
            try:
                rCount[i] = rCount[i]+1
            except:
                rCount[i] = 1
            kmers2.pop(kmers2.index(sharedKmers[i]))
        else:
            rCount[i] = 0
    # for i in range(len(sharedKmers)):
    #     if lCount[i]!= 0 and rCount[i] != 0:
    #         print(lCount[i],'\t',sharedKmers[i],'\t',rCount[i])
    # exit()
    if type(m) == int and m > 0:
        kmers1 = []
        kmers2 = []
        # for i in range(len(sharedKmers)):
        #     print(lCount[i],'\t',sharedKmers[i],'\t',rCount[i])
        # exit()
        for i in range(0,len(seq1)-k+1):
            kmers1.append(seq1[i:i+k])
        for i in range(0,len(seq2)-k+1):
            kmers2.append(seq2[i:i+k])
    
        for s in range(len(sharedKmers)):
            for i in range(len(sharedKmers)):
                if s != i and difference(sharedKmers[s],sharedKmers[i]) <= m and difference(sharedKmers[s],sharedKmers[i]) > 0:
                    try:
                        lMerged[sharedKmers[i]+sharedKmers[s]] = lCount[i] + lCount[s]
                        rMerged[sharedKmers[i]+sharedKmers[s]] = rCount[i] + rCount[s]
                    except:
                        print('error')
                        exit()
                    #     lMerged[sharedKmers[i]+sharedKmers[s]] = lMerged[sharedKmers[i]+sharedKmers[s]] + lCount[i] + lCount[s]
                    #     rMerged[sharedKmers[i]+sharedKmers[s]] = rMerged[sharedKmers[i]+sharedKmers[s]] + rCount[i] + rCount[s]
                        # print(lCount)
                        
                        # print(lMerged)
                        # print(sharedKmers[i]+sharedKmers[s])
                        # # lMerged[sharedKmers[i]+sharedKmers[s]]
                        # exit()
                    # lCount.pop(i)
                    # rCount.pop(i)
                    # lCount.pop(s)
                    # rCount.pop(s)
        merged = []
        lmerged = {}
        rmerged = {}
        # for i in lMerged.keys():
        #     if (lMerged[i] != 0 and rMerged[i]!=0):
        #         print(lMerged[i],'\t',rMerged[i])
        for key,value in lMerged.items():
            # print('loop triggered')
            for key1, value1 in rMerged.items():
                if lMerged[key1] == rMerged[key1] == 0:
                    # print('its 0')
                    pass
                if key[k:] == key1[:k]:
                    # print('keys merged')
                    merged.append(key[3:])
                    lmerged[merged[-1]] = lMerged[key] + lMerged[key1]
                    rmerged[merged[-1]] = rMerged[key] + rMerged[key1]
                    lMerged[key1] = 0
                    rMerged[key1] = 0
                    lMerged[key] = 0
                    rMerged[key] = 0
        lpop = []
        rpop = []
        for key,value in lmerged.items():
            if value == 0:
                lpop.append(key)
        for i in lpop:
            lmerged.pop(i)
        for key,value in rmerged.items():
            if value == 0:
                rpop.append(key)
        for i in rpop:
            rmerged.pop(i)
        # print(rmerged)
        # print(lmerged)
        # exit()
        lkeep = 0
        rkeep = 0
        for i in lmerged.keys():
            # print('SHAREDKMERS')
            try:
                lkeep += count(lmerged[i],rmerged[i])
            except:
                pass
        for j in rmerged.keys():
            try:
                rkeep += count(lmerged[i],rmerged[j])
            except:
                pass
        # print(rmerged)
        # print('returning value')
        if lkeep > rkeep:
            return lkeep
        else:
            return rkeep
    # if type(m) == int and m > 0:
    #     for s in range(len(sharedKmers)):
    #         # progress(s,len(sharedKmers),'kmer mismatch comparison')
    #         for i in range(len(kmers1)):
    #             if difference(sharedKmers[s],kmers1[i]) != 0 and difference(sharedKmers[s],kmers1[i]) <= m:
    #                 try:
    #                     merged[sharedKmers[s]+kmers1[i]] = count(lCount[i],rCount[i])
    #                     lCount[i] = 0
    #                     rCount[i] = 0
    #                 except:
    #                     merged[sharedKmers[s]+kmers1[i]] = merged[sharedKmers[s]+kmers1[i]] + count(lCount[i],rCount[i])
    #                     lCount[i] = 0
    #                     rCount[i] = 0
    #             for key,value in merged.items():
    #                 if lCount[i] == rCount[i] == 0:
    #                     break
    #                 if difference(key[:k],kmers1[i]) > 0 and difference(key[:k],kmers1[i]) <= m:
    #                     try:
    #                         merged[key] = count(lCount[i],rCount[i])
    #                         lCount[i] = 0
    #                         rCount[i] = 0
    #                     except:
    #                         merged[key] = merged[key] + count(lCount[i],rCount[i])
    #                         lCount[i] = 0
    #                         rCount[i] = 0
    #         for i in range(len(kmers2)):
    #             if difference(sharedKmers[s],kmers2[i]) != 0 and difference(sharedKmers[s],kmers2[i]) < m:
    #                 try:
    #                     merged[sharedKmers[s]+kmers2[i]] = count(lCount[i],rCount[i])
    #                     lCount[i] = 0
    #                     rCount[i] = 0
    #                 except:
    #                     merged[sharedKmers[s]+kmers2[i]] = merged[sharedKmers[s]+kmers2[i]] + count(lCount[i],rCount[i])
    #                     lCount[i] = 0
    #                     rCount[i] = 0
    #             for key,value in merged.items():
    #                 if lCount[i] == rCount[i] == 0:
    #                     break
    #                 if difference(key[:k],kmers2[i]) > 0 and difference(key[:k],kmers2[i]) <= m:
    #                     try:
    #                         merged[key] = count(lCount[i],rCount[i])
    #                         lCount[i] = 0
    #                         rCount[i] = 0
    #                     except:
    #                         merged[key] = merged[key] + count(lCount[i],rCount[i])
    #                         lCount[i] = 0
    #                         rCount[i] = 0
    # print(merged)
    # print('got to the end')
    # for i in range(len(sharedKmers)):
    #     print(lCount[i],'\t',sharedKmers[i],'\t',rCount[i])
    # exit()
    # if type(m) == int and m > 0:
    #     for i in range(len(sharedKmers)):
    #         progress(i,len(sharedKmers))
    #         for j in range(len(sharedKmers)):
    #             if differences(sharedKmers[i],sharedKmers[j]) <= m and differences(sharedKmers[i],sharedKmers[j]) > 0:
    #                 merged.append(sharedKmers[i])
    #                 merged.append(sharedKmers[j])
    #                 # mergedDict[i] = sharedKmers[i]
    #                 # mergedDict[j] = sharedKmers[i]
    #                 lCount[i] = lCount[i] + lCount[j]
    #                 rCount[i] = rCount[i] + rCount[j]
    #             if merged != []:
    #                 for mrgd in merged:
    #                     if len(merged) > 300:
    #                         print(merged)
    #                         exit()
    #                     if differences(sharedKmers[i],mrgd) <= m and differences(sharedKmers[i],mrgd) > 0:
    #                         merged.append(sharedKmers[i])
    #                         mergedDict[i] = sharedKmers[i]
    #                         lCount[merged.index(mrgd)] = lCount[merged.index(mrgd)] + lCount[i]
    #                         rCount[merged.index(mrgd)] = rCount[merged.index(mrgd)] + rCount[i]
    #                         break
    #                     if differences(sharedKmers[j],mrgd) <= m and differences(sharedKmers[j],mrgd) > 0:
    #                         merged.append(sharedKmers[j])
    #                         try:
    #                             lCount[merged.index(mrgd)] = lCount[merged.index(mrgd)] + lCount[j]
    #                         except:
    #                             print('error occurred')
    #                             print(merged)
    #                             exit()
    #                         rCount[merged.index(mrgd)] = rCount[merged.index(mrgd)] + rCount[j]
    #                         break
    
    kValue = 0
    for i in range(len(sharedKmers)):
        if lCount[i] * rCount[i] > 0:
            kValue += count(lCount[i],rCount[i])
    return(kValue)


import numpy as np
def mismatchKernel(seq1:str, seq2:str, k:int, m=' '):
    kmers1 = []
    kmers2 = []
    for i in range(0,len(seq1)-k+1):
        kmers1.append(seq1[i:i+k])
    for i in range(0,len(seq2)-k+1):
        kmers2.append(seq2[i:i+k])
    kmers1.sort()
    kmers2.sort()
    sharedKmers = sorted(set(kmers1+kmers2))
    merged = {}
    lMerged = {}
    rMerged = {}
    mergedDict = {}
    lCount = {}
    rCount = {}
    
    for i in range(len(sharedKmers)):
        if sharedKmers[i] in kmers1:
            try:
                lCount[i] = lCount[i]+1
            except:
                lCount[i] = 1
            kmers1.pop(kmers1.index(sharedKmers[i]))
        else:
            lCount[i] = 0
        if sharedKmers[i] in kmers2:
            try:
                rCount[i] = rCount[i]+1
            except:
                rCount[i] = 1
            kmers2.pop(kmers2.index(sharedKmers[i]))
        else:
            rCount[i] = 0
    # for i in range(len(sharedKmers)):
    #     if lCount[i]!= 0 and rCount[i] != 0:
    #         print(lCount[i],'\t',sharedKmers[i],'\t',rCount[i])
    # exit()
    if type(m) == int and m > 0:
        kmers1 = []
        kmers2 = []
        # for i in range(len(sharedKmers)):
        #     print(lCount[i],'\t',sharedKmers[i],'\t',rCount[i])
        # exit()
        for i in range(0,len(seq1)-k+1):
            kmers1.append(seq1[i:i+k])
        for i in range(0,len(seq2)-k+1):
            kmers2.append(seq2[i:i+k])
    
        for s in range(len(sharedKmers)):
            for i in range(len(sharedKmers)):
                if s != i and difference(sharedKmers[s],sharedKmers[i]) <= m and difference(sharedKmers[s],sharedKmers[i]) > 0:
                    lMerged[sharedKmers[i]+sharedKmers[s]] = lCount[i] + lCount[s]
                    rMerged[sharedKmers[i]+sharedKmers[s]] = rCount[i] + rCount[s]
        merged = []
        lmerged = {}
        rmerged = {}
        for key,value in lMerged.items():
            for key1, value1 in rMerged.items():
                if lMerged[key1] == rMerged[key1] == 0:
                    continue
                if key[k:] == key1[:k]:
                    # print('keys merged')
                    merged.append(key[3:])
                    lmerged[merged[-1]] = lMerged[key] + lMerged[key1]
                    rmerged[merged[-1]] = rMerged[key] + rMerged[key1]
                    lMerged[key1] = 0
                    rMerged[key1] = 0
                    lMerged[key] = 0
                    rMerged[key] = 0
        lpop = []
        rpop = []
        for key,value in lmerged.items():
            if value == 0:
                lpop.append(key)
        for i in lpop:
            lmerged.pop(i)
        for key,value in rmerged.items():
            if value == 0:
                rpop.append(key)
        for i in rpop:
            rmerged.pop(i)
        lkeep = 0
        rkeep = 0
        for i in lmerged.keys():
            # print('SHAREDKMERS')
            try:
                lkeep += count(lmerged[i],rmerged[i])
            except:
                pass
        for j in rmerged.keys():
            try:
                rkeep += count(lmerged[i],rmerged[j])
            except:
                pass
        # print(rmerged)
        # print('returning value')
        if lkeep > rkeep:
            return lkeep
        else:
            return rkeep
    kValue = 0
    for i in range(len(sharedKmers)):
        if lCount[i] * rCount[i] > 0:
            kValue += count(lCount[i],rCount[i])
    return(kValue)
def forcedMatchLogic(i:int,j:int,tempHumanS:list,tempHumanE:list,tempFligmaS:list,tempFligmaE:list):
    for m in range(len(tempHumanS)):
        if((tempFligmaS[m]<= i <= tempFligmaE[m]) and (tempHumanS[m]<= j <= tempHumanE[m])):
            return(1)
        else:
            return(0)


def BLOSUM62(BLOSUM:str,key:str,first:chr,second:chr):
    
    first = first.upper()
    second = second.upper()
    
    firstPos = 0
    secondPos = 0
    keyTranspose = key.replace(' ','')
    # print(keyTranspose)
    wrapery =len(keyTranspose)
    wraperx = len(key)+3
    # print(wraperx,wrapery)
    for i in range(len(key)):
        if(key[i] == first):
            firstPos = i
    for i in range(wrapery):
        if(keyTranspose[i] == second):
            secondPos = i
            # print('Second is: ',keyTranspose[i])
    # print(firstPos,'\t',secondPos)
    # print(secondPos)
    # print('BLOSUM[firstPos +wraper*secondPos-1]','\n',BLOSUM[firstPos +wraperx*secondPos-1])
    # exit()
    if(BLOSUM[firstPos +wraperx*secondPos-2] == '1'):
        return(11)
    if(BLOSUM[firstPos +wraperx*secondPos-2] == '-'):
        return(-1*int(BLOSUM[firstPos+wraperx*secondPos-1]))
    else:

        return(int(BLOSUM[firstPos+wraperx*secondPos-1]))


def nw(x, y, BLOSUM,key,gap, tempHumanS=[],tempHumanE=[],tempFligmaS=[],tempFligmaE = []):
    flagAligner = 0
    if(tempHumanS != []):
        flagAligner = 1
    lengX = len(x)
    lengY = len(y)
    # Optimal score at each possible pair of characters.
    Needleman = np.zeros((lengX + 1, lengY + 1))
    Needleman[:,0] = np.linspace(0, -lengX, lengX + 1)
    Needleman[0,:] = np.linspace(0, -lengY, lengY + 1)
    # Pointers to trace through an optimal aligment.
    P = np.zeros((lengX + 1, lengY + 1))
    P[:,0] = 3
    P[0,:] = 4
    # Temporary scores.
    t = np.zeros(3)
    for i in range(lengX):
        for j in range(lengY):
            # if((not(x[i] == y[j]))or((flagAligner == 1)and(forcedMatchLogic(i,j,tempHumanS,tempHumanE,tempFligmaS,tempFligmaE)==1))):
            #     print('seq force aligned, I: %i J: %i'%(i,j))
            if((x[i] == y[j])or((flagAligner == 1)and(forcedMatchLogic(i,j,tempHumanS,tempHumanE,tempFligmaS,tempFligmaE)==1))):
                t[0] = Needleman[i,j] + BLOSUM62(BLOSUM,key,x[i],y[j])
            else:
                t[0] = Needleman[i,j] + BLOSUM62(BLOSUM,key,x[i],y[j])
            t[1] = Needleman[i,j+1] - gap
            t[2] = Needleman[i+1,j] - gap
            tmax = np.max(t)
            Needleman[i+1,j+1] = tmax
            if((t[0] == tmax)or(forcedMatchLogic(i,j,tempHumanS,tempHumanE,tempFligmaS,tempFligmaE)==1)):
                P[i+1,j+1] += 2
            elif t[1] == tmax:
                P[i+1,j+1] += 3
            elif t[2] == tmax:
                P[i+1,j+1] += 4
    # Trace through an optimal alignment.
    i = lengX
    j = lengY
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i,j] in [2, 5, 6, 9]:
            rx.append(x[i-1])
            ry.append(y[j-1])
            i -= 1
            j -= 1
        elif P[i,j] in [3, 5, 7, 9]:
            rx.append(x[i-1])
            ry.append('-')
            i -= 1
        elif P[i,j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(y[j-1])
            j -= 1
    # Reverse the strings.
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    
    # print(rx)
    # print()
    # print(ry)
    alignment = ''
    
    for i in range(len(rx)):
        alignment = alignment + rx[i] + '\t'+ ry[i] + '\n'
    
    return(Needleman[lengX,lengY],alignment)
    
    # similarity = P[lengX,lengY]
    # return similarity


def readSCOP(search:str, fileIn:str):
    inFile = open(fileIn, 'r')
    SCOPentry = ''
    for line in inFile:
        if line[0:line.index(' ')] == search:
            SCOPentry = line.rstrip()
            break
    inFile.close()
    return(SCOPentry)
            
def getAllFamily(SCOP_entry:str, fileIn1:str, fileIn2:str):
    flag = 0
    inFile = open(fileIn1, 'r')
    familyMembers = ''
    allFamily = 'Family members of %s: \n'%SCOP_entry[:SCOP_entry.index(' ')]
    family = SCOP_entry[SCOP_entry.index('FA=')+3:]
    for line in inFile:
        if flag == 1:
            familyMembers +=line +'\n'
            flag = 0
        if family in line:
            allFamily += line.rstrip()+'\n'
            flag = 1
    inFile.close()
    flag = 0
    
    inFile = open(fileIn2, 'r')
    allSuperFamily = 'Super-Family members of %s: \n'%SCOP_entry[:SCOP_entry.index(' ')]
    superFamilyMembers=''
    superFamily = SCOP_entry[SCOP_entry.index('SF=')+3:SCOP_entry.index('FA=')+3]
    for line in inFile:
        if flag == 1:
            superFamilyMembers += line +'\n'
            flag = 0
        if superFamily in line:
            allSuperFamily += line.rstrip()
            flag = 1
    inFile.close()
    
    return(allFamily,familyMembers,allSuperFamily,superFamilyMembers)
    
        
        ##MORE GOES HERE
# # def generatePairwiseandCall(family, superfamily,sequences):
# #     output = np.zeros(len(sequences.keys()),len(sequences.keys()))
    
    

def main():
    names = readSCOP('8045611','scop-cla-latest.txt')
    names = getAllFamily(names,'scop_fa_represeq_lib_latest.fa','scop_sf_represeq_lib_latest.fa')
    json.dump(names, open('SCOP_entry.txt','w'))
    writeNames = open('Human-Readable SCOP_entry.txt', 'w')
    writeNames.write(names[0])
    writeNames.write(names[2])
    writeNames.close()
    
    
    family = names[0]
    superfamily = names[2]
    json.dump(family, open("tempFamilyStorage.txt",'w'))
    json.dump(superfamily, open("tempSuperFamilyStorage.txt",'w'))
    
    
    
    tempfamilyMembers = copy.deepcopy(names[1])
    familyMembers = tempfamilyMembers.split()
    tempsuperFamilyMembers = copy.deepcopy(names[3])
    superFamilyMembers = tempsuperFamilyMembers.split()
    
    
    # print(familyMembers[0])
    # print(familyMembers[1])
    # print(familyMembers)
    # exit()
    b = np.zeros((len(familyMembers),len(familyMembers)))
    # print(familyMembers)
    # exit()
    for i in (range(len(familyMembers))):
        for j in range(len(familyMembers)):
            b[i,j] = mismatchKernel(familyMembers[i],familyMembers[j],2)
    # print(b)
    
    key = '.  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *'
    BLOSUM = ''
    stepper = 0
    freeTemp = ''
    fileIn = open('BLOSUM62.txt','r')
    for line in fileIn:
        freeTemp = line
        BLOSUM = BLOSUM + line +'\n'
        
    # for i in range(len(list(names.values()))):
    #     nameIndex[list(names.values())[i]] = index
    #     revNameIndex[index] = list(names.values())[i]
    #     index += 1
    index = len(familyMembers)
    # print(familyMembers)
    # print(names)
    # exit()
    # seq1 = familyMembers[4]
    # seq2 = familyMembers[5]
    # string = nw(seq1,seq2,BLOSUM,key,5)
    # fileOut = open('testtest.txt','w')
    # fileOut.write(str(string))
    # fileOut.close()
    # exit()
    alignments = []
    init_i = -1
    # print(range(len(familyMembers)))
    # exit()
    nwMat = np.zeros((index,index))
    # print(familyMembers)
    # exit()
    # z=0
    # for i in familyMembers:
    #     print(i)
    #     z+=1
    #     print(z)
    #     print('\n')
    # exit()
    
    for i in tqdm(range(len(familyMembers))):
        for j in range(len(familyMembers)):
            seq1 = familyMembers[i]
            seq2 = familyMembers[j]
            temp =  nw(seq1,seq2,BLOSUM,key,5)
            new_temp = copy.deepcopy(temp)
            nwMat[i,j] = new_temp[0]
            alignments.append(temp[1])
            # if i == 4:
            #     print(familyMembers[i])
            #     print(familyMembers[2])
            #     print(nw(familyMembers[i],familyMembers[2],BLOSUM,key,5))
            #     exit()
    # print(nwMat)
    
    # exit()
    # exit()
    # print(alignments)
    # print(alignments[0])
    # exit()
    # json.dump(alignments, open('Alignments.json','w'))
    # print(type(alignments))
    # print(type(nwMat))
    # print(np.shape(nwMat))
    # print(len(alignments))
    # exit()
    with open("alignments.txt", "wb") as fp:   #Pickling
        pickle.dump(alignments, fp)
        
    
    np.save('nwMattemp',nwMat)
    nwMat = np.load('nwMattemp.npy')
    nwMat =nwMat*100
    img = Image.fromarray(nwMat)
    if img.mode != 'RGB':
        img = img.convert('RGB')
    img.save('NWA on SCOP proteins.png')
    nwMat = nwMat/100
    import csv
    f = open('HumanReadableFirstRowNWAonSCOP.csv',"w")
    writer = csv.writer(f)
    for row in nwMat:
        writer.writerow(row)

    f.close()
    
    # f = open('')
    exit()
    # np.set_printoptions(suppress=True)
    # mismatch - oldmismatch
    # np.save('tempFileMismatches',mismatch)
    # np.set_printoptions(suppress=True)
    # oldMismatch = np.load('tempFileMismatches.npy')
    # np.set_printoptions(suppress=True)
    # np.savetxt('humanReadableMismatch.txt',mismatch)    
    # np.set_printoptions(suppress=True)
    # mismatch = np.random.rand(index,index)
    # mismatch = mismatch *5
    # new_p = Image.fromarray(mismatch)
    # if new_p.mode != 'RGB': 
    #     new_p = new_p.convert('RGB')
    # new_p.save("jpgMismatch.jpeg")
    # import csv
    # mismatch = mismatch/5
    # f = open('HumanReadable.csv',"w")
    # writer = csv.writer(f)
    # for row in mismatch:
    #     writer.writerow(row)

    # f.close()
    # print(mismatch)
    exit()
    # k = spectrumkernel(seq1,seq2,3)
    # k = treeSearch(seq1,seq2,'A',3,False)
    # k = [x for x in k if x != None]
    # print(k)
    # mismatch = mismatchKernel(seq1,seq2,3,1)
    # print(mismatch)

# def mismatchKernel(seq1:str, seq2:str, k:int, m):
#     kmers1 = []
#     kmers2 = []
#     for i in range(0,len(seq1)-k+1):
#         kmers1.append(seq1[i:i+k])
#     for i in range(0,len(seq2)-k+1):
#         kmers2.append(seq2[i:i+k])
#     kmers1.sort()
#     kmers2.sort()
#     sharedKmers = sorted(set(kmers1+kmers2))
#     print(sharedKmers)
#     exit()
#     merged = {}
#     mergedDict = {}
#     lCount = {}
#     rCount = {}
    
#     for i in range(len(sharedKmers)):
#         if sharedKmers[i] in kmers1:
#             try:
#                 lCount[i] = lCount[i]+1
#             except:
#                 lCount[i] = 1
#             kmers1.pop(kmers1.index(sharedKmers[i]))
#         else:
#             lCount[i] = 0
#         if sharedKmers[i] in kmers2:
#             try:
#                 rCount[i] = rCount[i]+1
#             except:
#                 rCount[i] = 1
#             kmers2.pop(kmers2.index(sharedKmers[i]))
#         else:
#             rCount[i] = 0
#     print(kmers1)
#     print(kmers2)
#     exit()
#     if type(m) == int and m > 0:
#         for s in range(len(sharedKmers)):
#             print(sharedKmers)
#             print(kmers1)
#             exit()
#             for i in range(len(kmers1)):
#                 if difference(sharedKmers[s],kmers1[i]) != 0 and difference(sharedKmers[s],kmers1[i]) < m:
#                     print('seq1 threshold met')
#                     try:
#                         merged[sharedKmers[s]+kmers1[i]] = count(lCount[i],rCount[i])
#                         lCount[i] = 0
#                         rCount[i] = 0
#                     except:
#                         merged[sharedKmers[s]+kmers1[i]] = merged[sharedKmers[s]+kmers1[i]] + count(lCount[i],rCount[i])
#                         lCount[i] = 0
#                         rCount[i] = 0
#                 for key,value in merged.items():
#                     if lCount[i] == rCount[i] == 0:
#                         break
#                     if difference(key[:k],kmers1[i]) > 0 and difference(key[:k],kmers1[i]) <= m:
#                         try:
#                             merged[key] = count(lCount[i],rCount[i])
#                             lCount[i] = 0
#                             rCount[i] = 0
#                         except:
#                             merged[key] = merged[key] + count(lCount[i],rCount[i])
#                             lCount[i] = 0
#                             rCount[i] = 0
#             for i in range(len(kmers2)):
#                 print('seq2 threshold met')
#                 if difference(sharedKmers[s],kmers2[i]) != 0 and difference(sharedKmers[s],kmers2[i]) < m:
#                     try:
#                         merged[sharedKmers[s]+kmers2[i]] = count(lCount[i],rCount[i])
#                         lCount[i] = 0
#                         rCount[i] = 0
#                     except:
#                         merged[sharedKmers[s]+kmers2[i]] = merged[sharedKmers[s]+kmers2[i]] + count(lCount[i],rCount[i])
#                         lCount[i] = 0
#                         rCount[i] = 0
#                 for key,value in merged.items():
#                     if lCount[i] == rCount[i] == 0:
#                         break
#                     if difference(key[:k],kmers2[i]) > 0 and difference(key[:k],kmers2[i]) <= m:
#                         try:
#                             merged[key] = count(lCount[i],rCount[i])
#                             lCount[i] = 0
#                             rCount[i] = 0
#                         except:
#                             merged[key] = merged[key] + count(lCount[i],rCount[i])
#                             lCount[i] = 0
#                             rCount[i] = 0
#     print(merged)
#     exit()
# def THIS THING WORKS(seq1,seq2,letter,k,kernType):
#     # x = list(find_all(seq1[:-len(letter)],letter))  
#     x = list(find_all(seq1[:-1],letter))  
#     # x = list(find_all(seq1,letter))    
#     # y = list(find_all(seq2[:-len(letter)],letter))    
#     y = list(find_all(seq2[:-1],letter))    
#     # y = list(find_all(seq2,letter))
#     if x != [] or y != []:
#         memx = []
#         memy = []
#         for i in x:
#             memx.append(seq1[i+len(letter)])
#         for i in y:
#             memy.append(seq2[i+len(letter)])
#         memtemp = list(set(memx))+list(set(memy))
#         mem = [item for item, count in collections.Counter(memtemp).items() if count > 1]
#         if mem != [] and len(mem) == 1:
#             return(treeSearch(seq1,seq2,letter+mem[0],k,kernType))
#             # treeSearch(seq1,seq2,letter+mem[0],k,kernType)
#         elif mem !=[] and len(mem) > 1:
#             for l in mem:
#                 return(treeSearch(seq1,seq2,letter+l,k,kernType))
#         else:
#             if len(letter) == k:
#                 return(letter)
#             else:
#                 pass
#     elif x != [] or y != []:
#         # if len(letter) == 3:
#         #   return(letter)  
#         return []
    

#standard boilerplate to set 'main' as starting function
if __name__=='__main__':
    main()