# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 21:04:22 2022

@author: danie
"""


def main():
    print('TESTTESTTEST')
    key = '.  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *'
    BLOSUM = ''
    stepper = 0
    freeTemp = ''
    fileIn = open('BLOSUM62.txt','r')
    for line in fileIn:
        freeTemp = line
        BLOSUM = BLOSUM + line +'\n'
    seq1 = 'GCATGCG'
    seq2 = 'GATTACA'
    nwMat = nw(seq1,seq2,BLOSUM,key,5)
    print(nwMat)

import sys
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


import numpy as np
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
    print(Needleman)
    # exit()
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
    print(P)
    print(rx)
    print(ry)
    # exit()
    alignment = ''
    
    for i in range(len(rx)):
        alignment = alignment + rx[i] + '\t'+ ry[i] + '\n'
        
    print(alignment)
    print(Needleman)
    np.set_printoptions(threshold=sys.maxsize)
    print(P)
    # exit()
    similarity = P[lengX,lengY]
    print(similarity)
    print(lengX,lengY)
    # exit()
    return alignment
    # return similarity


#standard boilerplate to set 'main' as starting function
if __name__=='__main__':
    main()