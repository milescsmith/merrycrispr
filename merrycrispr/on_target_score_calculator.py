#!/usr/bin/env python3
'''
Calculates the on-target score for an sgRNA
Input: 30mer
Output: On-target score
Run as: python on_target_score_calculator.py <30mer> 
'''
import math


def calc_score(s):
    s_list = list(s)
    s_20mer = s[4:24]
    score = 0.597636154
    gc = s_20mer.count('G')+s_20mer.count('C')
    gc_low = -0.202625894
    gc_high = -0.166587752
    if gc < 10:
        gc_val = abs(gc-10)
        score = score+(gc_val*gc_low)
    elif gc > 10:
        gc_val = gc-10
        score = score+(gc_val*gc_high)
    sing_nuc_hash = {'G2':-0.275377128,'A3':-0.323887456,'C3':0.172128871,'C4':-0.100666209,'C5':-0.20180294, \
                    'G5':0.245956633,'A6':0.036440041,'C6':0.098376835,'C7':-0.741181291,\
                    'G7':-0.393264397,'A12':-0.466099015,'A15':0.085376945,'C15':-0.013813972,\
                    'A16':0.272620512,'C16':-0.119022648,'T16':-0.285944222,'A17':0.097454592,\
                    'G17':-0.17554617,'C18':-0.345795451,'G18':-0.678096426,'A19':0.22508903,\
                    'C19':-0.507794051,'G20':-0.417373597,'T20':-0.054306959,'G21':0.379899366,\
                    'T21':-0.090712644,'C22':0.057823319,'T22':-0.530567296,'T23':-0.877007428,\
                    'C24':-0.876235846,'G24':0.278916259,'T24':-0.403102218,'A25':-0.077300704,\
                    'C25':0.287935617,'T25':-0.221637217,'G28':-0.689016682,'T28':0.117877577,\
                    'C29':-0.160445304,'G30':0.386342585}
    dinuc_hash = {'GT2':-0.625778696,'GC5':0.300043317,'AA6':-0.834836245,'TA6':0.760627772,'GG7':-0.490816749,'GG12':-1.516907439,'TA12':0.7092612,'TC12':0.496298609,'TT12':-0.586873894,'GG13':-0.334563735,'GA14':0.76384993,'GC14':-0.53702517,'TG17':-0.798146133,'GG19':-0.66680873,'TC19':0.353183252,'CC20':0.748072092,'TG20':-0.367266772,'AC21':0.568209132,'CG21':0.329072074,'GA21':-0.836456755,'GG21':-0.782207584,'TC22':-1.029692957,'CG23':0.856197823,'CT23':-0.463207679,'AA24':-0.579492389,'AG24':0.649075537,'AG25':-0.077300704,'CG25':0.287935617,'TG25':-0.221637217,'GT27':0.117877577,'GG29':-0.697740024}
    for i,nuc in enumerate(s_list):
        key = nuc+str(i+1)
        if key in sing_nuc_hash:
            nuc_score = sing_nuc_hash[key]
        else:
            nuc_score = 0
        #nuc_score = score_mat[i,nuc_hash[nuc]]
        score = score+nuc_score
        if i<29:
            dinuc = nuc+s[i+1]+str(i+1)
            if dinuc in dinuc_hash:
                score = score+dinuc_hash[dinuc]
    partial_score = math.e**-score
    final_score = 1/(1+partial_score)
    return final_score
