#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 13:45:12 2021

@author: elio
"""

import statsmodels.api as sm
import statsmodels as sm
import numpy as np
import numpy.matlib as mtlb
from scipy.io import loadmat
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})

#%% define functions of interest


def cmpt_phase_locked(mat, subjdim=1, srate=30):
    
    if subjdim == 0:
        
        mat = mat.T
              
    lsignal = mat.shape[0]
    nsubjs = mat.shape[1]
    
    cmplx_longform = np.fft.fft(mat, axis=0)/lsignal
    cmplx_shortform = cmplx_longform[0:int(np.ceil(lsignal/2)), :]
    
    freqs = np.fft.fftfreq(lsignal, 1/srate)
    freqs = freqs[0:int(np.ceil(lsignal/2))]
    
    cmplx_shortform = cmplx_shortform
    phaselocksum = np.abs(cmplx_shortform.sum(axis=1))/nsubjs
    
    return phaselocksum, freqs, cmplx_shortform



def autoregress_model(mat, nperms=10000, subjdim=1, srate=30, ptile=95):
    
    if subjdim == 0:
        
        mat = mat.T
        subjdim = 1
    
    
    lsignal = mat.shape[0]
    nsubjs = mat.shape[1]
    
    perm_dat = np.zeros([lsignal, nsubjs, nperms])
    
    # loop to create AR1 data from actual data
    for isubj in range(nsubjs):
        
        this_subj_ts = mat[:, isubj]
        AR1_mdl = sm.tsa.arima.model.ARIMA(this_subj_ts, order=(1, 0, 0))
        fit_AR1mdl = AR1_mdl.fit()

        # generate surrogate
        for iperm in range(nperms):
            
            tsim = sm.tsa.arima_process.arma_generate_sample([1, -fit_AR1mdl.arparams[0]], [1], 
                                                             nsample=lsignal, 
                                                             scale=fit_AR1mdl.resid.std())

            perm_dat[:, isubj, iperm] = tsim
            
        print('(AR1) finished with subj ' + str(isubj) + '/' + str(nsubjs))
            

    # loop to perform phase locked sum on permutations
    PL_perm = np.zeros((int(np.ceil(lsignal/2)), nperms+1))
    # test Holmes correction as well
    PL_perm_holmes = np.zeros((int(np.ceil(lsignal/2)), nperms+1))

    for iperm in range(nperms +1):
        
        # add the real data into the permuted, in order not to have a complete 0 value as pval
        
        if iperm == 0:
            
            this_mat = mat
            phaselocksum, _, _ = cmpt_phase_locked(this_mat, subjdim=subjdim, srate=srate)
            real_dat = np.copy(phaselocksum)

        else:
                    
            this_mat = perm_dat[:, :, iperm-1]
            phaselocksum, _, _ = cmpt_phase_locked(this_mat, subjdim=subjdim, srate=srate)

     
        PL_perm[:, iperm] = phaselocksum
        PL_perm_holmes[:, iperm] = phaselocksum.max()
        
        
    bonf_ptile = 100 - 5/len(phaselocksum)   
    thresh_sign = {'uncorrected' : np.percentile(PL_perm, ptile, axis=1),
                   'bonferroni' : np.percentile(PL_perm, bonf_ptile, axis=1),
                   'holmes' : np.percentile(PL_perm_holmes, ptile, axis=1)
                   }
    
    temp_ = mtlb.repmat(real_dat.reshape(real_dat.shape[0], 1), 1, nperms+1)

    exceed_thresh = (PL_perm >= temp_)*1
    exceed_holmesthresh = (PL_perm_holmes >= temp_)*1


    pval_uncorrected = exceed_thresh.sum(axis=1)/(nperms+1)
    pval_bonfcorrected = pval_uncorrected * len(phaselocksum)   
    pval_bonfcorrected[pval_bonfcorrected>1] = 1
    pval_holmescorrected = exceed_holmesthresh.sum(axis=1)/(nperms+1)

    
    pvals = {'uncorrected' : pval_uncorrected,
             'bonferroni' : pval_bonfcorrected,
             'holmes' : pval_holmescorrected
              }
    
    return thresh_sign, perm_dat, PL_perm, pvals


def plotter(freqs, PL_dat, thresh_sign, pvals, cmplx_vals,
            correction='bonferroni', 
            title='this is title'):
    
    plt.figure()
    plt.subplot(1, 3, (1,2))
    plt.plot(freqs, PL_dat, linewidth=3, label='empirical data')
    plt.plot(freqs, thresh_sign['uncorrected'], '--', c=(0.90,0.60,0.00), 
         label='uncorrected threshold',  linewidth=2)
    plt.plot(freqs, thresh_sign[correction], '--', c=(1.00,0.16,0.00), 
         label=correction + ' correction',  linewidth=2)

    plt.title(title)

    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude (a.u.)')
    plt.legend(loc=8)

    vect_pvals_corrected = pvals[correction]
    
    mask_correction = vect_pvals_corrected < .05
    mask_uncorrected = pvals['uncorrected'] < .05
    
    print(title)
    nosignificance = False
    sig_corrected = False
    
    if np.any(mask_correction):
        
        sig_corrected = True
        
        print('Peaks surviving ' + correction + ' correction')
        
        freqs_surviving = freqs[mask_correction]
        pval_surviving = vect_pvals_corrected[mask_correction]
        
        for ipeak in range(len(freqs_surviving)):
            
            print('Freq: ' + str(np.round(freqs_surviving[ipeak], 3)) + 
                  '; p=' + str(np.round(pval_surviving[ipeak], 3)))
            
    else:
        
        print('No peaks surviving ' + correction + ' correction')
            
    if np.any(mask_uncorrected):
        
        print('Peaks BEFORE ' + correction + ' correction')
        
        freqs_surviving = freqs[mask_uncorrected]
        pval_surviving = pvals['uncorrected'][mask_uncorrected]
        
        for ipeak in range(len(freqs_surviving)):
            
            print('Freq: ' + str(np.round(freqs_surviving[ipeak], 3)) + 
                  '; p=' + str(np.round(pval_surviving[ipeak], 3)))
           
    else:
        
        print('No peaks at all here :(')
        nosignificance = True
        
    # find highest peak to show vectors there    
    if nosignificance:
            
        mask_showphaselock = PL_dat.argmax()
    
    else:
        
        if sig_corrected:
            
            mask_showphaselock = vect_pvals_corrected.argmin()
            
        else:
        
            mask_showphaselock = pvals['uncorrected'].argmin()          
            
    
    max_freq = freqs[mask_showphaselock]
    cmplx_vectors_max = cmplx_vals[mask_showphaselock, :]
    lims = np.max(np.abs(cmplx_vectors_max))
    
    
    plt.subplot(1, 3, 3)
    plt.scatter(cmplx_vectors_max.real, cmplx_vectors_max.imag, s=15, c='k')
    plt.xlim((-lims, lims))
    plt.ylim((-lims, lims))
    plt.xlabel('real plane')
    plt.ylabel('complex plane')

    PLvect = cmplx_vectors_max.sum() / cmplx_vectors_max.shape
    plt.scatter(PLvect.real, PLvect.imag, s=20, c='r')
    plt.plot([0, PLvect.real], [0, PLvect.imag], c='r', linewidth=3)
    plt.title('complex vectors at ' + str(np.round(max_freq, 2)) + ' Hz')
    
    plt.tight_layout()

    

#%% load data matrix


path_to_dat = '../general analyses/' # change accordingly

# load data already detrended
mat_detection = loadmat(path_to_dat + 'det_face.mat')
det_face_task = mat_detection['det_face_task']
det_face_task = det_face_task[0, :, :]

#%% perform analyses and relative permutations

PL_dat_face, freqs, cmplx_face = cmpt_phase_locked(det_face_task)
thresh_sign_face, perm_dat_face, PL_perm_face, pvals_face = autoregress_model(det_face_task)

#%% 

plotter(freqs, PL_dat_face, thresh_sign_face, pvals_face, cmplx_face,
            correction='bonferroni', 
            title='Face detection')

