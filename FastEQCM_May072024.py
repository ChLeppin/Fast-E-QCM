''' 
controls an MLA3 from Intermodulation Products AB 
(https://www.intermod.pro/products/multifrequency-lock-in-amplifier)
The purpose are fast EQCM measurements, exploiting potential modulation and averaging 
A trigger must be supplied to channel 1 (see the line: "trigsig = meta[3,:]", and the handbook)
Output parameters are the overtone-normalized complex frequency shifts, "Dfcbyns"
The average over the respective modulation cycle is subtracted before accumulation
    (slow drifts are subtracted)
    The averages (the slow drifts) are saved and displayed, as well 
The number of frequencies ("nfreqs_few") in a comb usually is less than 32 
    (32 is the maximum number supported by the MLA)
    nfreqs_few follows from the frequency spacing ("df_comb") and the width 
        of the resonance
    "df_comb" is made as large as possible to achieve good time resolution
    (df_comb must be less than Gamm) 

For an initial search of the resonances "broad" combs are used

Accumulation first occurs inside a "Loop". 
    The duration of a Loop is limited by the storage space of the MLA, which amounts 
    to 90 000 resonance curves (nDataSetsMax_MLA = 90 000).  With a time resolution of 
    5 ms, the accumulation time can be 90 000 * 0.005 seconds = 450 seconds, at most 
Accumulation also occurs between Loops with the line

Literature:

- Leppin, C.; Peschel, A.; Meyer, F. S.; Langhoff, A.; Johannsmann, D. 
  Kinetics of Viscoelasticity in the Electric Double Layer Following Steps in the 
  Electrode Potential Studied by a Fast Electrochemical Quartz Crystal Microbalance (EQCM). 
  Analyst 2021, 146 (7), 2160-2171. https://doi.org/10.1039/D0AN01965H.

- Leppin, C.; Langhoff, A.; Poggemann, H.-F.; Goedde, A. S.; Johannsmann, D. 
  Fast and Slow EQCM Response of Zwitterionic Weak Electrolytes to Changes in the 
  Electrode Potential: A pH-Mediated Mechanism. 
  Analyst 2021, 146 (19), 6005–6013. https://doi.org/10.1039/D1AN01306H.

- Leppin, C.; Langhoff, A.; Johannsmann, D. 
  Square-Wave Electrogravimetry Combined with Voltammetry Reveals Reversible 
  Submonolayer Adsorption of Redox-Active Ions. 
  Anal. Chem. 2022, 94 (28), 10227-10233. https://doi.org/10.1021/acs.analchem.2c01763.

- Leppin, C.; Langhoff, A.; Hoefft, O.; Johannsmann, D. A Modulation QCM Applied 
  to Copper Electrodeposition and Stripping. 
  Electroanalysis 2021, 33 (12), 2529-2538. https://doi.org/10.1002/elan.202100471.

The following people have contributed: 
 - Christian Leppin (christian.leppin(at)tu-clausthal.de)
 - Frederick Sebastian Meyer
 - Diethelm Johannsmann  
-------------------------------------
Procedures
Open the IMP MLA software supplied by intermodulation products and open this skript in the script panel.
1) Mount a resonator, set SearchResonancesOnly = True, 
    Go to Panels -> frequency sweep panel 
    Search the resonances in the dry state. 
    Enter the results into the lines for fress
2) Determine the reference state   
    Run the script and copy the output for fress and Gamms 
    into fress_ref and Gamms_ref. 
3) Add electrolyte onto the resonator and find the resonance again
    Go to Panels -> frequency sweep panel, 
    Search the resonances again. 
    Enter the results into the lines for fress
4) The measurement
    Set SearchResonancesOnly = False 
    Choose the time resolution 
    Start the potentiostat and set the variable "DataAcqRate_Potentiostat" and 
    "nPoints_per_ModCycle_Potentiostat"
    Adjust the trigger and define how many modulation cycles
    lie in the triggered interval (n_ModCycle_per_Trigger). 
    Run the script. 
 '''

import numpy as np
import time
import os
import datetime
from scipy.optimize import curve_fit

##################NEEDED TO BE SPECIFIED FOR MEASUREMENT#####################################
# 's' is for arrays (plural)
fress_ref = np.array([ 5009287.8, 14858217.5, 24765434.6, 34668387.2, 44571991.9, 54450400])
Gamms_ref = np.array([     120.3,      143.5,      156.0,      189.4,      210.3, 200     ])
fress     = np.array([ 5009287.8, 14870800.4, 24781711.7, 34691383.0, 44571991.9, 54450400])
Gamms     = np.array([      2000,      2000.,      2000.,      2000.,      2000., 2000    ])
# fress, Gamms : starting values, need to be determined before the measurment
measure_ovts = np.array([0, 1, 0, 0, 0, 0 ]) # include in measurement 
ovt_orders   = np.array([1, 3, 5, 7, 9, 11]) # ovt_orders differ from iovt (an index, 0,1,....)
n_ovt = len(fress)
Gmaxs     = np.zeros(n_ovt)
Phass     = np.zeros(n_ovt)
nfreqs_few_max = 32

colors = ['w','red','limegreen','deepskyblue','fuchsia','lightgrey']


Time_Resolution = 5e-3  #See comment below n_averages

'''
time resolution is 5 milliseconds, here 
Time_Resolution sets the minimum distance between the comb frequencies 
df_comb = 1 / Time_Resolution -> df_comb = 200 Hz  
df_comb must be less than the half bandwidth, Gamm.
'''

n_averages = int(1.01*Time_Resolution / 1.e-3) #please see comment below
'''
A bit of a tricky issue.  Because fitting may take longer than data aquisition 
on-board onbord integration of the admittance (before fitting) can be helpful.
An example: if n_averages is 2, the distance between the members of the combs is 
increased by 2.  In principle, the time resolution might inprove correspondingly, 
but the MLA in this case takes averages of 2 combs.  
The time resolution then stays the same, but the noise is improved
'''

nDataSetsMax_MLA = 90000
n_ModCycles_target = 10
n_Loops = 100  # n_Loops can be very large, measurement is then stopped from the keyboard
New_Folder_Interval   = 1 # changes the folder occasionally, so that intermediate results are not overwritten
DataAcqRate_Potentiostat = 1000.
nPoints_per_ModCycle_Potentiostat = 3000.
Period_Modulation = nPoints_per_ModCycle_Potentiostat/DataAcqRate_Potentiostat
n_Points_per_ModCycle = Period_Modulation / Time_Resolution

ThresholdFac_Noise_Discard = 10;  #results from failed fits do not enter accumulation

n_ModCycle_per_Trigger = 2
Name_of_Experiment = 'Test'
PlotResonances = True
SearchResonancesOnly = True
##########################################################################################

def Plot_versus_time(times,vals_ovt,title,location):
    ''' 
    plot values on the different overtones versus time 
    (part of the larger routine Plot_All)
    Parameters:
    ----------
    times (1d array of float): time
    vals_ovt (2d array of float):               
        1st index: time
        2nd index: overtone index
    title (string): title
    location (integer): a 3-digit-integer 
      1st digit: numper of rows
      2nd digit: numper of columns 
      3rd digit: where to place this particular plot (between 1 and total number of plots)
    '''
    ax = fig.add_subplot(location); ax.grid(False) 
    ax.tick_params(direction='in',right='on',top='on');ax.set_title(title)
    for iovt in range(n_ovt): 
        if measure_ovts[iovt] == 1:
            line = ax.plot(times,vals_ovt[:,iovt],'-',color = colors[iovt]) 

def Plot_Resonances(fs_ovt,Ys_meas_ovt,Ys_fit_ovt,title,location):
    '''
    part of the larger routine Plot_All
    Parameters: 
    ----------
    fs_ovt  (2d array of float)  : frequencies
    Ys_meas (2d array of complex) : admittances 
    Ys_fit  (2d array of complex) : fitted admittances
        1st index: overtone index
        2nd index: comb frequencies 
    Ys_meas and Ys_fit have units of inverse ohms (Siemens). 
        Values are multiplied by 1e3 because typical admittances are in the mS range
        Y is uncalibrated, though
    title and location: same as in Plot_versus_time  
    '''
    ax = fig.add_subplot(location); ax.grid(False) 
    ax.tick_params(direction='in',right='on',top='on')
    ax.set_title(title) 
    for iovt in range(n_ovt): 
        if measure_ovts[iovt] == 1:
            line = ax.plot(fs_ovt[     iovt,:nfreqs_few[iovt]]-fs_ovt[     iovt,0],\
                          (Ys_meas_ovt[iovt,:nfreqs_few[iovt]]-Ys_meas_ovt[iovt,0])*1e3,\
                               'x',color = colors[iovt]) 
            line = ax.plot(fs_ovt[     iovt,:nfreqs_few[iovt]]-fs_ovt[     iovt,0],\
                          (Ys_fit_ovt[ iovt,:nfreqs_few[iovt]]-Ys_fit_ovt[ iovt,0])*1e3,\
                               '-',color = colors[iovt]) 

def Plot_All():
    fig.clear(); wx.CallAfter(scriptplot.draw) 
    Plot_Resonances(freqs_few_ovt,Ys_few_ovt,Ys_few_fit_ovt,'Loop : '+str(int(i_Loop)),331)
    
    ax = fig.add_subplot(332); ax.grid(False) 
    ax.tick_params(direction='in',right='on',top='on');ax.set_title('triggers')
    xaxis_trig = 1.*np.arange(len(trigsig))/len(trigsig)
    line = ax.plot(xaxis_trig,trigsig, '-') 
    ax.set_title(Name_of_Experiment) 
    
    ax = fig.add_subplot(333); ax.grid(False) 
    ax.tick_params(direction='in',right='on',top='on');ax.set_title('fres,Gamm Slow')
    for iovt in range(n_ovt): 
        if measure_ovts[iovt] == 1: 
            line = ax.plot(times_slow,Dfcbyns_slow[:,iovt].real,'x-',color = colors[iovt])
            line = ax.plot(times_slow,Dfcbyns_slow[:,iovt].imag,'.-',color = colors[iovt]) 
    
    Plot_versus_time(times_in_ModCycle,Dfcbyns_ModResponse.real,'dfres/n',334)
    Plot_versus_time(times_in_ModCycle,Dfcbyns_ModResponse.imag,'dGamm/n',335)
    
    wx.CallAfter(scriptplot.draw); filename_gfx = folder + Name_of_Experiment + '.png' 
    wx.CallAfter(fig.savefig, filename_gfx)
    time.sleep(1.)

def ResCurve(freq,fres,Gamm,Gmax,Phas,OffG,OffB):
    ''' 
    phase-shifted Lorentzian (fit function)
    Parameters:
    ----------
    freq  (1d array of float)  : comb frequencies
    fres (float): resonance frequency
    Gamm (float): half bandwidth
    Gmax (float): prefactor (Amplitude)
    Phase (float): phase for rotation in the complex plance (assymmetry in resonance curve)
    OffG (float), OffB (float): real and imaginary offsets 
    Returns: 
    ----------
    complex admittance of the phase shifted Lorentzian
    '''
    return 1j* Gmax * np.exp(1j * Phas) *Gamm* 1. / \
        (fres+1j* Gamm - freq) +  OffG+1j*OffB

def Calc_GuessPars(freqs,Ys,iovt):
    ''' 
    creates an initial guess for the fit   
    Parameters:
    ----------
    freqs (1d array of float)   : comb frequencies 
    Ys    (1d array of complex) : complex admittances
    iovt  (integer)             : overtone iterator
    ---------
    Returns: 
    GuessPars (1d array of float): 
    guess for parameters of the resonance curve (as in ResCurve) 
    '''    
    Ys2 = Ys-1./2. * (Ys[0] + Ys[-1])
    fres = freqs[np.argmax(np.abs(Ys2))] 
    Gamm = np.abs(freqs[1]-freqs[2])
    Gmax = np.max(np.abs(Ys2))
    Phas    = np.angle(Ys[np.argmax(Ys)]-1./2. * (Ys[0] + Ys[-1]))
    OffG = np.average(Ys2.real)  
    OffB = np.average(Ys2.imag)
    GuessPars = [fres,Gamm,Gmax,Phas,OffG,OffB]
    return GuessPars

def Fit_ResCurve(freqs,Ys,iovt):
    ''' 
    performs the fit 
    Parameters: 
    ----------
    freqs (1d array of float)   : comb frequencies 
    Ys    (1d array of complex) : complex admittances
    iovt  (integer)             : overtone iterator
    Returns: 
    ---------
    ResParsFit (1d array of float) : fitted parmeters (as in ResCurve)
    Ys_fit   (1d array of complex) : complex admittances from the fit
    chi2                           : chi^2 (normalized sum of square deviations)
    '''    
    Ys_fit     = np.ndarray(len(Ys),dtype = complex)
    Ys_stacked = np.hstack((Ys.real,Ys.imag))
    freqs_stacked = np.hstack((freqs,freqs)) 
    # The minimizer in leastsq cannot handle complex numbers
    ResParsGuess = Calc_GuessPars(freqs,Ys,iovt) 
    ResParsBounds = ([np.min(freqs),10,          0     ,-2*np.pi      ,-np.inf,-np.inf],
                     [np.max(freqs),np.max(freqs)-np.min(freqs),np.inf,2*np.pi, np.inf, np.inf])

    def StackRealImag_ResCurve(freq,fres,Gamm,Gmax,Phase,OffG,OffB):
        '''
        ----------
        Parameters:
        freq (1d array of float): comb frequencies
        freq  (1d array of float)  : comb frequencies
        fres (float): resonance frequency
        Gamm (float): half bandwidth
        Gmax (float): prefactor (Amplitude)
        Phase (float): phase for rotation in the complex plance (assymmetry in resonance curve)
        OffG (float), OffB (float): real and imaginary offsets 
        ----------
        Returns:
        Stacked real and imaginary parts of the resonace curve (1d array of float)
        '''
        Ysmodel         = ResCurve(freqs,fres,Gamm,Gmax,Phase,OffG,OffB)  
        return np.hstack((Ysmodel.real, Ysmodel.imag))
    ResParsFit,cov = curve_fit(StackRealImag_ResCurve, 
                               xdata = freqs_stacked, 
                               ydata = Ys_stacked, 
                               p0 = ResParsGuess,
                               bounds = ResParsBounds,
                               method = 'dogbox') # curve_fit is an imported routine
    '''
    if ResParsFit[2] < 0:  
        #if the prefactor is negative, the sign is flipped and the phase is flipped by 180 degree 
        ResParsFit[2] = -ResParsFit[2]
        ResParsFit[3] = (ResParsFit[3] + np.pi)%(2.*np.pi)
    '''
    Ys_fit = ResCurve(freqs,*ResParsFit)
    chi2 = np.sum(np.abs(Ys_fit-Ys)**2) / len(Ys)
    return ResParsFit,Ys_fit,chi2

def Initialize_MLA():
    ''' 
    Returns: 
    nfreqs_broad (integer): maximum number of comb frequencies (often 32)
    '''   
    nfreqs_broad = mla.lockin.nr_output_freq; phases = np.zeros(nfreqs_broad) 
    out_1_mask = [1]*nfreqs_broad; mla.lockin.set_output_mask(out_1_mask,port=1) 
    out_2_mask = 0;                mla.lockin.set_output_mask(out_2_mask,port=2)
    mla.lockin.set_phases(phases,'degree') 
    mla.lockin.start_lockin(cluster_size=1)
    mla.lockin.set_input_multiplexer([i_in_port]*nfreqs_broad)
    return nfreqs_broad

def Set_freqs_and_amps_broad(fcenter,df_comb,nfreqs_broad):
    ''' 
    sets frequencies and amplitudes of a broad comb for an initial search of resonances
    Parameters: 
    ----------
    fcenter (float) : center of comb
    df_comb (float) : frequency spacing in the comb
    nfreqs_broad(integer)    : number of frequencies (often 32) 
    ---------
    Returns: 
    freqs (1d array of float): comb frequencies  
    '''    
    fcenter_tuned, df_comb_tuned = mla.lockin.tune1(fcenter,df_comb) 
        # tuning is critical, avoids Fourier leakage, see the manual of MLA
    n0       = int(np.around(fcenter_tuned/df_comb_tuned))
    ns_broad = (np.arange(nfreqs_broad)-np.around(nfreqs_broad/2.)) + n0 
    # ns_broad labels the frequencies inside the MLA
    mla.lockin.set_frequencies_by_n_and_df(ns_broad,df_comb_tuned,\
        wait_for_effect=False) 
    freqs_broad = mla.lockin.get_frequencies() 
    amplitudes = np.ones(nfreqs_broad)*peak_amp/nfreqs_broad
    mla.lockin.set_amplitudes(amplitudes)
    return freqs_broad

def Set_freqs_and_amps_few(fcenter,df_comb,nfreqs_few_ovt):
    ''' 
    Sets the frequencies and amplitudes for fast accumulation.
    Parameters: 
    ----------
    fcenter (float)          : center of comb
    df_comb (float)          : frequency spacing 
    nfreqs_few_ovt (integer) : number of comb frequencies for the overtone
    Returns: 
    ---------
    freqs_few (1d array of float): comb frequencies 
    '''    
    fcenter_tuned,df_comb_tuned = mla.lockin.tune1(fcenter,df_comb);
    n0 = int(np.around(fcenter_tuned/df_comb_tuned))
    print(df_comb_tuned,'df_comb_tuned')
    ns_few=n_averages*(np.arange(nfreqs_few_ovt)-np.around(nfreqs_few_ovt/2.))+n0
    mla.lockin.set_frequencies_by_n_and_df(ns_few,df_comb_tuned,\
            wait_for_effect=False) 
    freqs_few = mla.lockin.get_frequencies() 
    amplitudes = np.zeros(len(freqs_few))
    amplitudes[:nfreqs_few_ovt] = np.ones(nfreqs_few_ovt)*peak_amp/nfreqs_few_ovt
    mla.lockin.set_amplitudes(amplitudes)
    return freqs_few

def Search_Resonance(fcenter, df_comb_ini, iovt, n_avg_ini):
    '''  
    acquires a resonance using broad combs
    Parameters: 
    ---------
    fcenter (float)     : center of the resonance curve
    df_comb_ini (float) : frequency spacing in the comb 
    n_avg_ini (integer) : number of averages for acquiring the resonance
    ---------
    Returns: 
    ResParsFit (array of float):       : fit parameters
    freqs_braod (1d array of float)    : comb frequencies
    Ys_broad (1d array of complex)     : complex admittances
    Ys_broad_fit (1d array of complex) : complex admittances fit  
    '''
    freqs_broad = Set_freqs_and_amps_broad(fcenter,df_comb_ini,nfreqs_broad) 
    mla.lockin.wait_for_new_pixels(n_avg_ini + 1)
    if n_avg_ini == 1: 
        pixels,meta = mla.lockin.get_pixels(       n_avg_ini, data_format='IQreal')
    if n_avg_ini >  1: 
        pixels,meta = mla.lockin.get_pixel_average(n_avg_ini, data_format='IQreal')    
    pixI = -pixels[1:(len(pixels)):2]; pixQ = -pixels[0:(len(pixels)-1):2]
    Ys_broad = np.squeeze(np.asarray(pixQ+1j*pixI))
    ResParsFit, Ys_broad_fit,chi2 = Fit_ResCurve(freqs_broad, Ys_broad, iovt)
    fcenter = ResParsFit[0] 
    if PlotResonances:
        fig = scriptplot.fig; fig.clear(); wx.CallAfter(scriptplot.draw) 
        ax = fig.add_subplot(331); ax.grid(False)
        ax.tick_params(direction='in',right='on',top='on');ax.set_title(str(iovt))
        line = ax.plot(freqs_broad-freqs_broad[0], Ys_broad.real    , 'x')      
        line = ax.plot(freqs_broad-freqs_broad[0], Ys_broad_fit.real, '-') 
        wx.CallAfter(scriptplot.draw); time.sleep(0.5)
    return ResParsFit, freqs_broad, Ys_broad, Ys_broad_fit

def Acquire_Modulation_Response(iovt):
    ''' 
    Parameters:
    iovt (integer): overtone iterator
    ----------
    Returns: 
    trigsig (1d array of float): trigger signal (a set of steps)
    Dfcbyns_accu_ovt (1d array of complex): 
        overtone normalized frequency shift for one overtone 
    '''
    mla.lockin.set_trigger_flanks(trig1_flank='positive', trig2_flank='positive')
    freqs_few = Set_freqs_and_amps_few(ResParsFit_all_ovts[iovt,0],df_comb, nfreqs_few[iovt]) 
    # measure
    mla.lockin.wait_for_new_pixels(n_Points_per_Loop) 
    pixelsall,meta = mla.lockin.get_pixels(n_Points_per_Loop, data_format='IQreal')    
    trigsig = meta[3,:]
    i_Trigs = np.zeros(0, dtype = int)	# finds data points, where the trigger signal jumps
    for i in range(2,len(trigsig)-1):
        if (trigsig[i]-trigsig[i-1] > 0.5) and \
           (trigsig[i]-trigsig[i-1] > trigsig[i-1]-trigsig[i-2]) and \
           (trigsig[i]-trigsig[i-1] > trigsig[i+1]-trigsig[i]) :
            i_Trigs = np.append(i_Trigs,i)
    print('i_Trigs',i_Trigs)

    # retrieve and analyze
    Dfcbyns_Loop = np.zeros(n_Points_per_Loop,dtype = complex)
    chi2s_Loop   = np.zeros(n_Points_per_Loop)
    for i in range(n_Points_per_Loop):   
        pixels = pixelsall[:,i]; 
        pixI =-pixels[1:(len(pixels)  ):2]; 
        pixQ =-pixels[0:(len(pixels)-1):2]
        Ys_few = np.squeeze(np.asarray(pixQ+1j*pixI)) / (peak_amp / nfreqs_few[iovt])
        ResParsFit,Ys_few_fit,chi2s_Loop[i] = Fit_ResCurve(freqs_few[:nfreqs_few[iovt]],Ys_few[:nfreqs_few[iovt]],iovt)
        Dfcbyns_Loop[i] = ((ResParsFit[0]-fress_ref[iovt])+\
                        1j*(ResParsFit[1]-Gamms_ref[iovt]))/ovt_orders[iovt]
    Dfcbyns_accu, Dfcbyn_slow  = Accumulate(Dfcbyns_Loop,chi2s_Loop,i_Trigs,n_Points_per_ModCycle)
    for j in range(nfreqs_few[iovt]):
        freqs_few_ovt[iovt,j]  = freqs_few[j]
        Ys_few_ovt[iovt,j]     = Ys_few[j]
        Ys_few_fit_ovt[iovt,j] = Ys_few_fit[j]
    return trigsig, Dfcbyns_accu, Dfcbyn_slow

def Accumulate(vals,chi2s,i_Trigs,n_Points_per_ModCycle):
    ''' 
    Accumulates fres and Gamm in one Loop based on trigger events 
    Discards bad fits based on chi^2
    ----------
    Parameters:
    vals (1d array of complex)       : data to be accumulated (often Dfcbyns)
    chi2s (1d array of float)        : chi^2's 
    i_Trigs (1d array)               : data points corresponding to the trigger events
    n_ModCycle_per_Trigger (integer) : number of points per cycle 
    ----------
    Returns: 
    vals_accu (1d array of complex)  : accumulated data      
    val_avg   (float)                : average    
    '''
    vals_accu    = np.zeros(n_Points_per_ModCycle, dtype = complex)
    norm         = np.zeros(n_Points_per_ModCycle)
    chi2s_median = np.median(chi2s)
    for i in range(len(i_Trigs)):
        for j in range(1,n_ModCycle_per_Trigger):
            if len(vals) > i_Trigs[i] + (j+1)*n_Points_per_ModCycle :
                for k in range(n_Points_per_ModCycle):
                    if chi2s[i_Trigs[i] + j*n_Points_per_ModCycle + k] < \
                        chi2s_median * ThresholdFac_Noise_Discard : 
                        vals_accu[k] += vals[i_Trigs[i]+j*n_Points_per_ModCycle + k]
                        norm[k] += 1
    for k in range(n_Points_per_ModCycle):
        if norm[k] > 0 : vals_accu[k] /= norm[k]
    val_avg = np.nanmean(vals_accu)
    vals_accu -= val_avg 
    return vals_accu, val_avg

def Save_Modulation_Response(times_in_ModCycle,vals_real,vals_imag,label_real,label_imag,filename):
    ''' 
    Saves real and imaginary parts of some values (such as the complex frequency shift) 
    ----------
    Parameters:
    times_in_ ModCycle(2d array of float)
    vals_real (2d array of float): real part - 1st index: overtone - 2nd index: time
    vals_imag (2d array of float): imaginary part of complex array
    label_real (string): label for real part
    label_imag (stting): label for imagnary part
    filename (string): label for file name      
    '''
    rowlabels = [np.array(['time [ms]'])]
    for iovt in range(n_ovt):
        if measure_ovts[iovt] == 1:
            rowlabels_iovt = [np.array([label_real + '(' + str(int(fress[iovt]/1e6)) + ' MHz) [Hz]' ,\
                                        label_imag + '(' + str(int(fress[iovt]/1e6)) + ' MHz) [Hz]' ])]
            rowlabels = np.append(rowlabels, rowlabels_iovt)
    data = [times_in_ModCycle*1e3] # output in milliseconds
    for iovt in range(n_ovt):
        if measure_ovts[iovt] == 1:
            data = np.append(data,[vals_real[:,iovt] - np.average(vals_real[:,iovt])],axis = 0)
            data = np.append(data,[vals_imag[:,iovt] - np.average(vals_imag[:,iovt])],axis = 0)
    file = open(folder + Name_of_Experiment + filename + '.txt', 'w') 
    np.savetxt(file, [rowlabels], fmt='%s', delimiter='\t') 
    np.savetxt(file,np.transpose(data),fmt='%1.5e',delimiter='\t',newline='\r');file.close()

def Save_Slow(times_slow,Dfcbyns_slow,filename):
    rowlabels = [np.array(['time [s]'])]
    for iovt in range(n_ovt):
        if measure_ovts[iovt] == 1:
            rowlabels_iovt = [np.array(['Dfbyn(' + str(int(fress[iovt]/1e6)) + ' MHz) [Hz]' ,\
                                        'DGbyn(' + str(int(Gamms[iovt]/1e6)) + ' MHz) [Hz]' ])]
            rowlabels = np.append(rowlabels, rowlabels_iovt)
    data = [times_slow] # output in milliseconds
    for iovt in range(n_ovt):
        if measure_ovts[iovt] == 1:
            data = np.append(data,[Dfcbyns_slow[:,iovt].real],axis = 0)
            data = np.append(data,[Dfcbyns_slow[:,iovt].imag],axis = 0)
    file = open(folder + Name_of_Experiment + filename + '.txt', 'w') 
    np.savetxt(file,[rowlabels], fmt='%s', delimiter='\t') 
    np.savetxt(file,np.transpose(data),fmt='%1.5e',delimiter='\t',newline='\r');file.close()

def Save_Auxillary_Data(): 
    file = open(folder + Name_of_Experiment + "AuxillaryData.txt", 'w') 
    for i_ovt in range(n_ovt):
        if measure_ovts[i_ovt] == 1:
            file.write("i_ovt " + str(i_ovt) + '\n')
            file.write("fres " + str(ResParsFit_all_ovts[i_ovt,0]) + '\n')
            file.write("Gamm " + str(ResParsFit_all_ovts[i_ovt,1]) + '\n') 
            file.write("Gmax " + str(ResParsFit_all_ovts[i_ovt,2]) + '\n') 
            file.write("Phase "+ str(ResParsFit_all_ovts[i_ovt,3]) + '\n') 
            file.write("OffG " + str(ResParsFit_all_ovts[i_ovt,4]) + '\n') 
            file.write("OffB " + str(ResParsFit_all_ovts[i_ovt,5]) + '\n') 
            file.write("fres_ref " + str(fress_ref[i_ovt]) + '\n')
            file.write("Gamm_ref " + str(Gamms_ref[i_ovt]) + '\n')
    file.close()

def print_f():
    print     ('\n',\
              '1 f {:.1f} Hz'.format(ResParsFit_all_ovts[0,0]) + ', G {:.1f} Hz'.format(ResParsFit_all_ovts[0,1]),'\n',\
              '2 f {:.1f} Hz'.format(ResParsFit_all_ovts[1,0]) + ', G {:.1f} Hz'.format(ResParsFit_all_ovts[1,1]),'\n',\
              '3 f {:.1f} Hz'.format(ResParsFit_all_ovts[2,0]) + ', G {:.1f} Hz'.format(ResParsFit_all_ovts[2,1]),'\n',\
              '4 f {:.1f} Hz'.format(ResParsFit_all_ovts[3,0]) + ', G {:.1f} Hz'.format(ResParsFit_all_ovts[3,1]),'\n',\
              '5 f {:.1f} Hz'.format(ResParsFit_all_ovts[4,0]) + ', G {:.1f} Hz'.format(ResParsFit_all_ovts[4,1]),'\n',\
              '6 f {:.1f} Hz'.format(ResParsFit_all_ovts[5,0]) + ', G {:.1f} Hz'.format(ResParsFit_all_ovts[5,1]))
    if np.max(ResParsFit_all_ovts[:,1]) > 2*df_comb: 
        print('ATTENTION: For some overtone(s), the bandwidth is larger than 2*df_comb! \n The fit might fail.')
def Calc_rms_noise(x):
    '''
    Calculates a drift-corrected root-mean-square noise from the Hadamard variance
    ----------
    Parameters:
    x (1d array of float or complex): often overtone-normalized frequency shift
    ----------
    Returns:
    rms_noise (float)
    '''
    HadamardDeviations = np.zeros(len(x)-2,dtype = complex) 
    for i in range(1,len(x)-1):
        HadamardDeviations[i-1] = 1./6. * (x[i-1] - 2*x[i] + x[i+1])**2
    rms_noise = np.average(np.abs(HadamardDeviations))**0.5    
    return rms_noise
    
############ Main Program  ########## 

i_in_port = 1; nfreqs_broad = Initialize_MLA()
date = datetime.datetime.now().date(); timef = time.strftime('_%H-%M-%S_')
folder = os.path.dirname(scriptutils.generate_filename()) + \
    '/'+ str(date)+ timef+ '_' +Name_of_Experiment+'/'
try: os.makedirs(folder)
except OSError:
    if not os.path.isdir(folder): raise

'''
peak_amp can be choosen between 0 - 12 V on the MLA3 depending
on the selected output range (2 V or 12 V). Pay attention to the sensitivity 
of the input channel when changing the output range to 12 V.For most purpose, 
peak_amp = 1.95 is an appropriate choice combined with the sensitivity of 2.5 V 
on the input channel. For further information, please see the manual.
'''
peak_amp = 1.95  #MLA3 SEE THE COMMENT ABOVE
  
freqs_broad_all_ovts  = np.ones((n_ovt,nfreqs_broad))*np.nan
Ys_broad_all_ovts     = np.ones((n_ovt,nfreqs_broad),dtype=complex)*np.nan
Ys_fit_broad_all_ovts = np.ones((n_ovt,nfreqs_broad),dtype=complex)*np.nan
ResParsFit_all_ovts   = np.zeros((n_ovt,6))

# Search resonances 
for iovt in range(n_ovt): 
    if measure_ovts[iovt] == 1:
        df_comb_ini = Gamms[iovt] / 5; n_ini = 1  
        ResParsFit_all_ovts[iovt], freqs_broad_all_ovts[iovt], \
            Ys_broad_all_ovts[iovt], Ys_fit_broad_all_ovts[iovt] = \
            Search_Resonance(fress[iovt], df_comb_ini, iovt, n_ini)
        print ('ResParsFit_all_ovts[iovt]', ResParsFit_all_ovts[iovt])
        fress[iovt] = ResParsFit_all_ovts[iovt,0]
        Gamms[iovt] = ResParsFit_all_ovts[iovt,1]
        Gmaxs[iovt] = ResParsFit_all_ovts[iovt,2]
        Phass[iovt] = ResParsFit_all_ovts[iovt,3]

print_f()
if not(SearchResonancesOnly):
    # Initializes averages
    times_slow     = [0]
    time_slow_ini  = time.time()
    Dfcbyns_slow   = [((fress-fress_ref) + 1j*(Gamms-Gamms_ref)) / ovt_orders]

    # Initialize modulation response
    df_comb = 1 / Time_Resolution 
    n_Points_per_ModCycle = int(Period_Modulation / Time_Resolution)

    if n_ModCycles_target * n_Points_per_ModCycle < 90000 : 
        n_ModCycles = n_ModCycles_target
    else : 
        n_ModCycles = int(nDataSetsMax_MLA / n_Points_per_ModCycle)
        print ('n_ModCycles was too large, was lowered to ',n_ModCycles)

    n_Points_per_Loop = int(n_ModCycles * (Period_Modulation / Time_Resolution))
    times_in_ModCycle = 1. / df_comb * np.arange(n_Points_per_ModCycle)
    Dfcbyns_ModResponse      = np.zeros((n_Points_per_ModCycle,n_ovt),dtype=complex) 
    Dfcbyns_accu             = np.zeros((n_Points_per_ModCycle,n_ovt),dtype=complex) 
    nfreqs_few = np.zeros(n_ovt,dtype = int)
    for iovt in range(n_ovt) : 
        nfreqs_few[iovt] = 2 * int(3.* Gamms[iovt] / df_comb) + 1
        if nfreqs_few[iovt] < 3 : nfreqs_few[iovt] = 3
        if nfreqs_few[iovt] > nfreqs_few_max : nfreqs_few[iovt] = nfreqs_few_max

    freqs_few_ovt  = np.zeros((n_ovt,nfreqs_few_max)) 
    Ys_few_ovt     = np.zeros((n_ovt,nfreqs_few_max),dtype=complex) 
    Ys_few_fit_ovt = np.zeros((n_ovt,nfreqs_few_max),dtype=complex) 
    rms_noise_all_ovts     = np.zeros((n_ovt))
    rms_noise_all_ovts_raw = np.zeros((n_ovt))

    # long-time measurement
    i_Loop = 0 
    while i_Loop < n_Loops:
        Dfcbyns_avg  = np.zeros(n_ovt, dtype = complex)
        for iovt in range(n_ovt):
            if measure_ovts[iovt] == 1:
                trigsig,Dfcbyns_accu[:,iovt],Dfcbyns_avg[iovt] = Acquire_Modulation_Response(iovt)
                Dfcbyns_ModResponse[:,iovt] = Dfcbyns_ModResponse[:,iovt]*(i_Loop)/(i_Loop+1) + \
                                            Dfcbyns_accu[:,iovt]             *       1/(i_Loop+1)
                rms_noise_all_ovts[iovt]     = Calc_rms_noise(Dfcbyns_ModResponse[:,iovt])
                rms_noise_all_ovts_raw[iovt] = \
                    rms_noise_all_ovts[iovt] * (n_ModCycles * (i_Loop + 1))**0.5
                print ('rms_noise[iovt] {:.3f} Hz'.format(rms_noise_all_ovts[iovt]))
        Dfcbyns_slow = np.append(Dfcbyns_slow,[Dfcbyns_avg],axis=0)
        times_slow   = np.append(times_slow,time.time()-time_slow_ini)

        if (i_Loop % New_Folder_Interval == 0):
            date = datetime.datetime.now().date(); timef = time.strftime('_%H-%M-%S_')
            folder = os.path.dirname(scriptutils.generate_filename()) + \
                '/'+ str(date)+ timef+ '_' +Name_of_Experiment+'/'
            try: os.makedirs(folder)
            except OSError:
                if not os.path.isdir(folder): raise
        fig = scriptplot.fig 
        Plot_All(); 
        Save_Modulation_Response(times_in_ModCycle,Dfcbyns_accu.real,\
                                 Dfcbyns_accu.imag,'Df/n','DG/n','SingleCycle')
        Save_Modulation_Response(times_in_ModCycle,Dfcbyns_ModResponse.real,\
                                 Dfcbyns_ModResponse.imag,'Df/n','DG/n','ModResponse')
        Save_Slow(times_slow,Dfcbyns_slow,'SlowResponse')
        Save_Auxillary_Data()
        i_Loop += 1 
mla.lockin.stop_lockin()