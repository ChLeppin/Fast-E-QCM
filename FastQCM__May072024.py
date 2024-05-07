'''
controls an MLA3 from Intermodulation Products AB 
(https://www.intermod.pro/products/multifrequency-lock-in-amplifier)
The purpose are fast QCM measurements.
Output parameters are the overtone-normalized complex frequency shifts, "Dfcbyns"
The number of frequencies ("nfreqs_few") in a comb usually is less than 32 
    (32 is the maximum number supported by the MLA)
    nfreqs_few follows from the frequency spacing ("df_comb_comb") and the width 
        of the resonance
    "df_comb" is made as large as possible to achieve good time resolution
    (df_comb must be less than Gamm) 

Literature:

- Leppin, C.; Hampel, S.; Meyer, F. S.; Langhoff, A.; Fittschen, U. E. A.; Johannsmann, D. 
  A Quartz Crystal Microbalance, Which Tracks Four Overtones in Parallel with a Time Resolution of 10 Milliseconds: Application to Inkjet Printing. 
  Sensors 2020, 20 (20), 5915. https://doi.org/10.3390/s20205915.

- Wiegmann, J.; Leppin, C.; Langhoff, A.; Schwaderer, J.; Beuermann, S.; Johannsmann, D.; Weber, A. P. 
  Influence of the Solvent Evaporation Rate on the B-Phase Content of Electrosprayed PVDF Particles and Films Studied by a Fast Multi-Overtone QCM. 
  Advanced Powder Technology 2022, 33 (3), 103452. https://doi.org/10.1016/j.apt.2022.103452.

The following people have contributed:
Christian Leppin (christian.leppin(at)tu-clausthal.de)
Frederick Sebastian Meyer
Diethelm Johannsmann
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
    Run the script. 
'''

from scipy.optimize import curve_fit
import datetime
import os
import numpy as np
import time

##################NEEDED TO BE SPECIFIED FOR MEASUREMENT#####################################
# 's' is for arrays (plural)
fress_ref = np.array([ 5009287.8, 14852885.4, 24751460.8, 34648650.6, 44571991.9, 54450400])
Gamms_ref = np.array([     120.3,      119.3,      134.9,      167.3,      210.3, 200     ])
fress     = np.array([ 5009287.8, 14852885.4, 24751460.8, 34648650.6, 44571991.9, 54450400])
Gamms     = np.array([      2000,      2000.,      2000.,      2000.,      2000., 2000    ])
# fress, Gamms : starting values, need to be determined before the measurment
measure_ovts = np.array([0, 1, 1, 1, 0, 0 ]) # include in measurement 
ovt_orders   = np.array([1, 3, 5, 7, 9, 11]) # ovt_orders differ from iovt (an index, 0,1,....)
n_ovt = len(fress)
Gmaxs     = np.zeros(n_ovt)
Phass     = np.zeros(n_ovt)
nfreqs_few_max = 32

colors = ['w','red','limegreen','deepskyblue','fuchsia','lightgrey']


Time_Resolution = 5e-3  #Please see comment below

'''
time resolution is 5 milliseconds, here Time_Resolution sets the minimum 
distance between the comb frequencies df_comb = 1 / Time_Resolution -> df_comb = 200 Hz  
df_comb must be less than the half bandwidth, Gamm.
'''

n_averages = int(1.01*Time_Resolution / 1.e-3) #please see comment below
'''
the distance between comb frequencies can be increased to 
account for broader resonances. The time resolution remains unchained. It is comparable to some onbord integration 
during lockin amplification
'''
Duration_of_Experiment = 1#s
Name_of_Experiment   = 'Test'
PlotResonances = True
SearchResonancesOnly = False
ThresholdFac_Noise_Discard = 200000;
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
            line = ax.plot(times,vals_ovt[iovt],'-',color = colors[iovt]) 
    ax.set_xlabel('time [sec]')
    if title == 'dfres/n': 
        ax.set_ylabel('frequency shift /n [Hz]')
    elif title == 'dGamm/n': 
        ax.set_ylabel('bandwidth shift /n [Hz]')
    elif title == 'chisq': 
         ax.set_ylabel('chi-square')
        
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
                          (Ys_meas_ovt.real[iovt,:nfreqs_few[iovt]]-Ys_meas_ovt.real[iovt,0])*1e3,\
                               'x',color = colors[iovt]) 
            line = ax.plot(fs_ovt[     iovt,:nfreqs_few[iovt]]-fs_ovt[     iovt,0],\
                          (Ys_fit_ovt.real[ iovt,:nfreqs_few[iovt]]-Ys_fit_ovt.real[ iovt,0])*1e3,\
                               '-',color = colors[iovt]) 
            line = ax.plot(fs_ovt[     iovt,:nfreqs_few[iovt]]-fs_ovt[     iovt,0],\
                          (Ys_meas_ovt.imag[iovt,:nfreqs_few[iovt]]-Ys_meas_ovt.imag[iovt,0])*1e3,\
                               'o',color = colors[iovt]) 
            line = ax.plot(fs_ovt[     iovt,:nfreqs_few[iovt]]-fs_ovt[     iovt,0],\
                              (Ys_fit_ovt.imag[ iovt,:nfreqs_few[iovt]]-Ys_fit_ovt.imag[ iovt,0])*1e3,\
                                   '--',color = colors[iovt]) 
    ax.set_xlabel('frequency'); ax.set_xlabel('admittance');  

def Plot_All():
    fig.clear(); wx.CallAfter(scriptplot.draw); 
    #Plot_Resonances(freqs_broad_all_ovts,Ys_broad_all_ovts,Ys_fit_broad_all_ovts ,'',331)
    Plot_Resonances(freqs_few_ovt,Ys_few_ovt,Ys_few_fit_ovt,'',332)
    Plot_versus_time(time, Dfcbyns.real,'dfres/n',334);
    Plot_versus_time(time, Dfcbyns.imag,'dGamm/n',335);
    Plot_versus_time(time, chsq         ,'chisq'  ,336);
    wx.CallAfter(scriptplot.draw); filename_gfx = folder + Name_of_Experiment  + ".png"; 
    wx.CallAfter(fig.savefig, filename_gfx)

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
    # The minimizer in curve_fit cannot handle complex numbers
    #print(freqs,Ys,iovt)
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
                               method = 'trf') # curve_fit is an imported routine
    '''
    if ResParsFit[2] < 0:  
        #if the prefactor is negative, the sign is flipped and the phase is flipped by 180 degree 
        ResParsFit[2] = -ResParsFit[2]
        ResParsFit[3] = (ResParsFit[3] + np.pi)%(2.*np.pi)
    '''
    Ys_fit = ResCurve(freqs,*ResParsFit)
    chi2 = np.sum(np.abs(Ys_fit-Ys)**2) / len(Ys)
    return ResParsFit,Ys_fit,chi2

def Acquire_Kinetics():
    '''
    Configures the instument to acquire multiple overtones a certain time frame.   
    ----------------
    Parameters: 
    Returns: 
    t (1d numpy array of floats): time
    chsq (2d numpy array of floats): chi square
    dfc_by_n (2d numpy array of complex): overtone normalizedcomplex frequencyy shift 1/n(frequency shift + 1j* bandwidth)
    fs_ovt (1d numpy array of floats): number of comb frequencies for every overtone
    YAdm_ovt (2d array of complex): complex admittances on different overtones
    YAdm_ovt_fit (2d array of complex): fit complex admittances on different overtones 
    '''
    n_Points = int(Duration_of_Experiment / Time_Resolution);
    if n_Points > 90000 : n_Points = 90000
    Dfcbyns = np.zeros((n_ovt,n_Points),dtype = complex)
    chi2s     = np.zeros((n_ovt,n_Points))
    time = np.arange(n_Points)*Time_Resolution
    fs_ovt, freqs_few_ovt  = Set_freqs_and_amps_mult_ovt();
    Ys_few_ovt     = np.zeros((n_ovt,nfreqs_few_max),dtype=complex) 
    Ys_few_fit_ovt = np.zeros((n_ovt,nfreqs_few_max),dtype=complex) 
    
    mla.lockin.wait_for_new_pixels(n_Points + 1); 
    pixelsall,meta = mla.lockin.get_pixels(n_Points, data_format='IQreal')
    print ("Data done. Wait for data processing...")
    for i in range(n_Points):   
        pixels = pixelsall[:,i]; pixI = -pixels[1:(len(pixels)):2]; pixQ = -pixels[0:(len(pixels)-1):2]
        Y_arr_all = np.squeeze(np.asarray(pixQ+1j*pixI)) / (peak_amp/np.sum(nfreqs_few))
        i_count = 0
        for iovt in range(n_ovt):
            if measure_ovts[iovt] == 1:
                f_arr = freqs_few_ovt[iovt,:nfreqs_few[iovt]]
                Y_arr = Y_arr_all[i_count:i_count+nfreqs_few[iovt]]
                ResParsFit, Yfit_arr, chi2s[iovt,i] = Fit_ResCurve(f_arr,Y_arr,iovt)
                if i == 0:
                    fress_ref = ResParsFit[0];
                    Gamms_ref = ResParsFit[1];
                    for j in range(nfreqs_few_max):
                        if j < nfreqs_few[iovt]:
                            freqs_few_ovt[iovt,j] = f_arr[j]
                            Ys_few_ovt [iovt,j] = Y_arr[j]; 
                            Ys_few_fit_ovt[iovt,j] = Yfit_arr[j]; 
                Dfcbyns[iovt,i] =  ((ResParsFit[0] - fress[iovt]) + 1j*(ResParsFit[1] - Gamms[iovt]) )/ ovt_orders[iovt]; 
                i_count += nfreqs_few[iovt]
    for iovt in range(n_ovt):
        if measure_ovts[iovt] == 1:
            Dfcbyns[iovt] = remove_outliers(Dfcbyns[iovt],chi2s[iovt]) 
    return time,chi2s,Dfcbyns,freqs_few_ovt,Ys_few_ovt,Ys_few_fit_ovt

def Set_freqs_and_amps_mult_ovt():
    '''
    Sets frequencies and amplitudes for the comb during measurement. Several overtones can be included.   
    ----------------
    Parameters: 
    f0 (float): center of the resonance curve
    nfreqs_few (integer): number of comb frequencies for the overtone
    df_cmb (integer): frequency spacing in the comb
    Returns: 
    f_arr (1d numpy array of floats): comb frequencies
    fs_ovt (1d numpy array of floats): number of comb frequencies
    '''   
    n_array = []
    print('df_comb',df_comb)
    for iovt in range(n_ovt):
        if measure_ovts[iovt] == 1:
            f0_tuned, df_tuned = mla.lockin.tune1(fress[iovt],df_comb); 
            print(df_tuned , 'df_tuned ')
            n_center = int(np.around(f0_tuned/df_tuned))
            n_array = np.append(n_array,n_center - n_averages*int(nfreqs_few[iovt]/2.) + n_averages * np.arange(nfreqs_few[iovt]))
    mla.lockin.set_frequencies_by_n_and_df(n_array,df_tuned, wait_for_effect=False) 
    f_arr = mla.lockin.get_frequencies(); 
    fs_ovt = np.nan * np.ones((n_ovt,nfreqs_few_max))
    amplitudes = np.zeros(nfreqs_broad)
    i_count = 0
    for iovt in range(n_ovt):
        if measure_ovts[iovt] == 1:
            for i in range(nfreqs_few[iovt]): 
                fs_ovt[iovt,i] = f_arr[i + i_count]
                amplitudes[i + i_count] = peak_amp/np.sum(nfreqs_few)
            i_count += nfreqs_few[iovt]
    mla.lockin.set_amplitudes(amplitudes)
    return f_arr, fs_ovt


    

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


def remove_outliers(in_array, chsq): 
    '''
    Discards bad fits based on chi^2
    ----------------
    Parameters: 
    in_array (1d array of float): data to be discarded (oftn Dfcbyns)
    chsq (1d array of float): corresponding chi square to the Name_of_Experiment al data
    ----------------
    Returns: 
    out_array (1d numpy array of floats): outlier corrected Name_of_Experiment al data
    ''' 
    out_array = in_array
    for i in range(len(in_array)):
        if chsq[i] > ThresholdFac_Noise_Discard * np.average(chsq) : out_array[i] = np.nan
    return out_array

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

def print_f():
    print     ("\n",\
              "1 f {:.1f} Hz".format(ResParsFit_all_ovts[0,0]) + ", G {:.1f} Hz".format(ResParsFit_all_ovts[0,1]),"\n",\
              "2 f {:.1f} Hz".format(ResParsFit_all_ovts[1,0]) + ", G {:.1f} Hz".format(ResParsFit_all_ovts[1,1]),"\n",\
              "3 f {:.1f} Hz".format(ResParsFit_all_ovts[2,0]) + ", G {:.1f} Hz".format(ResParsFit_all_ovts[2,1]),"\n",\
              "4 f {:.1f} Hz".format(ResParsFit_all_ovts[3,0]) + ", G {:.1f} Hz".format(ResParsFit_all_ovts[3,1]),"\n",\
              "5 f {:.1f} Hz".format(ResParsFit_all_ovts[4,0]) + ", G {:.1f} Hz".format(ResParsFit_all_ovts[4,1]),"\n",\
              "6 f {:.1f} Hz".format(ResParsFit_all_ovts[5,0]) + ", G {:.1f} Hz".format(ResParsFit_all_ovts[5,1]))
    if np.max(ResParsFit_all_ovts[:,1]) > 2*df_comb: 
        print('ATTENTION: For some overtone(s), the bandwidth is larger than 2*df_comb! \n The fit might fail.')
    
def SaveData():
    rowlabels = [np.array(["t[ms]"])]
    for iovt in range(n_ovt):
        if measure_ovts[iovt] == 1:
            rowlabels_iovt = [np.array(["dfres/n(" + str(int(fress[iovt]/1E6)) + " MHz) [Hz]" ,"dGamm/n(" + str(int(Gamms[iovt]/1E6)) + " MHz) [Hz]" ])]
            rowlabels = np.append(rowlabels, rowlabels_iovt)
    data = [time*1e3]
    for iovt in range(n_ovt):
        if measure_ovts[iovt] == 1:
            data = np.append(data,[Dfcbyns[iovt,:].real],axis = 0)
            data = np.append(data,[Dfcbyns[iovt,:].imag],axis = 0)
    file = open(folder + Name_of_Experiment  + "_dfc_by_n.txt", 'w'); 
    np.savetxt(file,[rowlabels], fmt='%s', delimiter="\t"); np.savetxt(file, np.transpose(data), delimiter="\t",  newline='\r');  file.close();
    labels = np.array(["f [Hz]", "G [Hz]"])
    data   = np.zeros((6,2))
    data[:, 0] = ResParsFit_all_ovts[:,0]
    data[:, 1] = ResParsFit_all_ovts[:,1]
    file = open(folder + Name_of_Experiment  + "_ResPars.txt", 'w'); 
    np.savetxt(file,[labels], fmt='%s', delimiter="\t"); np.savetxt(file, data, delimiter="\t",  newline='\r');  file.close();

############ Main Program  ########## 

date = datetime.datetime.now().date(); timef = time.strftime(" %H-%M-%S ")
folder = os.path.dirname(scriptutils.generate_filename()) + "/"+ str(date)+ " " +Name_of_Experiment +"/"
try:
    os.makedirs(folder)
except OSError:
    if not os.path.isdir(folder):
        raise



i_in_port = 1; nfreqs_broad = Initialize_MLA();
Sum_Gamms = 0; 
for iovt in range(n_ovt):
    if measure_ovts[iovt] == 1:
        Sum_Gamms += Gamms[iovt] 
df_comb= 1. / Time_Resolution 

nf_wid_fac = (nfreqs_broad * df_comb / Sum_Gamms)*1.01
nfreqs_few = measure_ovts * (nf_wid_fac*Gamms/df_comb).astype(int)
for iovt in range(n_ovt):
    if measure_ovts[iovt] == 1:
        if nfreqs_few[iovt] < 3 : nfreqs_few[iovt] = 3;
        if nfreqs_few[iovt] > nfreqs_few_max : nfreqs_few[iovt] = nfreqs_few_max;
nfreqs_few_sum = np.sum(nfreqs_few);
if nfreqs_few_sum > nfreqs_broad: 
    imax = np.argmax(nfreqs_few)
    nfreqs_few[imax] -= nfreqs_few_sum - nfreqs_broad

peak_amp = 1.95; 
freqs_broad_all_ovts  = np.ones((n_ovt,nfreqs_broad))*np.nan
Ys_broad_all_ovts     = np.ones((n_ovt,nfreqs_broad),dtype=complex)*np.nan
Ys_fit_broad_all_ovts = np.ones((n_ovt,nfreqs_broad),dtype=complex)*np.nan
ResParsFit_all_ovts   = np.zeros((n_ovt,6))

#Search resonances
for iovt in range(n_ovt):
    if measure_ovts[iovt] == 1:
        df_comb_ini = Gamms[iovt] / 5.; n_ini = 1
        ResParsFit_all_ovts[iovt], freqs_broad_all_ovts[iovt], Ys_broad_all_ovts[iovt], Ys_fit_broad_all_ovts[iovt] = Search_Resonance(fress[iovt], df_comb_ini, iovt, n_ini)
        fress[iovt] =  ResParsFit_all_ovts[iovt,0]
        Gamms[iovt] =  ResParsFit_all_ovts[iovt,1]
        Gmaxs[iovt] = ResParsFit_all_ovts[iovt,2]
        Phass[iovt] = ResParsFit_all_ovts[iovt,3]
print_f()

if not(SearchResonancesOnly):
    
    time, chsq, Dfcbyns,freqs_few_ovt,Ys_few_ovt,Ys_few_fit_ovt= Acquire_Kinetics()
    fig = scriptplot.fig; 
    SaveData()
    Plot_All()

mla.lockin.stop_lockin()
