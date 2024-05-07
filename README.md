# General
The scripts contained in the project control an MLA3 from Intermodulation Products SE (https://www.intermod.pro/products/multifrequency-lock-in-amplifier).
To run the scripts, they need to be opened in the script panel of the IMP MLA software supplied by intermodulation products.
The purpose are fast QCM measurements (FastQCM.py) and fast EQCM measurements (FastEQCM.py). The output parameters in both scripts are the overtone-normalized complex frequency shifts, "Dfcbyns".

#  Fast QCM
The key feature of the MLA is the interrogation of resonances in one shot with a “comb”. The time resolution is equal to the inverse frequency spacing in the comb and is in the range of 1 ms for a resonator immersed in an aqueous solution (with the half bandwidth Gamm > 1000Hz). The number of frequencies ("nfreqs_few") in a comb usually is less than 32 (32 is the maximum number supported by the MLA). nfreqs_few follows from the frequency spacing ("df_comb") and the width of the resonance. For fitting the resonance curves with a phase-shifted Lorentzian, a minimum number of three is needed. "df_comb" is made as large as possible to achieve good time resolution (df_comb must be less than Gamm).  
For an initial search of the resonances "broad" combs are used. The frequency sweep panel in the IMP MLA software may also help when searching the resonances. For further information please see the further comments in the ReadMe.pdf and the python scripts.

# Fast EQCM
Please see the comments in Fast QCM. In addition, the script FastEQCM.py exploits potential modulation and averaging. Therfore, a trigger must be supplied to channel 1 (see the line: "trigsig = meta[3,:]", and the handbook).
In FastEQCM.py, the average over the respective modulation cycle is subtracted before accumulation (slow drifts are subtracted). The averages (the slow drifts) are saved and displayed, as well. Accumulation first occurs inside a "Loop". The duration of a Loop is limited by the storage space of the MLA, which amounts to 90 000 resonance curves (nDataSetsMax_MLA = 90 000). With a time resolution of 5 ms, the accumulation time can be 90 000 * 0.005 seconds = 450 seconds, at most. Accumulation also occurs between Loops with the line "Dfcbyns_ModResponse[:,iovt] = Dfcbyns_ModResponse[:,iovt]*(i_Loop)/(i_Loop+1) + Dfcbyns_accu[:,iovt] * 1/(i_Loop+1)".
For further information please see the further comments in the ReadMe.pdf and the python scripts.

# Some Literature on Fast QCM and Fast EQCM:
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
- Leppin, C.; Hampel, S.; Meyer, F. S.; Langhoff, A.; Fittschen, U. E. A.; Johannsmann, D. 
  A Quartz Crystal Microbalance, Which Tracks Four Overtones in Parallel with a Time Resolution of 10 Milliseconds: Application to Inkjet Printing. 
  Sensors 2020, 20 (20), 5915. https://doi.org/10.3390/s20205915.

- Wiegmann, J.; Leppin, C.; Langhoff, A.; Schwaderer, J.; Beuermann, S.; Johannsmann, D.; Weber, A. P. 
  Influence of the Solvent Evaporation Rate on the B-Phase Content of Electrosprayed PVDF Particles and Films Studied by a Fast Multi-Overtone QCM. 
  Advanced Powder Technology 2022, 33 (3), 103452. https://doi.org/10.1016/j.apt.2022.103452.

The following people have contributed: 
- Christian Leppin (christian.leppin(at)tu-clausthal.de)
- Frederick Sebastian Meyer
- Diethelm Johannsmann  
