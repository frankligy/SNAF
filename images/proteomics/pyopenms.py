#!/data/salomonis2/LabFiles/Frank-Li/proteomics_practice/pyopenms_env/bin/python3.7

import pyopenms
from pyopenms import *
from urllib.request import urlretrieve
import os,sys

'''
First tutorial: Get started
'''

# # # download small example file
# # gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-extra/master"
# # urlretrieve (gh + "/src/data/tiny.mzML", "tiny.mzML")

# exp = MSExperiment()
# MzMLFile().load("tiny.mzML", exp)
# # exp is subscriptable, each item is a spectrum
# for spectrum in exp: 
#     # spectrum is a pyopenms.pyopenms_1.MSSpectrum object, spectrum.size() return the number of peaks in this spectrum
#     # spectrum is iterable and subscriptable as well, each item will be pyopenms.pyopenms_4.Peak1D
#     mz, intensities = spectrum.get_peaks()
#     # mz is the ndarray of mz of peaks, intensities is the ndarray of intensity of peaks
#     ms_level = spectrum.getMSLevel()
#     # ms_level is the int of either 1 or 2
#     rt = spectrum.getRT()


'''
Second tutorial 2.1
'''
# spectrum = MSSpectrum()
# mz = range(1500, 500, -100)
# i = [0 for mass in mz]
# spectrum.set_peaks([mz, i])
# spectrum.sortByPosition()  # sort peaks based on mz in ascending order

'''
Second tutorial 2.2
'''
# # create spectrum and set properties
# spectrum = MSSpectrum()
# spectrum.setDriftTime(25) # 25 ms
# spectrum.setRT(205.2) # 205.2 s
# spectrum.setMSLevel(3) # MS3

# # Add peak(s) to spectrum
# spectrum.set_peaks( ([401.5], [900]) )

# # create precursor information
# p = Precursor()
# p.setMZ(600) # isolation at 600 +/- 1.5 Th
# p.setIsolationWindowLowerOffset(1.5)
# p.setIsolationWindowUpperOffset(1.5)
# p.setActivationEnergy(40) # 40 eV
# p.setCharge(4) # 4+ ion

# # and store precursor in spectrum
# spectrum.setPrecursors( [p] )

# # set additional instrument settings (e.g. scan polarity)
# InstrumentSettings = InstrumentSettings()
# InstrumentSettings.setPolarity(IonSource.Polarity.POSITIVE)
# spectrum.setInstrumentSettings(InstrumentSettings)

# # Optional: additional data arrays / peak annotations
# fda = FloatDataArray()
# fda.setName("Signal to Noise Array")
# fda.push_back(15)
# sda = StringDataArray()
# sda.setName("Peak annotation")
# sda.push_back("y15++")
# spectrum.setFloatDataArrays( [fda] )
# spectrum.setStringDataArrays( [sda] )

# # Add spectrum to MSExperiment
# exp = MSExperiment()
# exp.addSpectrum(spectrum)

# # Add second spectrum to the MSExperiment
# spectrum2 = MSSpectrum()
# spectrum2.set_peaks( ([1, 2], [1, 2]) )
# exp.addSpectrum(spectrum2)

# # store spectra in mzML file
# MzMLFile().store("testfile.mzML", exp)


'''
Second tutorial 2.3
'''
# import numpy as np

# def gaussian(x, mu, sig):
#     return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

# # Create new chromatogram
# chromatogram = MSChromatogram()

# # Set raw data (RT and intensity)
# rt = range(1500, 500, -100)
# i = [gaussian(rtime, 1000, 150) for rtime in rt]
# chromatogram.set_peaks([rt, i])

# # Sort the peaks according to ascending retention time
# chromatogram.sortByPosition()

# # Iterate over chromatogram of those peaks
# for p in chromatogram:
#     print(p.getRT(), p.getIntensity())

# # More efficient peak access with get_peaks()
# for rt, i in zip(*chromatogram.get_peaks()):
#     print(rt, i)

# # Access a peak by index
# print(chromatogram[2].getRT(), chromatogram[2].getIntensity())

# # Add meta information to the chromatogram
# chromatogram.setNativeID("Trace XIC@405.2")

# # Store a precursor ion for the chromatogram
# p = Precursor()
# p.setIsolationWindowLowerOffset(1.5)
# p.setIsolationWindowUpperOffset(1.5)
# p.setMZ(405.2) # isolation at 405.2 +/- 1.5 Th
# p.setActivationEnergy(40) # 40 eV
# p.setCharge(2) # 2+ ion
# p.setMetaValue("description", chromatogram.getNativeID())
# p.setMetaValue("peptide_sequence", chromatogram.getNativeID())
# chromatogram.setPrecursor(p)

# # Also store a product ion for the chromatogram (e.g. for SRM)
# p = Product()
# p.setMZ(603.4) # transition from 405.2 -> 603.4
# chromatogram.setProduct(p)

# # Store as mzML
# exp = MSExperiment()
# exp.addChromatogram(chromatogram)
# MzMLFile().store("testfile3.mzML", exp)

'''
Second tutorial 2.4
'''
# # The following examples creates an MSExperiment which holds six
# # MSSpectrum instances.
# exp = MSExperiment()
# for i in range(6):
#     spectrum = MSSpectrum()
#     spectrum.setRT(i)
#     spectrum.setMSLevel(1)
#     for mz in range(500, 900, 100):
#       peak = Peak1D()
#       peak.setIntensity(100 - 25*abs(i-2.5) )
#       spectrum.push_back(peak)
#     exp.addSpectrum(spectrum)

# # Iterate over spectra
# for spectrum in exp:
#     for peak in spectrum:
#         print (spectrum.getRT(), peak.getMZ(), peak.getIntensity())

# MzMLFile().store("testfile2.mzML", exp)

'''
Second tutorial 2.5
'''
# from urllib.request import urlretrieve
# gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-extra/master"
# urlretrieve (gh + "/src/data/PrecursorPurity_input.mzML", "PrecursorPurity_input.mzML")

# exp = MSExperiment()
# MzMLFile().load("PrecursorPurity_input.mzML", exp)
# # for this example, we check which are MS2 spectra and choose one of them
# for element in exp:
#     print(element.getMSLevel())
# # get the precursor information from the MS2 spectrum at index 3
# ms2_precursor = Precursor()
# ms2_precursor = exp[3].getPrecursors()[0]
# # get the previous recorded MS1 spectrum
# isMS1 = False
# i = 3 # start at the index of the MS2 spectrum
# while isMS1 == False:
#     if exp[i].getMSLevel() == 1:
#         isMS1 = True
#     else:
#         i -= 1
# ms1_spectrum = exp[i]
# # calculate the precursor purity in a 10 ppm precursor isolation window
# purity_score = PrecursorPurity().computePrecursorPurity(ms1_spectrum, ms2_precursor, 10, True)
# print(purity_score.total_intensity) # 9098343.890625
# print(purity_score.target_intensity) # 7057944.0
# print(purity_score.signal_proportion) # 0.7757394186070014
# print(purity_score.target_peak_count) # 1
# print(purity_score.residual_peak_count) # 4

'''
Second tutorial 2.6
'''
# from urllib.request import urlretrieve
# gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-extra/master"
# urlretrieve (gh + "/src/data/tiny.mzML", "test.mzML")
# inp = MSExperiment()
# MzMLFile().load("test.mzML", inp)
# # filter by MS level
# filtered = MSExperiment()
# for s in inp:
#     if s.getMSLevel() > 1:
#         filtered.addSpectrum(s)
# # filter by scan number
# scan_nrs = [0, 2, 5, 7]
# filtered = MSExperiment()
# for k, s in enumerate(inp):
#   if k in scan_nrs:
#     filtered.addSpectrum(s)
# # filter spctra and peaks
# mz_start = 6.0
# mz_end = 12.0
# filtered = MSExperiment()
# for s in inp:
#   if s.getMSLevel() > 1:
#     filtered_mz = []
#     filtered_int = []
#     for mz, i in zip(*s.get_peaks()):
#       if mz > mz_start and mz < mz_end:
#         filtered_mz.append(mz)
#         filtered_int.append(i)
#     s.set_peaks((filtered_mz, filtered_int))
#     filtered.addSpectrum(s)





