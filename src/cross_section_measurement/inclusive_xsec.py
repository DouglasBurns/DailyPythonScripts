from tools.file_utilities import read_data_from_JSON, make_folder_if_not_exists
from math import sqrt
import os.path as check
from ROOT import gROOT, gPad, gStyle, TFile, TMath, TGraph, TH1F, TH1

#  ufloat(value, std.dev) can just add ufloat together and errors work out.
from uncertainties import ufloat



def calcIncXSec(inputFile, diff_xsection, Type, Channel, Measurement, Variable):
	
	print "-"*30
	print "Channel : ", Channel, ", Measurement : ", Measurement, ", Variable : ", Variable

	METCorrections = '_patType1CorrectedPFMet'
	if (Variable == 'HT'): METCorrections = ''

	phaseSpaceInfo_electron_Hist = inputFile.Get( 'unfolding_' + Variable + '_analyser_electron_channel' + METCorrections + '/phaseSpaceInfoHist' )
	phaseSpaceInfo_muon_Hist = inputFile.Get( 'unfolding_' + Variable + '_analyser_muon_channel' + METCorrections + '/phaseSpaceInfoHist' )

	if (Channel == "electron"): 
		# BR = 0.1444 # (0.1071*0.6741)*2 from w decay percentages in pdg
		BR = 0.146
		e_efficiency = phaseSpaceInfo_electron_Hist.GetBinContent(4)
		e_correction = phaseSpaceInfo_electron_Hist.GetBinContent(5)

	if (Channel == "muon"):  
		# BR = 0.1433 
		BR = 0.146
		mu_efficiency = phaseSpaceInfo_muon_Hist.GetBinContent(4)
		mu_correction = phaseSpaceInfo_muon_Hist.GetBinContent(5)


	if (Channel == "combined"): 
		# BR = 0.2877
		BR = 0.292
		e_efficiency = phaseSpaceInfo_electron_Hist.GetBinContent(4)
		e_correction = phaseSpaceInfo_electron_Hist.GetBinContent(5)
		mu_efficiency = phaseSpaceInfo_muon_Hist.GetBinContent(4)
		mu_correction = phaseSpaceInfo_muon_Hist.GetBinContent(5)
		print "Total E Eff: ", e_efficiency
		print "Total E Cor: ", e_correction
		print "Total Mu Eff: ", mu_efficiency
		print "Total Mu Cor: ", mu_correction

	# efficiencies
	# Electron channel : 0.0895435238802
	# Muon channel : 0.105790628035
	e_efficiency = 0.0895435238802 #change this back to graphical input
	mu_efficiency = 0.105790628035

	# N_offline_SL / N_SL "Efficiency" : (NEvents = XSec*Lumi*Br*EFF)
	# number of semileptonic events that pass the event selection divided by the total number of semileptonic events.

	# N_offline_SL / N_offline_all "Correction" : (NDataEvents scaled using MC to more appropriate N)
	# number of semileptonic events that pass the selection divided by the total number of events that pass the selection.

	# Combined : 1/Lumi*2Br.(combined electron scalings + combined muon scalings)


	# lumi = 40.028
	# Better Calibration
	lumi = 41.629

	# Change lumi by 10% if lumi systematic
	if (Measurement == 'luminosity+'): lumi = 1.1*lumi
	if (Measurement == 'luminosity-'): lumi = 0.9*lumi

	no_events = no_events_error = ufloat(0,0)
	inc_xsec = inc_xsec_error = ufloat(0,0)

	# Actually adds up number of events
	for values in diff_xsection[Type]:
		# print (values)
		tmp = ufloat(values[0],values[1])
		no_events += tmp	


	print "Number of Measured Events", no_events.nominal_value
	# print "Number of Measured Events Error", no_events.std_dev

	# Applies efficiency and scaling to the number events, then finds luminosity
	if (Channel == "electron"):
		print "e correction : ",  e_correction
		print "e_eff : ", e_efficiency
		print "correction : ", e_correction/e_efficiency
		no_events *= (e_correction/e_efficiency)
		print "Number of Measured Events after corrections", no_events.nominal_value
		# print "Number of Measured Events after corrections Error", no_events.std_dev
		inc_xsec = no_events/(lumi*BR)
		print "crossection", inc_xsec.nominal_value
	if (Channel == "muon"): 
		no_events *= (mu_correction/mu_efficiency)
		# print "Number of Measured Events after corrections", no_events.nominal_value
		# print "Number of Measured Events after corrections Error", no_events.std_dev
		inc_xsec = no_events/(lumi*BR)
		print "crossection", inc_xsec.nominal_value
	if (Channel == "combined"): 
		no_events = (no_events*e_correction/e_efficiency)+(no_events*mu_correction/mu_efficiency) #is this correct???
		# print "Number of Measured Events after corrections", no_events.nominal_value
		# print "Number of Measured Events after corrections Error", no_events.std_dev
		inc_xsec = no_events/(lumi*BR)
		print "crossection", inc_xsec.nominal_value

	return (inc_xsec)


if __name__ == '__main__':

	# File containing efficinecies and scaling factor
	inputFile = TFile('/hdfs/TopQuarkGroup/run2/unfolding/13TeV/50ns/unfolding_TTJets_13TeV.root')
	
	# Files to output the calculated inclusive cross section for each variable and a list of the systematics (% of central)
	outputFile1 = open('Inc_XSec.txt', 'w')
	outputFile2 = open('Systematics.txt','w')

	Energy = '13TeV'

	# Lists of Variables, Phase Spacii, Channels and Measurements (Systematics)
	Variable = {
	1 : 'HT', 
	2 : 'MET',
	3 : 'WPT',
	4 : 'ST',
	5 : 'NJets',
	6 : 'abs_lepton_eta',
	7 : 'lepton_pt',
	}

	PhaseSpace = {
	1 : 'FullPS',
	# 2 : 'VisiblePS',
	}

	Channel = {
	1 : 'electron',
	2 : 'muon',
	3 : 'combined',
	}

	Measurement = {
	1 : 'central',
	2 : 'TTJets_scaleup',
	3 : 'TTJets_scaledown',
	4 : 'TTJets_massup',
	5 : 'TTJets_massdown',
	6 : 'TTJets_powheg_herwig',
	7 : 'TTJets_powheg_pythia',
	8 : 'TTJets_hadronisation',
	9 : 'TTJet_cross_section+',
	11 : 'TTJet_cross_section-',
	12 : 'SingleTop_cross_section+',
	13 : 'SingleTop_cross_section-',
	14 : 'V+Jets_cross_section+',
	15 : 'V+Jets_cross_section-',
	16 : 'QCD_cross_section+',
	17 : 'QCD_cross_section-',
	18 : 'QCD_shape',
	19 : 'JES_up',
	20 : 'JES_down',
	21 : 'JER_up',
	22 : 'JER_down',
	23 : 'Electron_up',
	24 : 'Electron_down',
	25 : 'Muon_up',
	26 : 'Muon_down',
	27 : 'luminosity+',
	28 : 'luminosity-',
	29 : 'ElectronEnUp',
	30 : 'ElectronEnDown',
	31 : 'MuonEnUp',
	32 : 'MuonEnDown',
	33 : 'TauEnUp',
	34 : 'TauEnDown',
	35 : 'UnclusteredEnUp',
	36 : 'UnclusteredEnDown',
	37 : 'TTJets_NLOgenerator',
	38 : 'BJet_up',
	39 : 'BJet_down',
	}


	for ps in PhaseSpace :
		print >> outputFile2, "Phase Space : ", PhaseSpace[ps]
		for ch in Channel :
			print >> outputFile2, "Channel : ", Channel[ch]
			for var in Variable :
				print >> outputFile2, "Variable : ", Variable[var]
				for mes in Measurement :

# No variation - lepton eta broken 
					filename = (
						'/hdfs/TopQuarkGroup/run2/dpsData/data/normalisation/background_subtraction/'
						+ Energy + '/'
						+ Variable[var] + '/'
						+ PhaseSpace[ps] + '/' 
						+'xsection_measurement_results/' 
						+ Channel[ch] + '/' 
						+ Measurement[mes] + '/' 
						+ 'normalisation_patType1CorrectedPFMet.txt'
						)

# Small Vairation between all
					# filename = (
					# 	'DougsDataUsage/'
					# 	+ Energy + '/'
					# 	+ Variable[var] + '/'
					# 	+ PhaseSpace[ps] + '/' 
					# 	+'xsection_measurement_results/' 
					# 	+ Channel[ch] + '/' 
					# 	+ Measurement[mes] + '/' 
					# 	+ 'normalisation_patType1CorrectedPFMet.txt'
					# 	)
					# print (filename)

					# Excludes MET related files in HT (which dont exist)
					if not (check.isfile(filename)): continue;

					# Get the differential cross section measurements
					diff_xsection = read_data_from_JSON( filename )

					#  Calculating central value and all available systematics. Syst Err = Syst.nominal_value - Central.nominal_value
					if (Measurement[mes] == 'central'): 
						inc_mes_xsec = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						print >> outputFile2, "TTJet_measured Central    : ", Measurement[mes], " : %.1f " %(inc_mes_xsec.nominal_value), "%"
						# print (PhaseSpace[ps], Channel[ch], Variable[var])
						# print "Central : ", inc_data_xsec.nominal_value
						# print "Std Dev : ", inc_data_xsec.std_dev

					if (Measurement[mes] == 'luminosity+'):
						inc_mes_lumi_up = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_lumi_up = (inc_mes_lumi_up.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_lumi_up/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'luminosity-'):
						inc_mes_lumi_down = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_lumi_down = (inc_mes_lumi_down.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_lumi_down/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'TTJets_scaleup'):
						inc_mes_TTJets_scaleup = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TTJets_scaleup = (inc_mes_TTJets_scaleup.nominal_value - inc_mes_xsec.nominal_value)
						print inc_mes_xsec.nominal_value
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TTJets_scaleup/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'TTJets_scaledown'):
						inc_mes_TTJets_scaledown = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TTJets_scaledown = (inc_mes_TTJets_scaledown.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TTJets_scaledown/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'TTJets_massup'):
						inc_mes_TTJets_massup = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TTJets_massup = (inc_mes_TTJets_massup.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TTJets_massup/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'TTJets_massdown'):
						inc_mes_TTJets_massdown = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TTJets_massdown = (inc_mes_TTJets_massdown.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TTJets_massdown/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'TTJets_powheg_herwig'):
						inc_mes_TTJets_powheg_herwig = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TTJets_powheg_herwig = (inc_mes_TTJets_powheg_herwig.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TTJets_powheg_herwig/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'TTJets_powheg_pythia'):
						inc_mes_TTJets_powheg_pythia = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TTJets_powheg_pythia = (inc_mes_TTJets_powheg_pythia.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TTJets_powheg_pythia/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'TTJets_hadronisation'):
						inc_mes_TTJets_hadronisation = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TTJets_hadronisation = (inc_mes_TTJets_hadronisation.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TTJets_hadronisation/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'TTJet_cross_section+'):
						inc_mes_TTJets_cross_section_up = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TTJets_cross_section_up = (inc_mes_TTJets_cross_section_up.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TTJets_cross_section_up/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'TTJet_cross_section-'):
						inc_mes_TTJets_cross_section_down = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TTJets_cross_section_down = (inc_mes_TTJets_cross_section_down.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TTJets_cross_section_down/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'SingleTop_cross_section+'):
						inc_mes_SingleTop_cross_section_up = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_SingleTop_cross_section_up = (inc_mes_SingleTop_cross_section_up.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_SingleTop_cross_section_up/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'SingleTop_cross_section-'):
						inc_mes_SingleTop_cross_section_down = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_SingleTop_cross_section_down = (inc_mes_SingleTop_cross_section_down.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_SingleTop_cross_section_down/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'V+Jets_cross_section+'):
						inc_mes_V_Jets_cross_section_up = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_V_Jets_cross_section_up = (inc_mes_TTJets_cross_section_up.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_V_Jets_cross_section_up/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'V+Jets_cross_section-'):
						inc_mes_V_Jets_cross_section_down = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_V_Jets_cross_section_down = (inc_mes_TTJets_cross_section_down.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_V_Jets_cross_section_down/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'QCD_cross_section+'):
						inc_mes_QCD_cross_section_up = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_QCD_cross_section_up = (inc_mes_QCD_cross_section_up.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_QCD_cross_section_up/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'QCD_cross_section-'):
						inc_mes_QCD_cross_section_down = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_QCD_cross_section_down = (inc_mes_QCD_cross_section_down.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_QCD_cross_section_down/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'QCD_shape'):
						inc_mes_QCD_shape = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_QCD_shape = (inc_mes_QCD_shape.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_QCD_shape/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'JES_up'):
						inc_mes_JES_up = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_JES_up = (inc_mes_JES_up.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_JES_up/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'JES_down'):
						inc_mes_JES_down = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_JES_down = (inc_mes_JES_down.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_JES_down/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'JER_up'):
						inc_mes_JER_up = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_JER_up = (inc_mes_JER_up.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_JER_up/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'JER_down'):
						inc_mes_JER_down = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_JER_down = (inc_mes_JER_down.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_JER_down/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'Electron_up'):
						inc_mes_Electron_up = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_Electron_up = (inc_mes_Electron_up.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_Electron_up/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'Electron_down'):
						inc_mes_Electron_down = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_Electron_down = (inc_mes_Electron_down.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_Electron_down/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'Muon_up'):
						inc_mes_Muon_up = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_Muon_up = (inc_mes_Muon_up.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_Muon_up/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'Muon_down'):
						inc_mes_Muon_down = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_Muon_down = (inc_mes_Muon_down.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_Muon_down/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'ElectronEnUp'):
						inc_mes_ElectronEnUp = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_ElectronEnUp = (inc_mes_ElectronEnUp.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_ElectronEnUp/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'ElectronEnDown'):
						inc_mes_ElectronEnDown = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_ElectronEnDown = (inc_mes_ElectronEnDown.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_ElectronEnDown/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'MuonEnUp'):
						inc_mes_MuonEnUp = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_MuonEnUp = (inc_mes_MuonEnUp.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_MuonEnUp/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'MuonEnDown'):
						inc_mes_MuonEnDown = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_MuonEnDown = (inc_mes_MuonEnDown.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_MuonEnDown/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'TauEnUp'):
						inc_mes_TauEnUp = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TauEnUp = (inc_mes_TauEnUp.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TauEnUp/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'TauEnDown'):
						inc_mes_TauEnDown = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TauEnDown = (inc_mes_TauEnDown.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TauEnDown/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'UnclusteredEnUp'):
						inc_mes_UnclusteredEnUp = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_UnclusteredEnUp = (inc_mes_UnclusteredEnUp.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_UnclusteredEnUp/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'UnclusteredEnDown'):
						inc_mes_UnclusteredEnDown = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_UnclusteredEnDown = (inc_mes_UnclusteredEnDown.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_UnclusteredEnDown/inc_mes_xsec.nominal_value*100), "%"

					
					if (Measurement[mes] == 'BJet_up'):
						inc_mes_BJet_up = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_BJet_up = (inc_mes_BJet_up.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_BJet_up/inc_mes_xsec.nominal_value*100), "%"

					if (Measurement[mes] == 'BJet_down'):
						inc_mes_BJet_down = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_BJet_down = (inc_mes_BJet_down.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_BJet_down/inc_mes_xsec.nominal_value*100), "%"


					if (Measurement[mes] == 'TTJets_NLOgenerator'):
						inc_mes_TTJets_NLOgenerator = calcIncXSec( inputFile, diff_xsection, 'TTJet_measured', Channel[ch], Measurement[mes], Variable[var] )
						error_mes_TTJets_NLOgenerator = (inc_mes_TTJets_NLOgenerator.nominal_value - inc_mes_xsec.nominal_value)
						print >> outputFile2, "TTJet_measured Systematic : ", Measurement[mes], " : %.1f " %(error_mes_TTJets_NLOgenerator/inc_mes_xsec.nominal_value*100), "%"

				# TTJET MEASURED
				# Take the largest (up/down) error as the systematic error
				error_mes_TTJets_scale = max(abs(error_mes_TTJets_scaleup), abs(error_mes_TTJets_scaledown))
				error_mes_TTJets_mass = max(abs(error_mes_TTJets_massup), abs(error_mes_TTJets_massdown))
				error_mes_TTJets_generator = max(abs(error_mes_TTJets_powheg_herwig), abs(error_mes_TTJets_powheg_pythia))
				error_mes_TTJets_hadronisation = abs(error_mes_TTJets_hadronisation)
				error_mes_TTJets_NLOgenerator = abs(error_mes_TTJets_NLOgenerator)
				error_mes_TTJets_cross_section = max(abs(error_mes_TTJets_cross_section_up), abs(error_mes_TTJets_cross_section_down))
				error_mes_SingleTop_cross_section = max(abs(error_mes_SingleTop_cross_section_up), abs(error_mes_SingleTop_cross_section_down))
				error_mes_V_Jets_cross_section = max(abs(error_mes_V_Jets_cross_section_up), abs(error_mes_V_Jets_cross_section_down))
				error_mes_QCD_cross_section = max(abs(error_mes_QCD_cross_section_up), abs(error_mes_QCD_cross_section_down))
				error_mes_QCD_shape = abs(error_mes_QCD_shape)
				error_mes_BJet = max(abs(error_mes_BJet_down), abs(error_mes_BJet_up))
				error_mes_JER = max(abs(error_mes_JER_down), abs(error_mes_JER_up))
				error_mes_JES = max(abs(error_mes_JES_down), abs(error_mes_JES_up))
				error_mes_Electron = max(abs(error_mes_Electron_up), abs(error_mes_Electron_down))
				error_mes_Muon = max(abs(error_mes_Muon_up), abs(error_mes_Muon_down))
				if (Variable[var] != 'HT'):
					error_mes_ElectronEn = max(abs(error_mes_ElectronEnUp), abs(error_mes_ElectronEnDown))
					error_mes_MuonEn = max(abs(error_mes_MuonEnUp), abs(error_mes_MuonEnDown))
					error_mes_TauEn = max(abs(error_mes_TauEnUp), abs(error_mes_TauEnDown))
					error_mes_UnclusteredEn = max(abs(error_mes_UnclusteredEnUp), abs(error_mes_UnclusteredEnDown))
				
				# Final error values
				error_mes_stat = inc_mes_xsec.std_dev
				error_mes_lumi = max(abs(error_mes_lumi_up), abs(error_mes_lumi_down))
				error_mes_syst = sqrt(error_mes_TTJets_scale*error_mes_TTJets_scale
									+ error_mes_TTJets_mass*error_mes_TTJets_mass
									+ error_mes_TTJets_generator*error_mes_TTJets_generator
									+ error_mes_TTJets_hadronisation*error_mes_TTJets_hadronisation
									+ error_mes_TTJets_cross_section*error_mes_TTJets_cross_section
									+ error_mes_TTJets_NLOgenerator*error_mes_TTJets_NLOgenerator
									+ error_mes_SingleTop_cross_section*error_mes_SingleTop_cross_section
									+ error_mes_V_Jets_cross_section*error_mes_V_Jets_cross_section
									+ error_mes_QCD_cross_section*error_mes_QCD_cross_section
									+ error_mes_QCD_shape*error_mes_QCD_shape
									+ error_mes_BJet*error_mes_BJet
									+ error_mes_JER*error_mes_JER
									+ error_mes_JES*error_mes_JES
									+ error_mes_Electron*error_mes_Electron
									+ error_mes_Muon*error_mes_Muon
									)
				if (Variable[var] != 'HT'):
					error_mes_syst = sqrt(error_mes_syst*error_mes_syst
										+ error_mes_ElectronEn*error_mes_ElectronEn
										+ error_mes_MuonEn*error_mes_MuonEn
										+ error_mes_TauEn*error_mes_TauEn
										+ error_mes_UnclusteredEn*error_mes_UnclusteredEn
										)

				print >> outputFile1, "In " + PhaseSpace[ps] + " using the " + Channel[ch] + " channel with the variable " + Variable[var] + ", the measured data inclusive cross section is %.2f " %inc_mes_xsec.nominal_value, "+/- %.2f" %inc_mes_xsec.std_dev, "(stat) +/- %.2f" %error_mes_syst, "(syst) +/- %.2f" %error_mes_lumi, "(lumi)."
			print >> outputFile1, "" 
