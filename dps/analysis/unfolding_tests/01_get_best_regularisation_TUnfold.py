'''
Created on 9 Mar 2015

@author: kreczko

This script creates, for each variable:
 - global correlation as function of log tau
 - optimal tau value

 
What it needs:
 - the response matrix file
 - the data to be unfolded
   
usage:
    python get_best_regularisation.py config.json
    # for 13 TeV in the visible phase space :
    python dps/analysis/unfolding_tests/01_get_best_regularisation_TUnfold.py config/unfolding/VisiblePS/*.json
'''
# imports
from __future__ import division
from optparse import OptionParser
import sys
# rootpy
from rootpy.io import File
from rootpy.plotting import Graph
# DailyPythonScripts
from dps.utils.file_utilities import read_data_from_JSON, make_folder_if_not_exists
from dps.utils.Unfolding import Unfolding, get_unfold_histogram_tuple, removeFakes
#from dps.analysis.xsection.lib import get_unfold_histogram_tuple
from dps.utils.ROOT_utils import set_root_defaults, get_histogram_from_file
from dps.config.xsection import XSectionConfig
from dps.config.variable_binning import reco_bin_edges_vis
from dps.utils.hist_utilities import value_error_tuplelist_to_hist
from matplotlib import rc
from dps.config import CMS
from ROOT import TGraph, TSpline3, Double, TUnfoldDensity, TUnfold, TDecompSVD, TMatrixD, TCanvas, gROOT
from rootpy import asrootpy
from dps.utils.pandas_utilities import read_tuple_from_file
from math import log

rc('font',**CMS.font)
rc( 'text', usetex = True )

class RegularisationSettings(object):
    '''
        Class for storing input configuration, and for getting and storing response matrix histograms
    '''
    n_toy = int( 1000 )
    n_tau_scan_points = int( 100 )
    
    def __init__( self, input_values ):
        self.centre_of_mass_energy = input_values['centre-of-mass energy']
        self.measurement_config = XSectionConfig( self.centre_of_mass_energy )
        self.channel = input_values['channel']
        self.variable = input_values['variable']
        self.phaseSpace = input_values['phaseSpace']
        self.output_folder = input_values['output_folder']
        self.output_format = input_values['output_format']
        self.truth = input_values['truth']
        self.gen_vs_reco = input_values['gen_vs_reco']
        self.measured = input_values['measured']
        self.data = input_values['data']

        # optional
        if 'n_tau_scan_points' in input_values:
            self.n_tau_scan_points = input_values['n_tau_scan_points']
        if 'n_toy' in input_values:
            self.n_toy = input_values['n_toy']
            
        self.__set_unfolding_histograms__()
        
    def __set_unfolding_histograms__( self ):
        # at the moment only one file is supported for the unfolding input
        files = set( [self.truth['file'],
                     self.gen_vs_reco['file'],
                     self.measured['file']]
                    )
        if len( files ) > 1:
            print "Currently not supported to have different files for truth, gen_vs_reco and measured"
            sys.exit()
            
        input_file = files.pop()

        visiblePS = False
        if self.phaseSpace == 'VisiblePS':
            visiblePS = True

        t, m, r, f = get_unfold_histogram_tuple( File(input_file),
                                              self.variable,
                                              self.channel,
                                              centre_of_mass = self.centre_of_mass_energy,
                                              ttbar_xsection=self.measurement_config.ttbar_xsection,
                                              luminosity=self.measurement_config.luminosity,
                                              load_fakes = True,
                                              visiblePS = visiblePS
                                            )
        self.h_truth = asrootpy ( t )
        self.h_response = asrootpy ( r )
        self.h_measured = asrootpy ( m )
        self.h_fakes = asrootpy ( f )

        data_file = self.data['file']

        edges = reco_bin_edges_vis[self.variable]
        self.h_data = asrootpy(value_error_tuplelist_to_hist( read_tuple_from_file(data_file)['TTJet'], edges ))
        self.h_data_no_fakes = asrootpy( removeFakes( self.h_measured, self.h_fakes, self.h_data ) )

        scale = self.h_data_no_fakes.integral(0,-1) / asrootpy(self.h_response.ProjectionY()).integral(0,-1)
        self.h_response.Scale(scale)

    def get_histograms( self ):
        return self.h_truth, self.h_response, self.h_measured, self.h_data, self.h_fakes
    
def main():
    options, input_values_sets, json_input_files = parse_options()
    results = {}
    for input_values, json_file in zip( input_values_sets, json_input_files ):
        # print 'Processing', json_file
        if 'combined' in json_file: continue
        regularisation_settings = RegularisationSettings( input_values )
        variable = regularisation_settings.variable
        channel = regularisation_settings.channel
        com = regularisation_settings.centre_of_mass_energy
        if not results.has_key(com): results[com] = {}
        if not results[com].has_key(channel): results[com][channel] = {}
        if not results[com][channel].has_key(variable): results[com][channel][variable] = {}
        print 'Variable = {0}, channel = {1}'.format(variable, channel)

        h_truth, h_response, h_measured, h_data, h_fakes = regularisation_settings.get_histograms()

        unfolding = Unfolding( 
            h_data, 
            h_truth, 
            h_measured, 
            h_response,
            fakes = None,
            method = 'TUnfold', 
            tau = 0. 
        )
        # print unfolding.getConditionNumber()

        tau_results = get_best_tau( regularisation_settings )
        results[com][channel][variable] = (tau_results)
    print_results_to_screen(results)

def parse_options():
    parser = OptionParser( __doc__ )

    ( options, args ) = parser.parse_args()
    
    input_values_sets = []
    json_input_files = []
    add_set = input_values_sets.append
    add_json_file = json_input_files.append
    for arg in args:
        input_values = read_data_from_JSON( arg )
        add_set( input_values )
        add_json_file( arg )

    return options, input_values_sets, json_input_files

def tau_from_L_curve( unfoldingObject, regularisation_settings ):
    '''
    Get best tau via l curve method
    Not tested
    '''
    variable = regularisation_settings.variable
    if variable == 'NJets': return 0

    nScan = 500
    minTau = 1.E-6
    maxTau = 1.E-0

    if variable == 'abs_lepton_eta':
        minTau = 1.E-8
        maxTau = 1.E-3
    elif variable == 'lepton_pt':
        minTau = 1.E-6
        maxTau = 1.E-2
    elif variable == 'NJets':
        minTau = 1.E-6
        maxTau = 1.E-2
    minTau = 0.
    maxTau = 0.

    lCurve = TGraph()
    logTauX = TSpline3()
    logTauY = TSpline3()
    logTauCurvature = TSpline3() # it should be a peaked function (similar to a Gaussian), the maximum corresponding to the final choice of tau
    iBest = unfoldingObject.ScanLcurve(nScan, minTau, maxTau, lCurve, logTauX, logTauY, logTauCurvature)

    # Additional info, plots
    t = Double(0)
    x = Double(0)
    y = Double(0)
    logTauX.GetKnot(iBest,t,x)
    logTauY.GetKnot(iBest,t,y)
    
    canvas = TCanvas()

    bestLcurve = Graph(1)
    bestLcurve.SetPoint(1,x,y)
    bestLcurve.markercolor = 'red'
    bestLcurve.SetMarkerSize(1.5)
    bestLcurve.SetMarkerStyle(34)
    bestLcurve.GetXaxis().SetTitle('log(#tau)')
    bestLcurve.GetYaxis().SetTitle('Something to do with the LCurve')
    bestLcurve.SetTitle('{0} {1}'.format(variable, regularisation_settings.channel))
    bestLcurve.Draw('AP')

    lCurve.SetMarkerColor(600)
    lCurve.SetMarkerSize(0.5)
    lCurve.SetMarkerStyle(20)
    lCurve.Draw('LPSAME')

    bestLcurve.Draw('PSAME')

    # Write to file
    output_dir = regularisation_settings.output_folder
    make_folder_if_not_exists(output_dir)
    canvas.SaveAs(output_dir + '/{0}_{1}_LCurve.pdf'.format(variable, regularisation_settings.channel) )

    canvas2 = TCanvas()
    bestLcurve = Graph(1)
    bestLcurve.SetPoint(1,x,y)
    bestLcurve.markercolor = 'red'
    bestLcurve.SetMarkerSize(1.5)
    bestLcurve.SetMarkerStyle(34)
    bestLcurve.GetXaxis().SetTitle('log(#tau)')
    bestLcurve.GetYaxis().SetTitle('Something to do with the LCurve')
    bestLcurve.SetTitle('{0} {1}'.format(variable, regularisation_settings.channel))
    bestLcurve.Draw('AP')

    logTauCurvature.SetMarkerColor(600)
    logTauCurvature.SetMarkerSize(0.5)
    logTauCurvature.SetMarkerStyle(20)
    logTauCurvature.Draw("LPSAME")

    bestLcurve.Draw('PSAME')

    # Write to file
    output_dir = regularisation_settings.output_folder
    make_folder_if_not_exists(output_dir)
    canvas2.SaveAs(output_dir + '/{0}_{1}_LCurveCurvature.pdf'.format(variable, regularisation_settings.channel) )

    return unfoldingObject.GetTau()

def tau_from_scan( unfoldingObject, regularisation_settings ):
    variable = regularisation_settings.variable

    # Plots that get outputted by the scan
    scanResult = TSpline3()
    d = 'signal'
    a = ''

    # Parameters of scan
    # Number of points to scan, and min/max tau
    nScan = 200
    minTau = 1.E-6
    maxTau = 1.E-0

    if variable == 'abs_lepton_eta':
        minTau = 1.E-8
        maxTau = 1.E-3
    elif variable == 'lepton_pt':
        minTau = 1.E-6
        maxTau = 1.E-2
    elif variable == 'NJets':
        minTau = 1.E-6
        maxTau = 1.E-2

    # Scan is performed here    
    iBest = unfoldingObject.ScanTau(nScan, minTau, maxTau, scanResult, TUnfoldDensity.kEScanTauRhoSquareAvg);
    t = Double(0)
    x = Double(0)
    scanResult.GetKnot(iBest,t,x);

    # Plot the scan result
    # Correlation as function of log tau
    canvas = TCanvas()

    bestTau = Graph(1)
    bestTau.SetPoint(1,t,x)
    bestTau.markercolor = 'red'
    bestTau.SetMarkerSize(1.5)
    bestTau.SetMarkerStyle(34)
    bestTau.GetXaxis().SetTitle('log(#tau)')
    bestTau.GetYaxis().SetTitle('Average global correlation coefficient squared')
    bestTau.SetTitle('{0} {1}'.format(variable, regularisation_settings.channel))
    bestTau.GetYaxis().SetRangeUser(x*0.8,0.95)
    bestTau.GetXaxis().SetLimits(log(minTau, 10), log(maxTau, 10))
    bestTau.Draw('AP')

    scanResult.SetMarkerColor(600)
    scanResult.SetMarkerSize(0.5)
    scanResult.SetMarkerStyle(20)
    scanResult.Draw('LPSAME')
    # Redraw to get it to appear on top of TSpline3...
    bestTau.Draw('PSAME')

    # Write to file
    output_dir = regularisation_settings.output_folder
    make_folder_if_not_exists(output_dir)
    canvas.SaveAs(output_dir + '/{0}_{1}.pdf'.format(variable, regularisation_settings.channel) )

    return unfoldingObject.GetTau()

def get_best_tau( regularisation_settings ):
    '''
        returns TODO
         - optimal_tau: TODO
    '''
    h_truth, h_response, h_measured, h_data, h_fakes = regularisation_settings.get_histograms()
    variable = regularisation_settings.variable

    unfolding = Unfolding( 
                            h_data, 
                            h_truth, 
                            h_measured, 
                            h_response,
                            fakes = None,
                            method = 'TUnfold', 
                            tau = -1
                        )

    # bestTau_LCurve = tau_from_L_curve( unfolding.unfoldObject, regularisation_settings  )
    # unfolding.tau = bestTau_LCurve
    bestTauScan = tau_from_scan( unfolding.unfoldObject, regularisation_settings )
    unfolding.tau = bestTauScan
    # print "Tau L Curve = {}, Min Ave Global Corr Coeff = {}".format(bestTau_LCurve, bestTauScan)

    return unfolding.tau

def get_condition_number( unfoldingObject ):

    probMatrix = unfoldingObject.GetProbabilityMatrix( 'ProbMatrix')
    # probabilityMatrix = unfoldingObject.GetProbabilityMatrix( probMatrix, TUnfold.kHistMapOutputVert)
    probMatrix.Print()
    # probMatrix.Draw('COLZ')
    raw_input()
    nBins = probMatrix.GetNbinsX()
    m = TMatrixD( nBins, nBins )
    for xbin in range(1, probMatrix.GetNbinsX() ):
        for ybin in range(1,probMatrix.GetNbinsY() ):
            m[xbin, ybin] = probMatrix.GetBinContent( xbin, ybin)
    svd = TDecompSVD( m )
    svd.Decompose()
    svd.Print()
    sig = svd.GetSig()
    sig.Print()
    nSig = len(sig)
    sigmaMax = sig[0]
    sigmaMin = sig[nSig-2]
    condition = sigmaMax / max(0.00001,sigmaMin)
    # condition = 1
    print condition
    return condition

def print_results_to_screen(result_dict):
    '''
        Print the results to the screen
        Can copy straight into config
    '''
    print "\n Tau Scan Outcomes: \n"
    for com in result_dict.keys():
        for channel in result_dict[com].keys():
            # Print in foprm such that neatly copy and paste into xsection.py
            print "\t\tself.tau_values_{ch} = {{".format(ch = channel)
            for variable in result_dict[com][channel].keys():
                print '\t\t\t"{0}" : {1},'.format(variable, result_dict[com][channel][variable])
            print "\t\t}"

if __name__ == '__main__':
    set_root_defaults( set_batch = True, msg_ignore_level = 3001 )
    main()
