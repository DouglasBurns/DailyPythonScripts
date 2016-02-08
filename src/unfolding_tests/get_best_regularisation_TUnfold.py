'''
Created on 9 Mar 2015

@author: kreczko

This script creates
 - global correlation plots for tau (float) & k (int) regularisation parameters
 - optimal tau & k values
 - d_i plots with fits for k values
 
What it needs:
 - a set of four histograms:
   - truth: quantity at generator level
   - gen_vs_reco: quantity after selection, generated vs reco
   - measured: measured (reco) quantity including fakes (background)
   - data
For example config files, please have a look at config/unfolding/*.json
   
usage:
    python get_best_regularisation.py config.json
    # for full 7 + 8 TeV stuff:
    python src/unfolding_tests/get_best_regularisation.py config/unfolding/*.json -c
'''
# imports
from __future__ import division
from math import log10, pow
from optparse import OptionParser
import sys
# rootpy
from rootpy.io import File
from rootpy.plotting import Graph
from rootpy.matrix import Matrix
# DailyPythonScripts
from tools.file_utilities import read_data_from_JSON, make_folder_if_not_exists
from tools.Unfolding import Unfolding, get_unfold_histogram_tuple
#from src.cross_section_measurement.lib import get_unfold_histogram_tuple
from tools.ROOT_utils import set_root_defaults, get_histogram_from_file
from config import XSectionConfig
from config.variable_binning import bin_edges, bin_edges_vis
from tools.hist_utilities import value_error_tuplelist_to_hist
from tools.table import PrintTable
import matplotlib.pyplot as plt
from tools.plotting import Histogram_properties
from matplotlib import rc
from config import CMS
from config.latex_labels import variables_latex
from ROOT import TGraph, TSpline3, Double, TUnfoldDensity, TUnfold, TDecompSVD, TMatrixD
from numpy.linalg import svd
rc('font',**CMS.font)
rc( 'text', usetex = True )

class RegularisationSettings():
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
        if input_values.has_key('n_tau_scan_points'):
            self.n_tau_scan_points = input_values['n_tau_scan_points']
        if input_values.has_key('n_toy'):
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

        t, m, r, _ = get_unfold_histogram_tuple( File(input_file),
                                              self.variable,
                                              self.channel,
                                              centre_of_mass = self.centre_of_mass_energy,
                                              ttbar_xsection=self.measurement_config.ttbar_xsection,
                                              luminosity=self.measurement_config.luminosity,
                                              load_fakes = False,
                                              visiblePS = visiblePS
                                            )
        self.h_truth = t
        self.h_response = r
        self.h_measured = m
        
        data_file = self.data['file']
        if data_file.endswith('.root'):
            self.h_data = get_histogram_from_file(self.data['histogram'], self.data['file'])
        elif data_file.endswith('.json') or data_file.endswith('.txt'):
            data_key = self.data['histogram']
            # assume configured bin edges
            edges = []
            if self.phaseSpace == 'FullPS':
                edges = bin_edges[self.variable]
            elif self.phaseSpace == 'VisiblePS':
                edges = bin_edges_vis[self.variable]
            json_input = read_data_from_JSON(data_file)

            if data_key == "": # JSON file == histogram
                self.h_data = value_error_tuplelist_to_hist(json_input, edges)
            else:
                self.h_data = value_error_tuplelist_to_hist(json_input[data_key], edges)
        else:
            print 'Unkown file extension', data_file.split('.')[-1]
            
    def get_histograms( self ):
        return self.h_truth, self.h_response, self.h_measured, self.h_data
    
def main():
    options, input_values_sets, json_input_files = parse_options()
    use_current_k_values = options.compare
    results = {}
    for input_values, json_file in zip( input_values_sets, json_input_files ):
        print 'Processing', json_file
        regularisation_settings = RegularisationSettings( input_values )
        variable = regularisation_settings.variable
        channel = regularisation_settings.channel
        com = regularisation_settings.centre_of_mass_energy
        if not results.has_key(com): results[com] = {}
        if not results[com].has_key(variable): results[com][variable] = {}
        print 'Variable = %s, channel = "%s", sqrt(s) = %d' % (variable, channel, com)

        h_truth, h_response, h_measured, h_data = regularisation_settings.get_histograms()

        unfolding = Unfolding( 
                                h_data, 
                                h_truth, 
                                h_measured, 
                                h_response,
                                fakes = None,
                                method = 'TUnfold', 
                                k_value = -1, 
                                tau = 0. 
                            )
        # h_response.Draw('COLZ')
        # raw_input('...')
        get_condition_number( unfolding.unfoldObject )

        # tau_results = get_best_tau( regularisation_settings )
        # results[com][variable][channel] = (tau_results)
        # plot(regularisation_settings, (k_results, tau_results), 
        #          use_current_k_values)
    print 'Results :',results
    # table(results, use_current_k_values, options.style)

def parse_options():
    parser = OptionParser( __doc__ )
    parser.add_option( "-c", "--compare", dest = "compare", action = "store_true",
                      help = "Compare to current values (k vs tau)", default = False )
    parser.add_option( "-t", "--table-style", dest = "style", default = 'simple',
                      help = "Style for table printing: simple|latex|twiki (default = simple)" )

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

def removeFakes( h_measured, h_data, h_response ):
    fakes = h_measured - h_response.ProjectionX()
    nonFakeRatio = 1 - fakes / h_measured
    h_measured *= nonFakeRatio
    h_data *= nonFakeRatio
    return h_measured, h_data

def tau_from_L_curve( unfoldingObject ):
    lCurve = TGraph()
    logTauX = TSpline3()
    logTauY = TSpline3()
    print 'Scanning'
    iBest = unfoldingObject.ScanLcurve(100,0.00001,100.0, lCurve, logTauX, logTauY);
    print 'Done'
    # Additional info, plots
    # t = Double(0)
    # x = Double(0)
    # y = Double(0)
    # logTauX.GetKnot(iBest,t,x);
    # logTauY.GetKnot(iBest,t,y);

    # bestLcurve = Graph(1);
    # bestLcurve.SetPoint(1,x,y);
    # bestLogTauLogChi2 = Graph(1);
    # bestLogTauLogChi2.SetPoint(1,t,x);

    # logTauX.Draw()
    # bestLogTauLogChi2.markercolor = 'red'
    # # bestLogTauLogChi2.Draw("*")
    # lCurve.Draw("AL");
    # raw_input('...')
    # bestLcurve.markercolor = 'red'
    # bestLcurve.Draw("*");

    return unfoldingObject.GetTau()

def tau_from_scan( unfoldingObject ):
    lCurve = TGraph()
    scanResult = TSpline3()
    # logTauY = TSpline3()
    d = 'signal'
    a = ''
    nScan = 20
    iBest = unfoldingObject.ScanTau(nScan,0.0001,10., scanResult, TUnfoldDensity.kEScanTauRhoAvg);
    # scanResult.Draw()
    # raw_input('...')

    # t = Double(0)
    # x = Double(0)
    # scanResult.GetKnot(iBest,t,x);
    # bestTau = TGraph(1)
    # print 'Best :',t,x
    # bestTau.SetPoint(1,t,x)

    # knots = TGraph(nScan)
    # for i in range( 0, nScan ) :
    #     t = Double(0)
    #     x = Double(0)
    #     y = Double(0)
    #     scanResult.GetKnot(i,t,x);
    #     knots.SetPoint(i,t,x)
    # knots.Draw('*')
    # bestTau.markercolor = 'red'
    # bestTau.Draw('*')
    # raw_input('...')
    # Additional info, plots
    # t = Double(0)
    # x = Double(0)
    # y = Double(0)
    # logTauX.GetKnot(iBest,t,x);
    # logTauY.GetKnot(iBest,t,y);

    # bestLcurve = Graph(1);
    # bestLcurve.SetPoint(1,x,y);
    # bestLogTauLogChi2 = Graph(1);
    # bestLogTauLogChi2.SetPoint(1,t,x);

    # logTauX.Draw()
    # bestLogTauLogChi2.markercolor = 'red'
    # # bestLogTauLogChi2.Draw("*")
    # lCurve.Draw("AL");
    # bestLcurve.markercolor = 'red'
    # bestLcurve.Draw("*");

    return unfoldingObject.GetTau()

def get_best_tau( regularisation_settings ):
    '''
        returns TODO
         - optimal_tau: TODO
    '''
    h_truth, h_response, h_measured, h_data = regularisation_settings.get_histograms()

    h_measured, h_data = removeFakes( h_measured, h_data, h_response)

    unfolding = Unfolding( 
                            h_data, 
                            h_truth, 
                            h_measured, 
                            h_response,
                            fakes = None,
                            method = 'TUnfold', 
                            k_value = -1, 
                            tau = 1.12332403298 
                        )

    bestTau_LCurve = tau_from_L_curve( unfolding.unfoldObject )
    unfolding.tau = bestTau_LCurve

    return unfolding.tau

def get_condition_number( unfoldingObject ):

    probMatrix = unfoldingObject.GetProbabilityMatrix( 'ProbMatrix')
    # probabilityMatrix = unfoldingObject.GetProbabilityMatrix( probMatrix, TUnfold.kHistMapOutputVert)
    probMatrix.Print()
    probMatrix.Draw('COLZ')
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
    condition = sigmaMax / max(0,sigmaMin)
    # condition = 1
    print condition
    return condition

if __name__ == '__main__':
    set_root_defaults( set_batch = False, msg_ignore_level = 3001 )
    main()
