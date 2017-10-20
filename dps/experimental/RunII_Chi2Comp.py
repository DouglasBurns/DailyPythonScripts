from dps.utils.ROOT_utils import set_root_defaults
from argparse import ArgumentParser
from dps.config.xsection import XSectionConfig
from dps.utils.pandas_utilities import file_to_df, matrix_from_df, read_tuple_from_file, dict_to_df, df_to_file, df_to_latexFile
from dps.utils.systematic import print_dictionary
from dps.config.latex_labels import variables_latex, measurements_latex
import numpy as np
from ROOT import TMath
import pandas as pd 
from dps.utils.file_utilities import make_folder_if_not_exists
import glob
import os
import gc

import matplotlib.pyplot as plt
from dps.config import CMS
np.set_printoptions(
	precision = 3,
	linewidth = 400,
)

systematics_latex = {
    'Electron'                : 'Electron efficiency',
    'Muon'                    : 'Muon efficiency',
    'ElectronEn'                : 'Electron energy',
    'MuonEn'                    : 'Muon energy',
    'TauEn'                     : 'Tau energy',
    'UnclusteredEn'             : 'Unclustered energy',
    'PileUp'                    : 'Pile-up',
    'BJet'                      : 'b-tagging efficiency',
    'JES'                       : 'JES',
    'JER'                       : 'JER',
    'luminosity'                : 'Luminosity',
    'V+Jets_cross_section'      : 'V+jets cross section',
    'SingleTop_cross_section'   : 'Single top cross section',
    'QCD_cross_section'         : 'QCD cross section',
    'QCD_shape'                 : 'QCD shape ',
    'PDF'                       : 'PDF ',
    'TTJets_topPt'              : 'Top pT',
    'TTJets_mass'               : 'Top mass',
    'TTJets_scale'              : 'Renormalization and factorization scales',
    'TTJets_matching'           : 'ME/PS matching',
    'TTJets_ue'                 : 'Underlying event tune',
    'TTJets_frag'               : 'Fragmentation',
    'TTJets_petersonFrag'       : 'Alternative fragmentation model',
    'TTJets_hdamp'              : 'hdamp',
    'TTJets_semiLepBr'          : 'B hadron decay semileptonic branching fraction',
    'TTJets_CR'                 : 'Colour reconnection',
    'TTJets_CR_erdOn'           : 'Colour reconnection (erdOn)',
    'TTJets_CR_QCDbased_erdOn'  : 'Colour reconnection (QCD-based erdOn)',
    'TTJets_CR_GluonMove'       : 'Colour reconnection (Gluon move)',
    'TTJets_CR_GluonMove_erdOn' : 'Colour reconnection (Gluon move erdOn)',
    'inputMC'                   : 'MC statistics',
    'central'                   : '',
    'systematic'                : '',
    'statistical'               : '',
    # 'Stat_normalisedXsection_inputMC' : 'MC statistics',
    'Stat' : 'statistics',
    # 'Stat_normalisedXsection' : 'MC statistics',
    'Total' : 'Total',
}

class chi2Info:
	def __init__(self, chi2, ndf, pValue):
		self.chi2 = chi2
		self.ndf = ndf
		self.pValue = pValue

def calculateChi2(measured_xsection,model_xsection,covariance):
	diff = (measured_xsection - model_xsection)
	print np.linalg.cond(covariance)

	# SVD decomposition doesnt help
	U, s, Vt = np.linalg.svd(covariance, full_matrices=False) 
	inv = np.transpose(Vt).dot( np.linalg.inv(np.diag(s)) ).dot(np.transpose(U))
	chi2 = np.transpose( diff ).dot( inv ).dot( diff )
	print covariance.dot( inv )

	# QR decomposition doesnt help
	# q, r = np.linalg.qr(covariance) 
	# qt = np.transpose(q).dot(diff)
	# chi2b = np.transpose( diff ).dot(np.linalg.solve( r, qt ) )

	# Numerical Inverse doesnt help
	# chi2 = np.transpose( diff ).dot( np.linalg.inv( covariance ) ).dot( diff )

	# Pseudo Inverse doesnt help
	# chi2 = np.transpose( diff ).dot( np.linalg.pinv( covariance ) ).dot( diff )

	ndf = len(measured_xsection)
	prob = TMath.Prob( chi2, len(measured_xsection) )
	return chi2Info( chi2, ndf, prob ), covariance.dot( inv )

def print_CovCovinv(m, t, utype, var, channel):

	fig = plt.figure( figsize = (16,16), dpi = CMS.dpi, facecolor = CMS.facecolor )
	ax = fig.add_subplot(1, 1, 1)
	ax.set_aspect('equal')
	im=plt.imshow(m, interpolation='nearest', cmap=plt.cm.ocean)

	plt.tick_params( **CMS.axis_label_major )
	plt.tick_params( **CMS.axis_label_minor )
	cbar=plt.colorbar(im,fraction=0.046, pad=0.04)
	cbar.ax.tick_params(labelsize=30) 
	title = '{} {} xsec in {} channel \n {}'.format(utype, variables_latex[var], channel, systematics_latex[t])
	plt.title( title , loc='right', **CMS.title )
	plt.tight_layout()
	output_folder = 'plots/chi2/'
	make_folder_if_not_exists(output_folder)
	file_template = output_folder + 'CovCovInv_{}_{}_{}_{}.pdf'.format(
		var,
		utype,
		channel,
		t,
	)
	fig.savefig(file_template)
	fig.clf()
	plt.close()
	gc.collect()
	return

def print_plot(utype, var, d_Plot, label=None):
	# tmp Labels

	names = [ systematics_latex[x] for x in d_Plot.keys() ]
	values = [ x for x in d_Plot.values() ]
	values_neg = [ -x for x in d_Plot.values() ]

	ints = range(1, len(names)+1)

	fig = plt.figure( figsize = (16,24), dpi = CMS.dpi, facecolor = CMS.facecolor )
	ax = fig.add_subplot(1, 1, 1)
	plt.bar(ints, values, align='center', color='blue')
	# plt.bar(ints, values_neg, align='center', color='red')
	plt.xticks(ints, names, rotation=90)
	plt.tick_params( **CMS.axis_label_major )
	plt.tick_params( **CMS.axis_label_minor )

	ax.set_yscale("symlog")
	# ax.set_yscale("log", nonposy='clip')
	# clip sets all -ve values to a very small positive one

	plt.xlabel( 'Systematic Variation', CMS.x_axis_title )
	plt.ylabel( 'Chi2', CMS.y_axis_title )
	plt.title( label.replace('_', ' '), loc='right', **CMS.title )

	plt.tight_layout()
	output_folder = 'plots/chi2/'
	make_folder_if_not_exists(output_folder)
	file_template = output_folder + '{var}_{type}'.format(
		var = var, 
		type = utype,
	)
	if label:
		file_template += '_{}.pdf'.format(label)
	else:
		file_template += '.pdf'
	fig.savefig(file_template)
	return

def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument( "-p", "--path", 
        dest    = "path", 
        default = 'data/normalisation/background_subtraction',
        help    = "set path to files containing dataframes" 
    )
    parser.add_argument( '--variable','-v', 
        dest    = "varToDO",
        default = None,
        help    = "DEBUG: Test a particular variable" 
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
	set_root_defaults()
	args = parse_arguments()

	measurement_config      = XSectionConfig( 13 )
	
	unc_type = [
		'normalised',
		'absolute',
	]
	# args.varToDO = 'HT'
	for utype in unc_type:
		for variable in measurement_config.variables:
			if args.varToDO and variable != args.varToDO: continue
			d_Chi2 = {}

			path_to_measurement = '{path}/13TeV/{variable}/VisiblePS/central'.format(
			    path = args.path, 
			    variable = variable,
			)
			electron_measurement = '{path}/xsection_{type}_electron_TUnfold.txt'.format(
				path=path_to_measurement, 
				type = utype, 
			)
			muon_measurement = '{path}/xsection_{type}_muon_TUnfold.txt'.format(
				path=path_to_measurement, 
				type = utype, 
			)
			combined_measurement = '{path}/xsection_{type}_combined_TUnfold.txt'.format(
				path=path_to_measurement, 
				type = utype, 
			)

			path_to_covariance = '{path}/covarianceMatrices/{type}'.format(
				path=path_to_measurement, 
				type = utype, 
			)
			electron_covariances = glob.glob(
				'{path}/*_Covariance_electron.txt'.format(
					path=path_to_covariance, 
				)
			)
			muon_covariances = glob.glob(
				'{path}/*_Covariance_muon.txt'.format(
					path=path_to_covariance, 
				)
			)
			combined_covariances = glob.glob(
				'{path}/*_Covariance_combined.txt'.format(
					path=path_to_covariance, 
				)
			)

			# Read Cross Sections
			e_XSec 		= read_tuple_from_file( electron_measurement )
			mu_XSec 	= read_tuple_from_file( muon_measurement )
			c_XSec 		= read_tuple_from_file( combined_measurement )

			# Data
			# Strip uncertainties
			e_XSec_Unf 	= np.array([ i[0] for i in e_XSec['TTJets_unfolded'] ])
			mu_XSec_Unf = np.array([ i[0] for i in mu_XSec['TTJets_unfolded'] ])
			c_XSec_Unf 	= np.array([ i[0] for i in c_XSec['TTJets_unfolded'] ])

			# Model
			# Strip uncertainties
			e_XSec_Mod 	= np.array([ i[0] for i in e_XSec['TTJets_powhegPythia8'] ])
			mu_XSec_Mod = np.array([ i[0] for i in mu_XSec['TTJets_powhegPythia8'] ])
			c_XSec_Mod 	= np.array([ i[0] for i in c_XSec['TTJets_powhegPythia8'] ])

			d_Chi2['electron_xsection_unfolded'] 	= e_XSec_Unf
			d_Chi2['electron_xsection_model'] 		= e_XSec_Mod
			d_Chi2['muon_xsection_unfolded'] 		= mu_XSec_Unf
			d_Chi2['muon_xsection_model'] 			= mu_XSec_Mod
			d_Chi2['combined_xsection_unfolded'] 	= c_XSec_Unf
			d_Chi2['combined_xsection_model'] 		= c_XSec_Mod

			d_Plot_diff = {}
			d_Plot_eChi2 = {}
			d_Plot_muChi2 = {}
			d_Plot_cChi2 = {}
			for e_Cov, mu_Cov, c_Cov in zip(electron_covariances, muon_covariances, combined_covariances):
				syst_Name_e = e_Cov.split('/')[-1].replace('_Covariance_electron.txt', '')
				syst_Name_mu = mu_Cov.split('/')[-1].replace('_Covariance_muon.txt', '')
				syst_Name_c = c_Cov.split('/')[-1].replace('_Covariance_combined.txt', '')

				if syst_Name_e != syst_Name_mu != syst_Name_c:
					print syst_Name_e
					print syst_Name_mu
					print syst_Name_c
					print "Electron and muon channel covariances are not syncronised"
					os.exit()
				else:
					syst_Name = syst_Name_e
					d_Chi2[syst_Name] = {}

				if variable in measurement_config.variables_no_met and syst_Name in measurement_config.systematic_group_met:
					continue
				if 'Stat_' in syst_Name: continue

				# Read Matrices
				e_Cov_Mat 		= matrix_from_df( file_to_df(e_Cov) )
				mu_Cov_Mat 		= matrix_from_df( file_to_df(mu_Cov) )
				c_Cov_Mat 		= matrix_from_df( file_to_df(c_Cov) )
				print '\n', syst_Name

				if syst_Name == 'Electron':
					e_Chi2_Info, e 	= calculateChi2(e_XSec_Unf, e_XSec_Mod, e_Cov_Mat)
					mu_Chi2_Info, mu 	= chi2Info( 0, 0, 0 ), [[0]]
					c_Chi2_Info, c 	= calculateChi2(c_XSec_Unf, c_XSec_Mod, c_Cov_Mat)
				elif syst_Name == 'Muon':
					e_Chi2_Info, e 	= chi2Info( 0, 0, 0 ), [[0]]
					mu_Chi2_Info, mu  	= calculateChi2(mu_XSec_Unf, mu_XSec_Mod, mu_Cov_Mat)
					c_Chi2_Info, c 	= calculateChi2(c_XSec_Unf, c_XSec_Mod, c_Cov_Mat)
				else:
					e_Chi2_Info, e 	= calculateChi2(e_XSec_Unf, e_XSec_Mod, e_Cov_Mat)
					mu_Chi2_Info, mu  	= calculateChi2(mu_XSec_Unf, mu_XSec_Mod, mu_Cov_Mat)
					c_Chi2_Info, c 	= calculateChi2(c_XSec_Unf, c_XSec_Mod, c_Cov_Mat)
				print_CovCovinv(e, syst_Name, utype, variable, 'electron')
				print_CovCovinv(mu, syst_Name, utype, variable, 'muon')
				print_CovCovinv(c, syst_Name, utype, variable, 'combined')


				d_Chi2[syst_Name]['electron_Covariance'] 	= e_Cov_Mat
				d_Chi2[syst_Name]['muon_Covariance'] 		= mu_Cov_Mat
				d_Chi2[syst_Name]['combined_Covariance'] 	= mu_Cov_Mat
				d_Chi2[syst_Name]['electron_Chi2'] 			= e_Chi2_Info.chi2
				d_Chi2[syst_Name]['muon_Chi2'] 				= mu_Chi2_Info.chi2
				d_Chi2[syst_Name]['combined_Chi2'] 			= c_Chi2_Info.chi2
				d_Chi2[syst_Name]['electron_ndf'] 			= e_Chi2_Info.ndf
				d_Chi2[syst_Name]['muon_ndf'] 				= mu_Chi2_Info.ndf
				d_Chi2[syst_Name]['combined_ndf'] 			= c_Chi2_Info.ndf
				d_Chi2[syst_Name]['electron_pValue'] 		= e_Chi2_Info.pValue
				d_Chi2[syst_Name]['muon_pValue'] 			= mu_Chi2_Info.pValue
				d_Chi2[syst_Name]['combined_pValue'] 		= c_Chi2_Info.pValue
				d_Chi2[syst_Name]['diff_chi2'] 				= abs(e_Chi2_Info.chi2 - mu_Chi2_Info.chi2)
				d_Chi2[syst_Name]['diff_pValue'] 			= abs(e_Chi2_Info.pValue - mu_Chi2_Info.pValue)
				d_Plot_diff[syst_Name] 						= d_Chi2[syst_Name]['diff_chi2']
				d_Plot_eChi2[syst_Name] 					= d_Chi2[syst_Name]['electron_Chi2']
				d_Plot_muChi2[syst_Name] 					= d_Chi2[syst_Name]['muon_Chi2']
				d_Plot_cChi2[syst_Name]						= d_Chi2[syst_Name]['combined_Chi2']

			# print_dictionary( 'Hi', d_Chi2 )
			print_plot(utype, variable, d_Plot_diff, label = 'Difference')
			print_plot(utype, variable, d_Plot_eChi2, label = 'Electron_Chi2')
			print_plot(utype, variable, d_Plot_muChi2, label = 'Muon_Chi2')
			print_plot(utype, variable, d_Plot_cChi2, label = 'Combined_Chi2')



