from dps.utils.ROOT_utils import set_root_defaults
from argparse import ArgumentParser
from dps.config.xsection import XSectionConfig
from dps.utils.pandas_utilities import file_to_df, matrix_from_df, read_tuple_from_file, dict_to_df, df_to_file, df_to_latexFile
from dps.utils.systematic import print_dictionary
from dps.utils.hist_utilities import value_error_tuplelist_to_hist
from dps.config.latex_labels import variables_latex, measurements_latex
from dps.config.variable_binning import bin_edges_vis, reco_bin_edges_vis
import numpy as np
from ROOT import TMath
import pandas as pd 
from dps.utils.file_utilities import make_folder_if_not_exists
import glob
import os
import gc

import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt
from dps.config import CMS
import matplotlib.gridspec as gridspec

np.set_printoptions(
	precision = 3,
	linewidth = 400,
)




if __name__ == '__main__':
	set_root_defaults()

	measurement_config      = XSectionConfig( 13 )
	path = '{data}/normalisation/background_subtraction/13TeV/{variable}/VisiblePS/{type}/normalisation_muon.txt'

	for variable in measurement_config.variables:
		if 'tau' in variable: continue

		# path_to_muon_qcd_nominal = path.format(data='data', variable=variable, type='central')
		path_to_muon_qcd_nominal = path.format(data='data_MuonIso_0p3_Inf', variable=variable, type='central')
		path_to_muon_qcd_alternate = path.format(data='data_MuonIso_0p3_Inf', variable=variable, type='QCD_other_control_region')
		path_to_muon_qcd_alternate1 = path.format(data='data_MuonIso_0p3_1p0', variable=variable, type='QCD_other_control_region')
		path_to_muon_qcd_alternate2 = path.format(data='data_MuonIso_1p0_Inf', variable=variable, type='QCD_other_control_region')
		path_to_muon_qcd_alternate3 = path.format(data='data_MuonIso_2p0_Inf', variable=variable, type='QCD_other_control_region')
		path_to_muon_qcd_alternate4 = path.format(data='data_MuonIso_3p0_Inf', variable=variable, type='QCD_other_control_region')

		qcd_nominal = read_tuple_from_file( path_to_muon_qcd_nominal )['QCD']
		qcd_alternate = read_tuple_from_file( path_to_muon_qcd_alternate )['QCD']
		qcd_alternate1 = read_tuple_from_file( path_to_muon_qcd_alternate1 )['QCD']
		qcd_alternate2 = read_tuple_from_file( path_to_muon_qcd_alternate2 )['QCD']
		qcd_alternate3 = read_tuple_from_file( path_to_muon_qcd_alternate3 )['QCD']
		qcd_alternate4 = read_tuple_from_file( path_to_muon_qcd_alternate4 )['QCD']

		h_qcd_nominal = value_error_tuplelist_to_hist(qcd_nominal, reco_bin_edges_vis[variable])
		h_qcd_alternate = value_error_tuplelist_to_hist(qcd_alternate, reco_bin_edges_vis[variable])
		h_qcd_alternate1 = value_error_tuplelist_to_hist(qcd_alternate1, reco_bin_edges_vis[variable])
		h_qcd_alternate2 = value_error_tuplelist_to_hist(qcd_alternate2, reco_bin_edges_vis[variable])
		h_qcd_alternate3 = value_error_tuplelist_to_hist(qcd_alternate3, reco_bin_edges_vis[variable])
		h_qcd_alternate4 = value_error_tuplelist_to_hist(qcd_alternate4, reco_bin_edges_vis[variable])


		fig = plt.figure( figsize = ( 20, 16 ), dpi = 400, facecolor = 'white' )

		gs = gridspec.GridSpec( 2, 1, height_ratios = [5, 2] )
		ax = plt.subplot( gs[0] )

		template = '%.1f fb$^{-1}$ (%d TeV)'
		label = template % ( measurement_config.new_luminosity/1000., measurement_config.centre_of_mass_energy)
		plt.title( label,loc='right', **CMS.title )

		ax.minorticks_on()
		ax.xaxis.labelpad = 12
		ax.yaxis.labelpad = 12
		plt.tick_params( **CMS.axis_label_major )
		plt.tick_params( **CMS.axis_label_minor )

		h_qcd_nominal.linewidth = 8
		h_qcd_nominal.markersize = 4
		h_qcd_nominal.color = 'red'
		h_qcd_nominal.linestyle = 'solid'

		h_qcd_alternate.linewidth = 8
		h_qcd_alternate.markersize = 4
		h_qcd_alternate.color = 'green'
		h_qcd_alternate.linestyle = 'solid'

		h_qcd_alternate1.linewidth = 6
		h_qcd_alternate1.markersize = 4
		h_qcd_alternate1.color = 'blue'
		h_qcd_alternate1.linestyle = 'dashed'

		h_qcd_alternate2.linewidth = 6
		h_qcd_alternate2.markersize = 4
		h_qcd_alternate2.color = 'orange'
		h_qcd_alternate2.linestyle = 'dashed'

		h_qcd_alternate3.linewidth = 6
		h_qcd_alternate3.markersize = 4
		h_qcd_alternate3.color = 'hotpink'
		h_qcd_alternate3.linestyle = 'dashed'

		h_qcd_alternate4.linewidth = 6
		h_qcd_alternate4.markersize = 4
		h_qcd_alternate4.color = 'cyan'
		h_qcd_alternate4.linestyle = 'dashed'

		rplt.hist(h_qcd_nominal, stacked=False, axes = ax, label = 'QCD Prediction from Nominal Muon CR (0.15$<$Iso$<$0.3)')
		rplt.hist(h_qcd_alternate, stacked=False, axes = ax, label = 'QCD Prediction from Alternate Muon CR (0.3$<$Iso)')
		rplt.hist(h_qcd_alternate1, stacked=False, axes = ax, label = 'QCD Prediction from Alternate Muon CR (0.3$<$Iso$<$1.0)')
		rplt.hist(h_qcd_alternate2, stacked=False, axes = ax, label = 'QCD Prediction from Alternate Muon CR (1.0$>$Iso)')
		rplt.hist(h_qcd_alternate3, stacked=False, axes = ax, label = 'QCD Prediction from Alternate Muon CR (2.0$>$Iso)')
		rplt.hist(h_qcd_alternate4, stacked=False, axes = ax, label = 'QCD Prediction from Alternate Muon CR (3.0$>$Iso)')
		# rplt.errorbar(qcd_from_data_in_other, axes = ax, label = '')
		# rplt.errorbar(qcd_estimate_from_central, axes = ax, label = '')
		# rplt.errorbar(other_control_histograms['QCD'], axes = ax, label = '')

		ax.set_ylim( ymin = 0.)
		if 'NJets' in variable:
			ax.set_xlim( xmin = 3.5, xmax = 9.5)
		if variable in ['HT', 'ST']:
			ax.set_xlim( xmin = 100., xmax = 1000)
		if variable in ['WPT', 'lepton_pt']:
			ax.set_xlim( xmin = 0., xmax = 500)
		if 'MET' in variable:
			ax.set_xlim( xmin = 0., xmax = 150)
		if 'abs_lepton_eta' in variable:
			ax.set_xlim( xmin = 0., xmax = 2.4)

		leg = plt.legend(
			loc='best',
			prop = {'size':26},	
		)
		leg.draw_frame(False)	
		plt.ylabel( 'Events', CMS.y_axis_title)

		### RATIO
		ratio = h_qcd_alternate / h_qcd_nominal
		ratio.linewidth = 4
		ratio.color = 'green'
		ratio.linestyle = 'solid'

		ratio1 = h_qcd_alternate1 / h_qcd_nominal
		ratio1.linewidth = 4
		ratio1.color = 'blue'
		ratio1.linestyle = 'solid'

		ratio2 = h_qcd_alternate2 / h_qcd_nominal
		ratio2.linewidth = 4
		ratio2.color = 'orange'
		ratio2.linestyle = 'solid'

		ratio3 = h_qcd_alternate3 / h_qcd_nominal
		ratio3.linewidth = 4
		ratio3.color = 'hotpink'
		ratio3.linestyle = 'solid'

		ratio4 = h_qcd_alternate4 / h_qcd_nominal
		ratio4.linewidth = 4
		ratio4.color = 'cyan'
		ratio4.linestyle = 'solid'

		plt.setp( ax.get_xticklabels(), visible = False )
		ax2 = plt.subplot( gs[1] )
		ax2.grid( True, 'major', linewidth = 0.5 )
		ax2.axhline(y=1, linewidth = 2)
		plt.ylabel( r'$\frac{\mathrm{data}}{\mathrm{pred.}}$', CMS.y_axis_title)
		ax2.yaxis.set_label_coords(-0.115, 0.8)
		
		rplt.hist(ratio, axes = ax2)
		rplt.hist(ratio1, axes = ax2)
		rplt.hist(ratio2, axes = ax2)
		rplt.hist(ratio3, axes = ax2)
		rplt.hist(ratio4, axes = ax2)

		ax2.set_ylim( ymin = 0., ymax = 2.)
		if 'NJets' in variable:
			ax2.set_xlim( xmin = 3.5, xmax = 9.5)
		if variable in ['HT', 'ST']:
			ax2.set_xlim( xmin = 100., xmax = 1000)
		if variable in ['WPT', 'lepton_pt']:
			ax2.set_xlim( xmin = 0., xmax = 500)
		if 'MET' in variable:
			ax2.set_xlim( xmin = 0., xmax = 150)
		if 'abs_lepton_eta' in variable:
			ax2.set_xlim( xmin = 0., xmax = 2.4)

		plt.tick_params( **CMS.axis_label_major )
		plt.tick_params( **CMS.axis_label_minor )

		x_title = r''
		x_title += variables_latex[variable]
		if variable in ['HT', 'MET', 'WPT', 'ST', 'lepton_pt']:
			x_title += ' [GeV]'
		plt.xlabel( x_title, CMS.x_axis_title )

		plt.tight_layout()

		# outputFileName = 'plots/QCDvalidation/{var}_{channel}.pdf'.format( var=variable, channel=ch )
		outputFileName = '{var}.pdf'.format( var=variable )
		fig.savefig(outputFileName)