'''
	This module unfolds the data with different response matrices
'''
from rootpy.io import File

from dps.config.variable_binning import reco_bin_edges_vis, bin_edges_vis, bin_widths_visiblePS
from dps.utils.Unfolding import Unfolding, get_unfold_histogram_tuple, removeFakes
from dps.utils.hist_utilities import hist_to_value_error_tuplelist, value_error_tuplelist_to_hist, values_and_errors_to_hist
from dps.utils.Calculation import calculate_normalised_xsection
from dps.config.xsection import XSectionConfig
from dps.utils.plotting import compare_measurements, Histogram_properties
from dps.config import latex_labels
from rootpy import asrootpy
from collections import OrderedDict
from dps.utils.latex import setup_matplotlib
from dps.utils.pandas_utilities import file_to_df, read_tuple_from_file


# latex, font, etc
setup_matplotlib()
def main():
	config = XSectionConfig(13)

	file_for_powhegPythia  		= File(config.unfolding_powheg_pythia8, 'read')
	file_for_amcatnlo_pythia   	= File(config.unfolding_amcatnlo_pythia8, 'read')
	file_for_powhegHerwig       = File(config.unfolding_powheg_herwig, 'read')
	file_for_madgraph 		    = File(config.unfolding_madgraphMLM, 'read')
	file_for_data       		= '/storage/ec6821/DailyPythonScripts/new/DailyPythonScripts/data/normalisation/background_subtraction/13TeV/{var}/VisiblePS/central/normalisation_{ch}.txt'

	samples_and_files_to_compare = {}
	samples_and_files_to_compare['powhegPythia'] 	= { 'file' : file_for_powhegPythia }
	samples_and_files_to_compare['amcatnloPythia'] 	= { 'file' : file_for_amcatnlo_pythia }
	samples_and_files_to_compare['powhegHerwig'] 	= { 'file' : file_for_powhegHerwig }
	samples_and_files_to_compare['madgraph'] 		= { 'file' : file_for_madgraph }
	samples_and_files_to_compare['data'] 			= { 'file' : file_for_data }
	
	channels = [
		'electron',
		'muon',
	]

	for channel in channels:
		print 'Channel : {}'.format(channel)
		for variable in config.variables:
			tau_value = get_tau_value(config, channel, variable) 
			tau_value = 0
			print 'Variable : {}'.format(variable)

			d_unfoldingInfo = get_unfolding_info(config, samples_and_files_to_compare, variable, channel)
			d_unfoldingInfo = unfolding_shenanigans(d_unfoldingInfo, tau_value)
			plotting_shenanigans(d_unfoldingInfo, variable, channel)
	return


def get_unfolding_info(config, samples_and_files_to_compare, variable, channel):
	'''
	Return Unfolding Information
	'''

	for sample, unfolding_info in samples_and_files_to_compare.iteritems():
		if sample == 'data': continue

		h_truth, h_smeared, h_response, h_fakes = get_unfold_histogram_tuple(
			inputfile=unfolding_info['file'],
			variable=variable,
			channel=channel,
			centre_of_mass=config.centre_of_mass_energy,
			ttbar_xsection=config.ttbar_xsection,
			luminosity=config.luminosity,
			load_fakes=True,
			visiblePS=True,
		)

		samples_and_files_to_compare[sample]['response'] 	= h_response
		samples_and_files_to_compare[sample]['smeared'] 	= h_smeared
		samples_and_files_to_compare[sample]['truth'] 		= h_truth
		samples_and_files_to_compare[sample]['fakes'] 		= h_fakes
		samples_and_files_to_compare[sample]['integral'] 	= asrootpy(h_response.ProjectionY()).integral(0,-1)

	edges = reco_bin_edges_vis[variable]
	samples_and_files_to_compare['data']['hist'] 		= asrootpy(value_error_tuplelist_to_hist( read_tuple_from_file( samples_and_files_to_compare['data']['file'].format(var=variable, ch=channel))['TTJet'] , edges ))
	samples_and_files_to_compare['data']['integral'] 	= asrootpy(samples_and_files_to_compare['data']['hist'])
	return samples_and_files_to_compare


def get_tau_value(config, channel, variable):
    if channel == 'electron':
        return config.tau_values_electron[variable]
    if channel == 'muon':
        return config.tau_values_muon[variable]


def unfolding_shenanigans(d_unfoldingInfo, tau):
	'''
	Add unfolded histograms to d_unfoldingInfo
	'''
	d_unfoldingInfo = unfold_with_powhegpythia(d_unfoldingInfo, tau)
	d_unfoldingInfo = unfold_with_alternate(d_unfoldingInfo, tau)
	d_unfoldingInfo = unfold_the_data(d_unfoldingInfo, tau)
	return d_unfoldingInfo

def unfold_with_powhegpythia(d_unfoldingInfo, tau):
	'''
	Unfold the alternate smeared distributions with powhegPythia response matrix
	'''
	for sample, unfolding_info in d_unfoldingInfo.iteritems():
		if sample == 'data': continue
		scale = unfolding_info['integral'] / d_unfoldingInfo['powhegPythia']['integral']
		scale = 1

		smeared = unfolding_info['smeared']
		truth = unfolding_info['truth']
		response = d_unfoldingInfo['powhegPythia']['response']
		response.Scale(scale)

		# Unfold, and set 'data' to 'smeared' 
		unfolding = Unfolding( 
			smeared,
			truth,
			smeared,
		    response, 
		    fakes=None,
		    method='TUnfold',
		    tau=tau
		)
		d_unfoldingInfo[sample]['unfolded_with_powhegPythia'] = unfolding.unfold()
	return d_unfoldingInfo

def unfold_with_alternate(d_unfoldingInfo, tau):
	'''
	Unfold the smeared powhegPythia distribution with alternate response matrices
	'''
	for sample, unfolding_info in d_unfoldingInfo.iteritems():
		if sample == 'data': continue
		scale = d_unfoldingInfo['powhegPythia']['integral'] / unfolding_info['integral']
		scale = 1

		smeared = d_unfoldingInfo['powhegPythia']['smeared']
		truth = d_unfoldingInfo['powhegPythia']['truth']
		response = unfolding_info['response']
		response.Scale(scale)

		# Unfold, and set 'data' to 'smeared' 
		unfolding = Unfolding( 
			smeared,
			truth,
			smeared,
		    response, 
		    fakes=None,
		    method='TUnfold',
		    tau=tau
		)
		d_unfoldingInfo['powhegPythia']['unfolded_with_{}'.format(sample)] = unfolding.unfold()
	return d_unfoldingInfo

def unfold_the_data(d_unfoldingInfo, tau):
	'''
	Unfold the smeared powhegPythia distribution with alternate response matrices
	'''
	for sample, unfolding_info in d_unfoldingInfo.iteritems():
		if sample == 'data': continue

		h_data_no_fakes = removeFakes( d_unfoldingInfo[sample]['smeared'], d_unfoldingInfo[sample]['fakes'],  d_unfoldingInfo['data']['hist'] )
		d_unfoldingInfo['data']['hist'] = asrootpy(h_data_no_fakes)
		d_unfoldingInfo['data']['integral'] = d_unfoldingInfo['data']['hist'].integral(0,-1)

		scale = d_unfoldingInfo['data']['integral'] / unfolding_info['integral']
		scale = 1

		smeared = unfolding_info['smeared']
		truth = unfolding_info['truth']
		response = unfolding_info['response']
		response.Scale(scale)

		# Unfold, and set 'data' to 'smeared' 
		unfolding = Unfolding( 
			d_unfoldingInfo['data']['hist'],
			truth,
			smeared,
		    response, 
		    fakes=None,
		    method='TUnfold',
		    tau=tau
		)
		d_unfoldingInfo['data']['unfolded_with_{}'.format(sample)] = unfolding.unfold()
	return d_unfoldingInfo


def plotting_shenanigans(d_unfoldingInfo, variable, channel):
	output_folder = 'plots/unfolding/closure_with_generators/'

	hp = Histogram_properties()
	v_latex = latex_labels.variables_latex[variable]
	unit = ''
	if variable in ['HT', 'ST', 'MET', 'WPT', 'lepton_pt']: unit = ' [GeV]'

	hp.x_axis_title = v_latex + unit
	hp.y_axis_title = 'Number of unfolded events'  
	hp.title = 'Generator tests for {variable}'.format(variable=v_latex)

	# Alternate with PowhegPythia
	hp.name = 'nUnfEvent_{channel}_{variable}_Alternate_Unfolded_With_PowhegPythia'.format(
	    channel=channel,
	    variable=variable,
	)

	models = OrderedDict()
	models['powhegPythia'] 				=   d_unfoldingInfo['powhegPythia']['truth']
	models['amcatnloPythia'] 			=   d_unfoldingInfo['amcatnloPythia']['truth']
	models['powhegHerwig'] 				=   d_unfoldingInfo['powhegHerwig']['truth']
	models['madgraph'] 					=   d_unfoldingInfo['madgraph']['truth']

	measurements = OrderedDict()
	measurements['powhegPythia'] 		=   d_unfoldingInfo['powhegPythia']['unfolded_with_powhegPythia']
	measurements['amcatnloPythia'] 		=   d_unfoldingInfo['amcatnloPythia']['unfolded_with_powhegPythia']
	measurements['powhegHerwig'] 		=   d_unfoldingInfo['powhegHerwig']['unfolded_with_powhegPythia']
	measurements['madgraph'] 			=   d_unfoldingInfo['madgraph']['unfolded_with_powhegPythia']

	compare_measurements(
		models = models,
		measurements = measurements,
		show_measurement_errors=True,
		histogram_properties=hp,
		save_folder=output_folder,
		line_styles_for_models=["solid", "solid", "solid", "solid"],
		save_as=['pdf'],
		match_models_to_measurements = True
	)

	# PowhegPythia with Alternate
	hp.name = 'nUnfEvent_{channel}_{variable}_PowhegPythia_Unfolded_With_Alternate'.format(
	    channel=channel,
	    variable=variable,
	)

	models = OrderedDict()
	models['powhegPythia'] 				=   d_unfoldingInfo['powhegPythia']['truth']
	models['amcatnloPythia'] 			=   d_unfoldingInfo['amcatnloPythia']['truth']
	models['powhegHerwig'] 				=   d_unfoldingInfo['powhegHerwig']['truth']
	models['madgraph'] 					=   d_unfoldingInfo['madgraph']['truth']

	measurements = OrderedDict()
	measurements['powhegPythia'] 		=   d_unfoldingInfo['powhegPythia']['unfolded_with_powhegPythia']
	measurements['amcatnloPythia'] 		=   d_unfoldingInfo['powhegPythia']['unfolded_with_amcatnloPythia']
	measurements['powhegHerwig'] 		=   d_unfoldingInfo['powhegPythia']['unfolded_with_powhegHerwig']
	measurements['madgraph'] 			=   d_unfoldingInfo['powhegPythia']['unfolded_with_madgraph']

	compare_measurements(
		models = models,
		measurements = measurements,
		show_measurement_errors=True,
		histogram_properties=hp,
		save_folder=output_folder,
		line_styles_for_models=["solid", "dashed", "dashed", "dashed"],
		save_as=['pdf'],
		match_models_to_measurements = True
	)

	# Data with Generators
	hp.name = 'nUnfEvent_{channel}_{variable}_Data_Unfolded_With_Generator'.format(
	    channel=channel,
	    variable=variable,
	)

	models = OrderedDict()
	models['powhegPythia'] 				=   d_unfoldingInfo['powhegPythia']['truth']
	models['amcatnloPythia'] 			=   d_unfoldingInfo['amcatnloPythia']['truth']
	models['powhegHerwig'] 				=   d_unfoldingInfo['powhegHerwig']['truth']
	models['madgraph'] 					=   d_unfoldingInfo['madgraph']['truth']

	measurements = OrderedDict()
	measurements['powhegPythia'] 		=   d_unfoldingInfo['data']['unfolded_with_powhegPythia']
	measurements['amcatnloPythia'] 		=   d_unfoldingInfo['data']['unfolded_with_amcatnloPythia']
	measurements['powhegHerwig'] 		=   d_unfoldingInfo['data']['unfolded_with_powhegHerwig']
	measurements['madgraph'] 			=   d_unfoldingInfo['data']['unfolded_with_madgraph']

	compare_measurements(
		models = models,
		measurements = measurements,
		show_measurement_errors=True,
		histogram_properties=hp,
		save_folder=output_folder,
		line_styles_for_models=["solid", "solid", "solid", "solid"],
		save_as=['pdf'],
		match_models_to_measurements = True
	)
	return

if __name__ == '__main__':
	main()
