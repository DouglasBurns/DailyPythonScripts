import gc
from argparse import ArgumentParser

from dps.utils.Unfolding import Unfolding, get_unfold_histogram_tuple
from dps.utils.file_utilities import make_folder_if_not_exists
from dps.utils.ROOT_utils import set_root_defaults

from dps.config import CMS, latex_labels
from dps.config.xsection import XSectionConfig

import ROOT as r
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import plot

from rootpy import asrootpy
from rootpy.io import File, Directory
from rootpy.plotting import Hist, Hist2D, Canvas, Legend
import rootpy.plotting.root2matplotlib as rplt

mpl.use('agg')
from dps.utils.latex import setup_matplotlib
setup_matplotlib()

def get_tau_value(config, channel, variable):
    if channel == 'electron':
        return config.tau_values_electron[variable]
    if channel == 'muon':
        return config.tau_values_muon[variable]
    if channel == 'combined':
        return config.tau_values_combined[variable]

def makeToyResponse( response ):
	for binx in range(0,response.GetNbinsX() + 2 ):
		for biny in range(0,response.GetNbinsY() + 2 ):
			binContent = response.GetBinContent( binx, biny )
			binError = response.GetBinError( binx, biny )
			if binContent <= 0:
				newBinContent = 0
			else:
				newBinContent = np.random.poisson(binContent)
			if newBinContent < 0: newBinContent = 0

			newBinError = binError

			newBinError = np.sqrt( newBinContent )
			response.SetBinContent( binx, biny, newBinContent )
			response.SetBinError( binx, biny, newBinError)
	return response

def parse_arguments():
    parser = ArgumentParser(__doc__)
    parser.add_argument( "-v", "--variable", dest = "variable", default = None,
                      help = "set the variable to analyse (MET, HT, ST, MT)" )
    parser.add_argument( "--channel", dest = "channel", default = None,
                      help = "set the channel (electron, muon, combined)" )

    args = parser.parse_args()
    return args

def main():
	args = parse_arguments()

	config = XSectionConfig(13)
	method = 'TUnfold'

	files_for_response = [
		File(config.unfolding_central, 'read')
	]

	files_for_toys = [
		File(config.unfolding_central, 'read')
	]

	channels = [
		'electron', 
		'muon',
		# 'combined',
	]

	for channel in channels:
		for variable in config.variables:
			if args.variable and variable != args.variable: continue
			if args.channel and channel != args.channel: continue

			tau_value = get_tau_value(config, channel, variable)
			print "Running pull test on {} in the {} channel. Tau value retreived is {}".format(variable, channel, tau_value)
			
			pullHistogram = None

			for file_for_response in files_for_response:

				_, _, h_response, _ = get_unfold_histogram_tuple(
				    inputfile=file_for_response,
				    variable=variable,
				    channel=channel,
				    centre_of_mass=config.centre_of_mass_energy,
				    ttbar_xsection=config.ttbar_xsection,
				    luminosity=config.luminosity,
				    load_fakes=False,
				    visiblePS=True,
				)

				if pullHistogram is None:
					pullHistogram = Hist2D( h_response.GetNbinsY(), 1, h_response.GetNbinsY()+1, 1000, -10, 10 )
					pullHistogram.SetDirectory(0)

				for file_for_toys in files_for_toys:

					_, _, h_response_for_toys, _ = get_unfold_histogram_tuple(
					    inputfile=file_for_toys,
					    variable=variable,
					    channel=channel,
					    centre_of_mass=config.centre_of_mass_energy,
					    ttbar_xsection=config.ttbar_xsection,
					    luminosity=config.luminosity,
					    load_fakes=False,
					    visiblePS=True,
					)

					n=5000
					for i in range(0,n):

						if i % 500 == 0: print 'Toy number :',i

						toy_response = makeToyResponse( h_response_for_toys.Clone() )
						toy_measured = asrootpy(toy_response.ProjectionX('px',1))
						toy_truth = asrootpy(h_response_for_toys.ProjectionY())

						toy_response_unfolding = makeToyResponse( h_response.Clone() )
						toy_response_unfolding.Scale( toy_response.integral(overflow=True) / toy_response_unfolding.integral(overflow=True) )

						# Unfold toy data with independent toy response
						unfolding = Unfolding( toy_measured,
							toy_truth, toy_measured, toy_response_unfolding, None,
							method='TUnfold', tau=tau_value)

						unfolded_results = unfolding.unfold()

						cov, cor, mc_cov = unfolding.get_covariance_matrix()
						total_statistical_covariance = cov + mc_cov
						for i in range(0,total_statistical_covariance.shape[0] ):
							unfolded_results.SetBinError(i+1, np.sqrt( total_statistical_covariance[i,i] ) )


						for bin in range(1,unfolded_results.GetNbinsX() + 1 ):
							diff = unfolded_results.GetBinContent(bin) - toy_truth.GetBinContent(bin)
							pull = diff / unfolded_results.GetBinError( bin )
							pullHistogram.Fill( bin, pull )

			plot_pull_distribution(pullHistogram, variable, channel, n)
			plot_pull_quantities(pullHistogram, variable, channel, n)
	return

def plot_pull_distribution(pullHistogram, var, ch, n):
	pullDistribution = asrootpy(pullHistogram.ProjectionY())
	pullDistribution.Fit('gaus', 'WWSQ')
	fit_pull = pullDistribution.GetFunction('gaus')

	A 			= fit_pull.GetParameter(0)
	Aerror 		= fit_pull.GetParError(0)
	mean 		= fit_pull.GetParameter(1)
	meanError 	= fit_pull.GetParError(1)
	sigma 		= fit_pull.GetParameter(2)
	sigmaError 	= fit_pull.GetParError(2)

	fig = plt.figure(figsize=(16, 16), dpi=200, facecolor='white')
	axes = fig.add_subplot(1,1,1)
	pullDistribution.SetMarkerSize(CMS.data_marker_size)

	rplt.errorbar(pullDistribution, xerr=True, emptybins=True, axes=axes, zorder=1)
    # *4 for a very smooth curve
	x = np.linspace(fit_pull.GetXmin(), fit_pull.GetXmax(), fit_pull.GetNpx() * 4)
	function_data = np.frompyfunc(fit_pull.Eval, 1, 1)
	plot(x, function_data(x), axes=axes, color='red', linewidth=4, zorder=2)

	axes.set_xlim([-5,5])

	plt.xlabel('$\\frac{N^{\mathrm{unfolded}} - N^{\mathrm{true}}}{\sigma}$', CMS.x_axis_title)
	plt.ylabel('$N^{\\text{pseudo experiments}}$', CMS.y_axis_title)
	plt.tick_params(**CMS.axis_label_major)
	plt.tick_params(**CMS.axis_label_minor)

	title_template = 'Total Pull Distribution for {variable}\n'
	title_template += '$\sqrt{{s}}$ = 13 TeV, {channel}, $N^{{\\text{{pseudo experiments}}}}_{{\\text{{bin}}}}$ = {exp}'
	title = title_template.format(
	    variable=latex_labels.variables_latex[var],
	    channel=latex_labels.channel_latex[ch],
	    exp=n
	)
	plt.title(title, CMS.title)

	text_template = 'mean = ${mean} \pm  {mean_error}$\n'
	text_template += '$\sigma = {sigma} \pm  {sigma_error}$'
	text = text_template.format(
	    mean=round(mean, 2),
	    mean_error=round(meanError, 3),
	    sigma=round(sigma, 2),
	    sigma_error=round(sigmaError, 3),
	)
	axes.text(0.6, 0.8, text,
		verticalalignment='bottom', horizontalalignment='left',
		transform=axes.transAxes,
		color='black', fontsize=40, bbox=dict(facecolor='white', edgecolor='none', alpha=0.5)
	)
	plt.tight_layout()

	outputDir = 'plots/unfolding/pulls/new/'
	make_folder_if_not_exists(outputDir)
	outputName = '{dir}/{variable}_{channel}_PullDist.pdf'.format( dir = outputDir, variable = var, channel = ch)
	plt.savefig(outputName)
	fig.clf()
	plt.close()
	gc.collect()
	return


def plot_pull_quantities(pullHistogram, var, ch, n):

	nBins = pullHistogram.GetNbinsX()

	plots = r.TObjArray()
	pullHistogram.FitSlicesY(0,0,-1,0,'QNR',plots)
	means, widths = None, None

	for p in plots:
		if p.GetName()[-2:] == '_1':
			means = asrootpy(p)
		elif p.GetName()[-2:] == '_2':
			widths = asrootpy(p)

	fig = plt.figure(figsize=(16, 16), dpi=200, facecolor='white')
	axes = fig.add_subplot(1,1,1)

	means.SetMarkerColor(2)
	means.SetMarkerSize(4)
	means.SetLineColor(2)
	widths.SetMarkerColor(4)
	widths.SetMarkerSize(4)
	widths.SetLineColor(4)

	plt.axhline(y=1, linewidth=2, color='blue', linestyle='--')
	plt.axhline(y=0, linewidth=2, color='red', linestyle='--')
	rplt.errorbar(means, xerr=True, yerr=None, emptybins=True, axes=axes, capsize=0, elinewidth=4, zorder=3, label="Pull Mean")
	rplt.errorbar(widths, xerr=True, yerr=None, emptybins=True, axes=axes, capsize=0, elinewidth=4, zorder=4, label="Pull Width")
	axes.set_ylim([-2,2])
	axes.set_xlim([1,nBins+1])

	plt.xlabel(latex_labels.variables_latex[var]+' Bin Number', CMS.x_axis_title)
	plt.tick_params(**CMS.axis_label_major)
	plt.tick_params(**CMS.axis_label_minor)

	title_template = 'Pull Distribution Mean and Width for {variable}\n'
	title_template += '$\sqrt{{s}}$ = 13 TeV, {channel}, $N^{{\\text{{pseudo experiments}}}}_{{\\text{{bin}}}}$ = {exp}'
	title = title_template.format(
	    variable=latex_labels.variables_latex[var],
	    channel=latex_labels.channel_latex[ch],
	    exp=n
	)
	plt.title(title, CMS.title)
	leg = plt.legend(loc='lower right', prop={'size':50}, numpoints=1)
	leg.draw_frame(False)	

	plt.tight_layout()
 
	outputDir = 'plots/unfolding/pulls/new/'
	make_folder_if_not_exists(outputDir)
	outputName = '{dir}/{variable}_{channel}.pdf'.format( dir = outputDir, variable = var, channel = ch)
	plt.savefig(outputName)
	fig.clf()
	plt.close()
	gc.collect()
	return

if __name__ == '__main__':
    set_root_defaults( set_batch = True, msg_ignore_level = 3001 )
    main()