from argparse import ArgumentParser
from dps.config.latex_labels import samples_latex, channel_latex, variables_latex
from dps.config.variable_binning import fit_variable_bin_edges, control_plots_bins
from dps.config.histogram_colours import histogram_colours as colours
from dps.config.xsection import XSectionConfig
from dps.utils.file_utilities import make_folder_if_not_exists
from dps.utils.plotting import make_data_mc_comparison_plot, Histogram_properties
from dps.utils.hist_utilities import prepare_histograms, clean_control_region
from dps.utils.ROOT_utils import get_histograms_from_trees, set_root_defaults
from dps.utils.pandas_utilities import dict_to_df, df_to_file
from dps.utils.latex import setup_matplotlib
from uncertainties import ufloat
import pandas as pd 
import glob

# latex, font, etc
setup_matplotlib()
title_template = '%.1f fb$^{-1}$ (%d TeV)'

def getUnmergedDirectory( f ) :
    baseDir = f.split('combined')[0]
    sampleName = f.split('combined')[-1].strip('/').split('_tree.root')[0]
    print baseDir
    print sampleName
    new_f = baseDir + '/' + sampleName + '/analysis_central_job_*/*root'
    return new_f
