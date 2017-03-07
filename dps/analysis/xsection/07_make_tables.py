from argparse import ArgumentParser
from dps.config.xsection import XSectionConfig
from dps.utils.pandas_utilities import file_to_df
from dps.config.latex_labels import variables_latex
from dps.config.variable_binning import bin_edges_vis
from dps.utils.file_utilities import make_folder_if_not_exists

def makeLatexTable( xsections, outputPath, variable, crossSectionType ):
	'''
	Generate and write the Latex table for the cross sections
	'''
    #########################################################################################################
    ### Table Header
    #########################################################################################################
	latexHeader = '\\begin{table}\n'
	latexHeader += '\t\centering\n'
	latexHeader += '\t\\begin{tabular}{|cccc|}\n'
	latexHeader += '\t\t\hline\n'

	labelHeader = '\t\t'
	labelHeader += '\\textbf{{{var}}} \t& \\textbf{{Central}} \t& \\textbf{{Statistical Uncertainty}} \t& \\textbf{{Systematic Uncertainty}}\t \\\\ \n'.format(var=variables_latex[variable])
	# labelHeader += '\t\t\hline\n'

	fullTable = latexHeader
	fullTable += labelHeader

	channel = ''
	for ch in range(len(xsections)):
	    #########################################################################################################
	    ### Table Sub Header
	    #########################################################################################################
		if   ch == 0: channel = 'electron'
		elif ch == 1: channel = 'muon'
		elif ch == 2: channel = 'combined'

		subHeader = '\t\t\\textbf{{{ch}}} \t& \t& \t& \t \\\\ \n'.format(ch=channel)
		subHeader += '\t\t\hline\n'
		fullTable += subHeader

	    #########################################################################################################
	    ### Table Content
	    #########################################################################################################
		for bin in range (len(bin_edges_vis[variable])-1):
			line_for_bin = '\t\t{edge_down}-{edge_up} \t& {val} \t& {stat} \t& {sys} \t \\\\ \n'.format(
				edge_down = bin_edges_vis[variable][bin],
				edge_up = bin_edges_vis[variable][bin+1],
				val = xsections[ch]['central'][bin],
				stat = xsections[ch]['statistical'][bin],
				sys = xsections[ch]['systematic'][bin],
			)
			fullTable += line_for_bin

    #########################################################################################################
    ### Table Footer
    #########################################################################################################
	tableFooter = '\t\end{tabular}\n'
	tableFooter += '\t\caption{{Results of the {type} differential cross sections with respect to {var}.}}\n'.format(
		type 	=crossSectionType,
		var 	=variable,
	)
	tableFooter += '\t\label{{tb:xsection_{type}_{var}}}\n'.format(
		type 	=crossSectionType,
		var 	=variable,
	)	
	tableFooter += '\\end{table}\n'

	fullTable += tableFooter

    #########################################################################################################
    ### Write Table
    #########################################################################################################
	make_folder_if_not_exists(outputPath)
	file_template = outputPath + '/{var}_XSections.tex'.format(var=variable)
	output_file = open(file_template, 'w')
	output_file.write(fullTable)
	output_file.close()

def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument( "-p", "--path", 
        dest    = "path", 
        default = 'data/normalisation/background_subtraction/',
        help    = "set path to files containing dataframes" 
    )
    parser.add_argument( '--visiblePS', 
        dest    = "visiblePS", 
        action  = "store_true",
        help    = "Unfold to visible phase space" 
    )
    parser.add_argument( '--outputTablePath','-o', 
        dest    = "outputTablePath",
        default = 'tables/xsections/',
        help    = "Output path for chi2 tables" 
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
	args = parse_arguments()

	measurement_config      = XSectionConfig( 13 )
	visiblePS               = args.visiblePS
	outputTablePath 		= args.outputTablePath

	phase_space = 'FullPS'
	if visiblePS:
		phase_space = 'VisiblePS'

	channels = [
		'electron', 
		'muon', 
		'combined', 
		# 'combinedBeforeUnfolding',
	]
	unc_type = [
		'normalised',
		'absolute',
	]

	input_file_template 	= "xsection_{type}_{channel}_TUnfold_summary_absolute.txt"
	path_to_input_template  = '{path}/{com}TeV/{variable}/{phase_space}/central/'
	path_to_output_template = '{path}/{crossSectionType}/'
	for utype in unc_type:
		for variable in measurement_config.variables:
			print "Writing the {type} {var} cross sections to Latex Tables".format(type = utype, var=variable)
			path_to_output = path_to_output_template.format(
				path=outputTablePath, 
				crossSectionType=utype,
			)
			path_to_input = path_to_input_template.format(
			    path = args.path, 
			    com = 13,
			    variable = variable,
			    phase_space = phase_space,
			)

			# Read cross sections and write tables
			xsections = [file_to_df(path_to_input+input_file_template.format(type = utype, channel = ch)) for ch in channels]
			makeLatexTable( xsections=xsections, outputPath=path_to_output, variable=variable, crossSectionType=utype )
