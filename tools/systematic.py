from __future__ import division, print_function
from tools.file_utilities import read_data_from_JSON, write_data_to_JSON
from tools.Calculation import combine_errors_in_quadrature
from config import XSectionConfig
from copy import deepcopy
from math import sqrt
import numpy as np

def write_normalised_xsection_measurement(options, measurement, measurement_unfolded, summary = '' ):
    '''
    Writes the list of normalised measurements and normalised unfolded measurements of the form: 
    [Central Value, Lower Systemtic, Upper Systematic] to a json. Different combinations of 
    systematic uncertainty are stored as different json by appending different 'summary'
    '''
    path_to_JSON=options['path_to_JSON']
    method=options['method']
    channel=options['channel']

    output_file = '{path_to_JSON}/TEST/central/normalised_xsection_{channel}_{method}_with_errors.txt'
    output_file = output_file.format(
                    path_to_JSON = path_to_JSON,
                    channel = channel,
                    method = method,
                    )
    
    if not summary == '':
        output_file = output_file.replace( 'with_errors', summary + '_errors' )
    
    output = {'TTJet_measured':measurement, 'TTJet_unfolded': measurement_unfolded}
    
    write_data_to_JSON( output, output_file )

    print("Data written to : ", output_file)

def append_PDF_uncertainties(all_systematics):
    '''
    Replace 'PDF' entry in list of all systematics with actual PDF variations
    '''
    variation = []
    for index in xrange (0, 100):
        variation.append('PDFWeights_'+str(index))
    all_systematics['PDF'] = variation
    return all_systematics

def read_normalised_xsection_measurement(options, category):
    '''
    Returns the normalised measurement and normalised unfolded measurement for 
    the file associated with the variable under study
    '''
    variable=options['variable']
    variables_no_met=options['variables_no_met']
    met_specific_systematics=options['met_specific_systematics']
    path_to_JSON=options['path_to_JSON']
    method=options['method']
    channel=options['channel']
    filename = '{path}/{category}/normalised_xsection_{channel}_{method}.txt'
    # Disregarding Met Uncertainties if variable does not use MET
    if (category in met_specific_systematics) and (variable in variables_no_met):
        filename = filename.format(
            path = path_to_JSON,
            channel = channel,
            category = 'central',
            method = method,
        )
    else:
        filename = filename.format(
            path = path_to_JSON,
            channel = channel,
            category = category,
            method = method
        )
    normalised_xsection = read_data_from_JSON( filename )
    measurement = normalised_xsection['TTJet_measured']#should this be measured without fakes???
    measurement_unfolded = normalised_xsection['TTJet_unfolded']
    return measurement, measurement_unfolded  

def read_normalised_xsection_systematics(options, variation, is_multiple_sources=False):
    '''
    Returns the list of normalised measurements and normalised unfolded measurements 
    for each systematic category
    variation: current systematic (BJet, PDF etc)
    is_multiple_sources: is svariation composed of multiple sources? (PDF : PDFWeight1, ...)
    '''
    systematics = {}
    systematics_unf = {}

    if is_multiple_sources:
        for source in variation:
            src, src_unf = read_normalised_xsection_measurement(options, source)
            systematics[source] = src
            systematics_unf[source] = src_unf  
    else:

        upper_variation = variation[0]
        lower_variation = variation[1]
        systematic_up, systematic_up_unf = read_normalised_xsection_measurement(options, upper_variation)
        systematic_down, systematic_down_unf = read_normalised_xsection_measurement(options, lower_variation)
        systematics['lower'] = systematic_down
        systematics['upper'] = systematic_up        
        systematics_unf['lower'] = systematic_down_unf
        systematics_unf['upper'] = systematic_up_unf

    systematics = get_systematic_measured_values_only(options, systematics)   
    systematics_unf = get_systematic_measured_values_only(options, systematics_unf)   
    return systematics, systematics_unf

def get_systematic_measured_values_only(options, syst_x_secs_and_errs):
    '''
    When retreiveing systematic uncertainties they come in form [Value, Error]. 
    This strips the Errors.
    '''
    number_of_bins=options['number_of_bins']
    for variation, x_sec_and_error in syst_x_secs_and_errs.iteritems():
        for bin_i in xrange(number_of_bins):
            x_sec = x_sec_and_error[bin_i][0]
            x_sec_and_error[bin_i] = x_sec
    return syst_x_secs_and_errs


def get_normalised_cross_sections(options, list_of_systematics):
    '''
    Gets the normalised cross sections in the form:
    List of Systematics in Group = {[[central], [upper], [lower]]}
        [central] = [[x sec, unc], [x sec, unc]...[x sec, unc]] For N Bins
        [upper]   = [[x sec], [x sec]...[x sec]] For N Bins
        [lower]   = [[x sec], [x sec]...[x sec]] For N Bins
    '''
    # Copy structure of the systmatic lists in order to replace systematic name with their normailsed x sections
    normalised_systematic_uncertainty_x_sections = deepcopy(list_of_systematics)
    unfolded_normalised_systematic_uncertainty_x_sections = deepcopy(list_of_systematics)

    central_measurement, central_measurement_unfolded = read_normalised_xsection_measurement(options, 'central')
    # central measurement
    normalised_systematic_uncertainty_x_sections['central'] = central_measurement
    unfolded_normalised_systematic_uncertainty_x_sections['central'] = central_measurement

    for systematic, variation in list_of_systematics.iteritems():
        if (systematic == 'PDF'):
            syst_unc_x_sec, unf_syst_unc_x_sec = read_normalised_xsection_systematics(options, variation, is_multiple_sources=True)
            normalised_systematic_uncertainty_x_sections[systematic] = []
            unfolded_normalised_systematic_uncertainty_x_sections[systematic] = []
            for weight, vals in syst_unc_x_sec.iteritems():
                normalised_systematic_uncertainty_x_sections[systematic].append([weight, vals])
            for weight, vals in unf_syst_unc_x_sec.iteritems():
                unfolded_normalised_systematic_uncertainty_x_sections[systematic].append([weight, vals])

        elif (systematic == 'TTJets_envelope'):
            syst_unc_x_sec, unf_syst_unc_x_sec = read_normalised_xsection_systematics(options, variation, is_multiple_sources=True)
            env_lower, env_upper = get_scale_envelope(options, syst_unc_x_sec)
            unf_env_lower, unf_env_upper = get_scale_envelope(options, unf_syst_unc_x_sec)
            normalised_systematic_uncertainty_x_sections[systematic] = [
                env_lower,
                env_upper,
            ]
            unfolded_normalised_systematic_uncertainty_x_sections[systematic] = [
                unf_env_lower,
                unf_env_upper,
            ]
        # elif (systematic == 'TTJets_alphaS'): continue
        else :
            syst_unc_x_sec, unf_syst_unc_x_sec = read_normalised_xsection_systematics(options, variation)
            normalised_systematic_uncertainty_x_sections[systematic] = [
                syst_unc_x_sec['upper'], 
                syst_unc_x_sec['lower'],
            ]
            unfolded_normalised_systematic_uncertainty_x_sections[systematic] = [
                unf_syst_unc_x_sec['upper'], 
                unf_syst_unc_x_sec['lower'],
            ]

    return normalised_systematic_uncertainty_x_sections, unfolded_normalised_systematic_uncertainty_x_sections


def get_scale_envelope(options, d_syst):
    '''
    Calculate the scale envelope for the renormalisation/factorisation/combined systematic uncertainties
    Do we need to store syst-central or just syst?
    '''
    down_variations = []
    up_variations = []
    envelope_up = []
    envelope_down = []

    # Separate into up/down scale variations
    for scale_variation in d_syst:
        if 'down' in scale_variation:
            down_variations.append(d_syst[scale_variation])
        elif 'up' in scale_variation:
            up_variations.append(d_syst[scale_variation])

    # find min/max
    for v1, v2, v3 in zip (up_variations[0], up_variations[1], up_variations[2]):
        envelope_up.append(max(v1, v2, v3))
    for v1, v2, v3 in zip (down_variations[0], down_variations[1], down_variations[2]):
        envelope_down.append(min(v1, v2, v3))
    return envelope_down, envelope_up




def calculate_total_PDFuncertainty(options, central_measurement, pdf_uncertainty_values):
    '''
    Calculates the appropriate lower and upper PDF uncertainty
    @param central_measurement: measurement from central PDF weight
    @param pdf_uncertainty_values: dictionary of measurements with different PDF weights; 
                                    format {PDFWeights_%d: measurement}
    '''
    number_of_bins = options['number_of_bins']
    pdf_min = []
    pdf_max = []
    # print_dictionary("Central",central_measurement)
    # print_dictionary("PDF",pdf_uncertainty_values)

    # split PDF uncertainties into downwards (negative) and upwards (positive) components
    for bin_i in xrange(number_of_bins):
        negative = []
        positive = []
        central = central_measurement[bin_i][0]
        for pdf_variation in pdf_uncertainty_values:
            # Get PDF weight index from its name
            index = int(pdf_variation[0].rsplit('_', 1)[1])
            pdf_uncertainty = pdf_variation[1][bin_i]
            if index % 2 == 0:  # even == negative
                negative.append(pdf_uncertainty)
            else:
                positive.append(pdf_uncertainty)

        # now to calculate the RMS (sigma)
        rms_up, rms_down = 0, 0
        for n,p in zip(negative, positive):
            rms_up += (p-central)**2
            rms_down += (n-central)**2   
        pdf_max.append( sqrt( rms_up / (len(positive)-1) ))
        pdf_min.append( sqrt( rms_down / (len(negative)-1) ))

        pdf_sym = max(pdf_min, pdf_max)

    # return pdf_min, pdf_max  
    return pdf_sym  

def get_symmetrised_systematic_uncertainty(norm_syst_unc_x_secs, options):
    '''
    Gets the symmetrised normalised cross sections uncertainties in the form:
    Group of Systematics : { List of Systematics in Group : [[central], [symmetrised uncertainty], [signed uncertainty]]}
        [central]                   = [[x sec, unc], [x sec, unc]...[x sec, unc]] For N Bins
        [symmetrised uncertainty]   = [sym unc, sym unc...sym unc] For N Bins
        [signed uncertainty]        = [sign, sign...sign] For N Bins
    '''
    normalised_x_sections_with_symmetrised_systematics = deepcopy(norm_syst_unc_x_secs)
    central_measurement = norm_syst_unc_x_secs['central']
    for systematic, variation in norm_syst_unc_x_secs.iteritems():
        if (systematic == 'PDF'):
            # Replace all PDF weights with full PDF combination
            pdf_sym = calculate_total_PDFuncertainty(
                options, 
                central_measurement, 
                variation,
            )
            # TODO Find signs etc... i.e. do proper covariance for PDF
            normalised_x_sections_with_symmetrised_systematics[systematic] = [
                pdf_sym, 
                [0]*len(central_measurement),
            ]  

        else:
            upper_measurement = variation[0]
            lower_measurement = variation[1]

            isTopMassSystematic = True if systematic == 'TTJets_mass' else False

            symmetrised_uncertainties, signed_uncertainties = get_symmetrised_errors(
                central_measurement, 
                upper_measurement, 
                lower_measurement, 
                options, 
                isTopMassSystematic,
            )

            normalised_x_sections_with_symmetrised_systematics['central'] = norm_syst_unc_x_secs['central']
            normalised_x_sections_with_symmetrised_systematics[systematic] = [
                symmetrised_uncertainties, 
                signed_uncertainties,
            ]         



    # Combine PDF with alphaS variations
    alphaS = normalised_x_sections_with_symmetrised_systematics['TTJets_alphaS'][0]
    pdf = normalised_x_sections_with_symmetrised_systematics['PDF'][0]
    pdf_tot = [combine_errors_in_quadrature([e1, e2]) for e1, e2 in zip(alphaS, pdf)]
    normalised_x_sections_with_symmetrised_systematics['PDF'][0] = pdf_tot
    # TODO combine the signs....

    # Now alphaS is combined with pdfs dont need it in dictionary anymore.
    del normalised_x_sections_with_symmetrised_systematics['TTJets_alphaS']
    return normalised_x_sections_with_symmetrised_systematics           

def get_symmetrised_errors(central_measurement, upper_measurement, lower_measurement, options, isTopMassSystematic=False ):
    '''
    Returns the symmetric error in each bin for a specific systematic and also the sign of the systematic.
    Sign is used for calculating the covariance matrices. Returns of the form:
    [Symmetric Uncertainties],[Signed Uncertainties]
        [symmetrised uncertainty]   = [sym unc, sym unc...sym unc] For N Bins
        [signed uncertainty]        = [sign, sign...sign] For N Bins
    '''
    number_of_bins = len(central_measurement)
    symm_uncerts = []
    sign_uncerts = []


    for c, u, l in zip(central_measurement, upper_measurement, lower_measurement):
        central = c[0] # Getting the measurement, not the error [xsec, unc]
        upper = u
        lower = l

        upper_uncertainty = abs(central - upper)
        lower_uncertainty = abs(central - lower)

        if isTopMassSystematic:
            upper_uncertainty, lower_uncertainty = scaleTopMassSystematic( upper_uncertainty, lower_uncertainty, options['topMasses'], options['topMassUncertainty'] )

        symmetrised_uncertainty = max(upper_uncertainty, lower_uncertainty)
        #  Getting the sign of the uncertainty
        if u == l:
            sign = np.sign( upper - central )
        else:
            sign = np.sign( upper - lower )

        symm_uncerts.append(symmetrised_uncertainty)
        sign_uncerts.append(sign)
    return symm_uncerts, sign_uncerts

def scaleTopMassSystematic( upper_uncertainty, lower_uncertainty, topMasses, topMassUncertainty ):
    lowMassDifference = topMasses[1] - topMasses[0]
    upMassDifference = topMasses[2] - topMasses[1]

    upper_uncertainty *= topMassUncertainty / upMassDifference
    lower_uncertainty *= topMassUncertainty / lowMassDifference

    return upper_uncertainty, lower_uncertainty

def get_sign(central, upper, lower, upper_variation, lower_variation):
    '''
    Returns the sign of the uncertainty - i.e. was the upper variation bigger than the lower variation
    Returns 0 if the systematic is symmetrical
    '''
    if ((upper_variation > lower_variation) and (upper > central) ): return 1
    elif ((upper_variation > lower_variation) and (upper < central) ): return -1
    elif ((lower_variation > upper_variation) and (lower > central) ): return -1
    elif ((lower_variation > upper_variation) and (lower < central) ): return 1
    else: return 0

def get_measurement_with_total_systematic_uncertainty(options, x_sec_with_symmetrised_systematics):
    '''
    Returns the measurement with the total symmetrised systematic uncertainties of the form:
    Group of Systematics : [central value][central uncertainty][symmetrised systematic uncertainty+][symmetrised systematic uncertainty-]
        [central value]                         = [[x sec, unc], [x sec, unc],... [x sec, unc]] For N Bins
        [symmetrised systematic uncertainty+]   = [central unc, central unc,... central unc] For N Bins
        [symmetrised systematic uncertainty-]   = [sys unc, sys unc,... sys unc] For N Bins
    '''
    number_of_bins = options['number_of_bins']
    sys_unc = 0
    measurement_with_total_uncertainty = []
    for bin_i in range( 0, number_of_bins ):
        central = x_sec_with_symmetrised_systematics['central'][bin_i] # Still [Value, Error]
        for systematic, measurement in x_sec_with_symmetrised_systematics.iteritems():
            if (systematic == 'central'): continue
            sys_unc += measurement[0][bin_i]**2
        measurement_with_total_uncertainty.append( [central[0], sqrt(sys_unc), sqrt(sys_unc)] )
    return measurement_with_total_uncertainty

def print_dictionary(title, dictionary_to_print):
    '''
    Prints dictionaries in a nicer form
    '''
    import pprint
    pp = pprint.PrettyPrinter(indent=4)
    print ()
    print ('-'*100)
    print (title)
    print ('-'*100)
    pp.pprint(dictionary_to_print)
    print ()
    return

def generate_covariance_matrices(options, x_sec_with_symmetrised_systematics):
    '''
    Iterates through each group of systematics generating and plotting the covariance matrix of the systematic
    '''
    number_of_bins=options['number_of_bins']

    for systematic, measurement in x_sec_with_symmetrised_systematics.iteritems():
        if systematic == 'central': continue
        covariance_matrix, correlation_matrix = generate_covariance_matrix(number_of_bins, systematic, measurement)
        make_covariance_plot(options, systematic, covariance_matrix)
        make_covariance_plot(options, systematic, correlation_matrix, label='Correlation')
    return

def generate_covariance_matrix(number_of_bins, systematic, measurement):
    '''
    Uses the symmetrised normalised cross sections uncertainties in the form:
    Group of Systematics : { List of Systematics in Group : [[central], [symmetrised uncertainty], [signed uncertainty]]}
    Covariance_ij = (Sign_i*Unc_i) * (Sign_j*Unc_j)
    Variance_ii = (Unc_i) * (Unc_i)
    Returns the matrix in the form [[Bin_i, Bin_j], Cov_ij]
    '''
    covariance_matrix = []
    correlation_matrix = []
    for bin_i in xrange(number_of_bins):
        for bin_j in xrange(number_of_bins):
            if (bin_j < bin_i): continue
            uncertainty_i = measurement[0][bin_i]
            uncertainty_j = measurement[0][bin_j]
            sign_i = measurement[1][bin_i]
            sign_j = measurement[1][bin_j]
            cov_ij = (sign_i*uncertainty_i)*(sign_j*uncertainty_j)
            cor_ij = ( sign_i * sign_j )
            # Bins when plotting Histogram start from 1 not 0
            bin_and_value_covariance = [[bin_i+1, bin_j+1], cov_ij]
            covariance_matrix.append(bin_and_value_covariance)
            bin_and_value_correlation = [[bin_i+1, bin_j+1], cor_ij]
            correlation_matrix.append(bin_and_value_correlation)
            if not bin_i == bin_j:
                bin_and_value = ([[bin_j+1, bin_i+1], cov_ij])
                covariance_matrix.append(bin_and_value)
                bin_and_value_correlation = [[bin_j+1, bin_i+1], cor_ij]
                correlation_matrix.append(bin_and_value_correlation)
    return covariance_matrix, correlation_matrix

def make_covariance_plot( options, systematic, matrix, label='Covariance' ):
    '''
    Take the matrix in list form and bin edges in list form to create a TH2F of the covariance matrix
    Saves to plots/covariance_matrices/{PhaseSpace}/{Channel}/{Variable}/
    '''
    from config.variable_binning import bin_edges_vis
    from ROOT import TH2F, TCanvas, TPad, gROOT, gStyle
    from array import array
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)

    variable = options['variable']
    channel = options['channel']
    covariance_matrix_output_path = options['covariance_matrix_output_path']

    x_binning = array ( 'f' , bin_edges_vis[variable] )
    y_binning = array ( 'f', bin_edges_vis[variable] )
    n_xbins = len( x_binning ) - 1
    n_ybins = len( y_binning ) - 1

    hist = TH2F('name', 'title', n_xbins, x_binning, n_ybins, y_binning )
    set_bin_value = hist.SetBinContent
    for entry in matrix:
        bin_i = entry[0][0]
        bin_j = entry[0][1]
        cov_ij = entry[1]
        set_bin_value( bin_i, bin_j, cov_ij )
    # Easy access to .pngs 
    canvas = TCanvas("canvas_name","canvas_title", 0, 0, 1000, 800)
    hist.SetTitle(systematic+" "+label+" matrix for "+variable+" in channel "+channel+" ; Bin_i; Bin_j")
    hist.Draw("colz text")
    canvas.Update()
    canvas.SaveAs(covariance_matrix_output_path+systematic+'_'+label+'_matrix.pdf')