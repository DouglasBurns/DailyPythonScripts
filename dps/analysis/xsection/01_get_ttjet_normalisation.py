from __future__ import division
from argparse import ArgumentParser
from dps.utils.logger import log
from dps.config.xsection import XSectionConfig
from dps.utils.file_utilities import get_files_in_path, read_data_from_JSON, make_folder_if_not_exists, get_filename_without_extension
from dps.utils.measurement import Measurement
from dps.utils.ROOT_utils import set_root_defaults
from dps.utils.pandas_utilities import file_to_df, dict_to_df, df_to_file

# define logger for this module
mylog = log["01b_get_ttjet_normalisation"]

def main():
    '''
    1 - Read Config file for normalisation measurement
    2 - Run measurement
    3 - Combine measurement before unfolding
    '''
    results = {}

    # config file template
    input_template          = 'config/measurements/background_subtraction/{com}TeV/{ch}/{var}/{ps}/'

    ps = 'FullPS'
    if args.visiblePS:
        ps = 'VisiblePS'

    channels = [
        'electron', 
        'muon',
    ]

    for ch in channels:
        for var in measurement_config.variables:
            if args.variable and args.variable != var: continue
            qcd_transfer_factor = {}

            # Create measurement_filepath
            measurement_filepath = input_template.format(
                com = args.CoM,
                ch = ch,
                var = var,
                ps = ps,
            )
            
            # Get all config files in measurement_filepath
            measurement_files = get_files_in_path(measurement_filepath, file_ending='.json')

            for f in sorted(measurement_files):
                if args.test:
                    if 'central' not in f: continue
                print('Processing file ' + f)
                sample = get_filename_without_extension(f)

                # Read in Measurement JSON
                config = read_data_from_JSON(f)

                # Create Measurement Class using JSON
                if 'electron' in ch:
                    electron_measurement = Measurement(config)
                    if args.forControlPlots:
                        electron_measurement.output_folder = electron_measurement.output_folder.replace('data', 'data_for_01')
                    electron_measurement.calculate_normalisation()
                    electron_measurement.save(ps)
                    qcd_transfer_factor[sample] = [electron_measurement.t_factor]

                elif 'muon' in ch:
                    muon_measurement = Measurement(config)
                    if args.forControlPlots:
                        muon_measurement.output_folder = muon_measurement.output_folder.replace('data', 'data_for_01')

                    muon_measurement.calculate_normalisation()
                    muon_measurement.save(ps)
                    qcd_transfer_factor[sample] = [muon_measurement.t_factor]

            output_folder = electron_measurement.output_folder.format(
                com = args.CoM,  var = var,
                ps  = ps,   cat = 'central',
            )
            store_transfer_factor(qcd_transfer_factor, output_folder, ch)
    return

def store_transfer_factor(tf, output_file, channel):
    make_folder_if_not_exists(output_file)
    f = output_file+'table_of_transfer_factors_'+channel+'.txt'
    df = dict_to_df(tf)
    df_to_file(f, df)
    return

def parse_arguments():
    parser = ArgumentParser(__doc__)
    parser.add_argument("-v", "--variable", dest="variable", default=None,
                            help="set the variable to analyse (MET, HT, ST, MT, WPT). Default is MET.")
    parser.add_argument("-c", "--centre-of-mass-energy", dest="CoM", default=13, type=int,
                            help="set the centre of mass energy for analysis. Default = 13 [TeV]")
    parser.add_argument('--visiblePS', dest="visiblePS", action="store_true",
                            help="Unfold to visible phase space")
    parser.add_argument('--forControlPlots', dest="forControlPlots", action="store_true",
                            help="Store 01 for use in make_controlPlots_from01")
    parser.add_argument('--test', dest="test", action="store_true",
                            help="test on central only.")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    set_root_defaults()
    args = parse_arguments()
    measurement_config = XSectionConfig(args.CoM)
    main()




