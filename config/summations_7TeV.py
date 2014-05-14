'''
Created on 19 Jan 2013

@author: kreczko
'''

# 7 TeV
electron_qcd_samples = [ 'QCD_Pt-20to30_BCtoE',
                 'QCD_Pt-30to80_BCtoE',
                 'QCD_Pt-80to170_BCtoE',
                 'QCD_Pt-20to30_EMEnriched',
                 'QCD_Pt-30to80_EMEnriched',
                 'QCD_Pt-80to170_EMEnriched',
                 'GJets_HT-40To100',
                 'GJets_HT-100To200',
                 'GJets_HT-200']
singleTop_samples = [ 'T_tW-channel',
             'T_t-channel',
             'T_s-channel',
             'Tbar_tW-channel',
             'Tbar_t-channel',
             'Tbar_s-channel']
wplusjets_samples = [ 'W1Jet', 'W2Jets', 'W3Jets', 'W4Jets']
vplusjets_samples = wplusjets_samples
vplusjets_samples.append('DYJetsToLL')
diboson_samples = [ 'WWtoAnything', 'WZtoAnything', 'ZZtoAnything']
signal_samples = [ 'TTJet', 'SingleTop']

wplusjets_matchingup_samples = [ 'WJets-matchingup' ]
dyplusjets_matchingup_samples = [ 'ZJets-matchingup' ]
vplusjets_matchingup_samples = wplusjets_matchingup_samples + dyplusjets_matchingup_samples

wplusjets_matchingdown_samples = [ 'WJets-matchingdown' ]
dyplusjets_matchingdown_samples = [ 'ZJets-matchingdown' ]
vplusjets_matchingdown_samples = wplusjets_matchingdown_samples + dyplusjets_matchingdown_samples

wplusjets_scaledown_samples = [ 'WJets-scaledown' ]
dyplusjets_scaledown_samples = [ 'ZJets-scaledown' ]
vplusjets_scaledown_samples = wplusjets_scaledown_samples + dyplusjets_scaledown_samples

wplusjets_scaleup_samples = [ 'WJets-scaleup' ]
dyplusjets_scaleup_samples = [ 'ZJets-scaleup' ]
vplusjets_scaleup_samples = wplusjets_scaleup_samples + dyplusjets_scaleup_samples

ttjets_unfolding_samples = ['TTJets']
ttjets_mcatnlo_unfolding_samples = ['TTJets']
ttjets_powheg_unfolding_samples = ['TTJets']
ttjets_matchingup_unfolding_samples = ['TTJets-matchingup']
ttjets_matchingdown_unfolding_samples = ['TTJets-matchingdown']
ttjets_scaleup_unfolding_samples = ['TTJets-scaleup']
ttjets_scaledown_unfolding_samples = ['TTJets-scaledown']

sample_summations = {
                  'QCD-Electron':electron_qcd_samples,
                  'SingleTop' : singleTop_samples,
                  'WJets' : wplusjets_samples,
                  'VJets' : vplusjets_samples,
#                  'DiBoson': diboson_samples,
                  'Signal': signal_samples,
                  'VJets-matchingup' : vplusjets_matchingup_samples,
                  'VJets-matchingdown' : vplusjets_matchingdown_samples,
                  'VJets-scaledown' : vplusjets_scaledown_samples,
                  'VJets-scaleup' : vplusjets_scaleup_samples,
                  'unfolding-merged' : ttjets_unfolding_samples,
                  'unfolding-TTJets-7TeV-mcatnlo' : ttjets_mcatnlo_unfolding_samples,
                  'unfolding-TTJets-7TeV-powheg' : ttjets_powheg_unfolding_samples,
                  'unfolding-TTJets-7TeV-matchingup' : ttjets_matchingup_unfolding_samples,
                  'unfolding-TTJets-7TeV-matchingdown' : ttjets_matchingdown_unfolding_samples,
                  'unfolding-TTJets-7TeV-scaleup' : ttjets_scaleup_unfolding_samples,
                  'unfolding-TTJets-7TeV-scaledown' : ttjets_scaledown_unfolding_samples,
                  }
