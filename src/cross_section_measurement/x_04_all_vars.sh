#!/bin/bash
python 04_make_plots_matplotlib.py &> logs/met_plot.log &
python 04_make_plots_matplotlib.py -v HT &> logs/ht_plot.log &
python 04_make_plots_matplotlib.py -v ST &> logs/st_plot.log &
python 04_make_plots_matplotlib.py -v MT &> logs/mt_plot.log &
python 04_make_plots_matplotlib.py -v WPT &> logs/wpt_plot.log &