"""
Copyright 2020-2022 Francesco Di Lauro. All Rights Reserved.
See LICENCE file for details
"""

import os as os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import seaborn as sns

def my_plot_configs():
    plt.style.use('seaborn-bright')
    plt.rcParams["figure.frameon"] = False
    plt.rcParams['ps.useafm'] = True
    plt.rcParams['pdf.use14corefonts'] = True
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Helvetica'
    plt.rcParams['axes.labelweight'] = 'bold'

my_plot_configs()
today = pd.to_datetime('today')

#whole india
url = 'https://api.covid19india.org/csv/latest/case_time_series.csv'
fname = 'COVID19_timeseries_India_' + today.strftime('%m%d') + '.csv'

#state by state
#url = 'https://api.covid19india.org/csv/latest/state_wise_daily.csv'
#fname = 'COVID19_timeseries_India_state_wise' + today.strftime('%m%d') + '.csv'


date_fields = ["Date", "Date_YMD"]
full_data = pd.read_csv(url, parse_dates=date_fields)

data_folder = 'India/'

full_data.to_csv(os.path.join(data_folder,fname), index = False)

ifPlot = False

if ifPlot:
    fig = plt.figure()
    plt.plot(full_data["Date_YMD"][::3], full_data["DL"][1::3], '-', color='gray', lw=3)
    plt.xlabel('Date')
    plt.ylabel('Daily Confirmed')
    ax = plt.gca()
    # date_form = DateFormatter("%m-%d")
    # ax.xaxis.set_major_formatter(date_form)
    sns.despine()
    fig.show()
    fname = 'COVID19_India_Cases' + today.strftime('%m%d') + '.pdf'
    fig.savefig(os.path.join(plot_folder, fname), format='pdf', transparent=True)
