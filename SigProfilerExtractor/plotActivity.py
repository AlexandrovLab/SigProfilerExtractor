import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import json
from matplotlib import rc
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages

colors = ['tab:pink',
 'tab:orange',
 'tab:purple',
 'tab:olive',
 'tab:brown',
 'deeppink',
 'tab:red',
 'rebeccapurple',
 'blue',
 'magenta',
 'dodgerblue',
 'salmon',
 'lime',
 'cyan',
 'lightgray',
 'silver',
 'darkgrey',
 'gray',
 'dimgrey',
 'black']

color_code = pd.DataFrame(np.array([['SBS1', 'forestgreen'],
       ['SBS5', 'limegreen'],
       ['SBS2', 'lightseagreen'],
       ['SBS13', 'turquoise'],
       ['SBS7a', 'darkgoldenrod'],
       ['SBS7b', 'goldenrod'],
       ['SBS7c', 'gold'],
       ['SBS7d', 'khaki'],
       ['SBS10a', 'indianred'],
       ['SBS10b', 'lightcoral'],
       ['SBS17a', 'deepskyblue'],
       ['SBS17b', 'skyblue']]),columns=['signature', 'color'])


########## If read from file ############
#color_code = pd.read_table('./color_dic.txt')
#with open('./additional_colors.txt', 'r') as f:
#    colors = f.read().splitlines()
##########################################

def plotActivity(activity_file, output_file, bin_size = 50):
    size = int(bin_size)
    inputDF = pd.read_table(activity_file,index_col=0)
    inputDF = inputDF.loc[:, (inputDF != 0).any(axis=0)]
    inputDF["sum"]=inputDF.sum(axis=1)
    inputDF.sort_values(by="sum",inplace=True,ascending=False)
    inputDF.drop(columns="sum",inplace=True)
    list_of_dfs = [inputDF.iloc[i:i+size,:] for i in range(0, len(inputDF),size)]
    all_sig = list(inputDF.columns.values)
    s1 = color_code[color_code['signature'].isin(all_sig)]["signature"].tolist()
    s2=[x for x in all_sig if x not in s1]
    signature_list=s1+s2
    color_list = color_code[color_code['signature'].isin(all_sig)]["color"].tolist()+colors[:len(s2)]
    pp = PdfPages(output_file)
    for j in range(0,len(list_of_dfs)):
        nsamples = len(list_of_dfs[j])
        plot_length = nsamples/50*10
        figure_length = plot_length + 5
        Lmargin = 1.5/figure_length
        names = list_of_dfs[j].index.values.tolist()
        x_pos = list(range(0, len(names)))
        barWidth=1
        sig_activity_list=[]
        for i in range(0,len(signature_list)):
            sig_activity_list.append(list_of_dfs[j][signature_list[i]].tolist())
        plot = plt.figure(figsize=(figure_length,7))
        plt.rc('axes', edgecolor = 'lightgray')
        panel1 = plt.axes([Lmargin, 0.25, plot_length / figure_length , 0.6])
        bottom_bars = []
        plot_list = []
        plt.xlim([-0.5, len(names)-0.5])
        for i in range(0,len(signature_list)):
            if i == 0:
                plot_list.append(plt.bar(x_pos, sig_activity_list[i], color = color_list[i], edgecolor = 'white', width = barWidth))
            elif i == 1:
                bottom_bars = sig_activity_list[0]
                plot_list.append(plt.bar(x_pos, sig_activity_list[i], bottom = bottom_bars, color = color_list[i], edgecolor = 'white', width = barWidth))
            else:
                bottom_bars = np.add(bottom_bars,sig_activity_list[i-1]).tolist()
                plot_list.append(plt.bar(x_pos, sig_activity_list[i], bottom = bottom_bars, color = color_list[i], edgecolor = 'white', width = barWidth))
        plt.xticks(x_pos, names, rotation = 40,ha = "right")
        plt.legend(plot_list, signature_list, fontsize="small", bbox_to_anchor=(1 + 0.5 / plot_length, 1), loc='upper left', borderaxespad=0.)
        pp.savefig(plot) 
        plt.close()
    pp.close()
