#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 20:41:28 2020

@author: mishugeb
"""


import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def plotTMB(inputDF , output, scale="genome", Yrange="adapt"):
    if type(scale) == int:
        scale = scale
    elif scale == "genome":
        scale  = 2800
    elif scale == "exome":
        scale = 55
    else:
        print("Please input valid scale values: \"exome\", \"genome\" or a numeric value")
        return
    inputDF.columns = ['Types', 'Mut_burden']
    df=inputDF[inputDF["Mut_burden"]!=0]
    df['log10BURDENpMB'] = df.apply(lambda row: np.log10(row.Mut_burden/scale), axis = 1)
    groups = df.groupby(["Types"])
    means = groups.mean()["log10BURDENpMB"].sort_values(ascending=True)
    names = groups.mean()["log10BURDENpMB"].sort_values(ascending=True).index
    counts = groups.count()["log10BURDENpMB"][names]
    ngroups = groups.ngroups
    if Yrange == "adapt":
        ymax = math.ceil(df['log10BURDENpMB'].max())
        ymin = math.floor(df['log10BURDENpMB'].min())
    elif Yrange == "cancer":
        ymax = 3
        ymin = -3
    elif type(Yrange) == list:
        print("Yrange is a list")
        ymax = int(math.log10(Yrange[1]))
        ymin = int(math.log10(Yrange[0]))
    else:
        print("Please input valid scale values: \"adapt\", \"cancer\" or a list of two power of 10 numbers")
        return
    #plotting
    fig, ax = plt.subplots(figsize=(1+0.4*ngroups,(ymax-ymin)*0.5+3))
    plt.xlim(0,2*ngroups)
    plt.ylim(ymin,ymax)
    yticks_loc = range(ymin,ymax+1,1)
    plt.yticks(yticks_loc,list(map((lambda x: 10**x), list(yticks_loc)))) 
    plt.xticks(np.arange(1, 2*ngroups+1, step = 2),counts) 
    plt.tick_params(axis = 'both', which = 'both',length = 0)
    plt.hlines(yticks_loc,0,2*ngroups,colors = 'black',linestyles = "dashed",linewidth = 0.5,zorder = 1)
    for i in range(0,ngroups,2):
        greystart = [(i)*2,ymin]
        rectangle = mpatches.Rectangle(greystart, 2, ymax-ymin, color = "lightgrey",zorder = 0)
        ax.add_patch(rectangle)
    for i in range(0,ngroups,1):
        X_start = i*2+0.1
        X_end = i*2+2-0.1
        rg = 1.8
        y_values = groups.get_group(names[i])["log10BURDENpMB"].sort_values(ascending = True).values.tolist()
        x_values = list(np.linspace(start = X_start, stop = X_end, num = counts[i]))
        plt.scatter(x_values,y_values,color = "black",s=1.5)
        plt.hlines(means[i],X_start,X_end,colors='red',zorder=2)
    axes2 = ax.twiny()
    plt.tick_params(axis = 'both', which = 'both',length = 0)
    plt.xticks(np.arange(1, 2*ngroups+1, step = 2),names,rotation = -35,ha = 'right');
    fig.subplots_adjust(top=0.8,left=1/(1+0.4*ngroups), right=1-(0.3/(1+0.4*ngroups)),bottom=0.05)
    plt.savefig(output) 






