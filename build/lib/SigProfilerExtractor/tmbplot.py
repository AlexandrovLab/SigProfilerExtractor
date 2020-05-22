import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.axis as axis

def prepend(list, str): 
    str += '{0}'
    list = [str.format(i) for i in list] 
    return(list) 

def plotTMB(inputDF, scale, Yrange = "adapt", cutoff = 1, output = "TMB_plot.pdf"):
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
    df=inputDF[inputDF["Mut_burden"] >= cutoff]
    df['log10BURDENpMB'] = df.apply(lambda row: np.log10(row.Mut_burden/scale), axis = 1)
    groups = df.groupby(["Types"])
    means = groups.mean()["log10BURDENpMB"].sort_values(ascending=True)
    names = groups.mean()["log10BURDENpMB"].sort_values(ascending=True).index
    counts = groups.count()["log10BURDENpMB"][names]
    ngroups = groups.ngroups
    
    #second row of bottom label
    input_groups = inputDF.groupby(["Types"])
    input_counts = input_groups.count()["Mut_burden"][names]
    list1 = counts.to_list()
    list2 = input_counts.to_list()
    str1 = ''
    list3 = prepend(list1, str1)
    str2 = '\n'
    list4 = prepend(list2, str2)
    result = [None]*(len(list3)+len(list4))
    result[::2] = list3
    result[1::2] = list4
    tick_labels = result
    new_labels = [ ''.join(x) for x in zip(tick_labels[0::2], tick_labels[1::2]) ]
    
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
    fig_width = 2.3 + 0.4 * ngroups
    fig_length = (ymax - ymin) * 0.7 + 3
    fig, ax = plt.subplots(figsize=(fig_width, fig_length))
    plt.xlim(0,2*ngroups)
    plt.ylim(ymin,ymax)
    yticks_loc = range(ymin,ymax+1,1)
    plt.yticks(yticks_loc,list(map((lambda x: 10**x), list(yticks_loc)))) 
    plt.xticks(np.arange(1, 2*ngroups+1, step = 2),new_labels) 
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
        plt.hlines(means[i], X_start, X_end, colors='red', zorder=2)
        plt.text((2.2 + i * 0.4) / fig_width , 0.85 / fig_length , "___",  horizontalalignment='center',transform=plt.gcf().transFigure)
        #plt.hlines(ymin+0.2, bottom_start, bottom_end, colors='blue', zorder=2)
    axes2 = ax.twiny()
    plt.tick_params(axis = 'both', which = 'both',length = 0)
    plt.xticks(np.arange(1, 2*ngroups+1, step = 2),names,rotation = -35,ha = 'right');
    fig.subplots_adjust(top = ((ymax - ymin) * 0.7 + 1) / fig_length, bottom = 1 / fig_length, left = 2 / fig_width, right=1 - 0.3 / fig_width)
    plt.text(1.7 / fig_width, 0.2 / fig_length, "*Excluding samples with less than %d mutations" % cutoff, transform=plt.gcf().transFigure)
    plt.savefig(output)