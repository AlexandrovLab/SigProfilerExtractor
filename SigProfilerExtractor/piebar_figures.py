import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.patches as mpatches
import pylab
import numpy as np
import matplotlib.font_manager as font_manager

# Example input to program:
#cosmic_contribs = [['SBS1', '11.28%'], ['SBS5', '42.58%'], ['SBS15', '16.90%'], ['SBS20', '29.24%']]

# Color Palettes - we will want to look at this in more details
colors = [[237/255,201/255,81/255],[79/255, 55/255, 45/255], [235/255, 104/255, 65/255], [0/255, 160/255, 176/255], [204/255, 42/255, 54/255]]

# context input determines output Path


# parse a given string of format "12.34%" to ".1234"
def contribs_str_to_num(num):
    num_percentile = float(num[:len(num)-1])
    num_dec = num_percentile/100

    return num_dec


# The format of the input being passed
#cosmic_contribs[i][0] = signature name
#cosmic_contribs[i][1] = contribution
def gen_bar_png(cosmic_contribs,output_dir):
    figData = plt.figure(figsize=(.2, 5.5))
    stacked_bar = figData.add_subplot(111)
    xpos = 0
    bottom = 0
    width = .2

    # Generate the stacked percentile bar plot
    for j in range(len(cosmic_contribs)):
        height = contribs_str_to_num(cosmic_contribs[j][1])
        stacked_bar.bar(xpos, height, width, bottom=bottom, color=colors[j], edgecolor='none')
        ypos = bottom + stacked_bar.patches[j].get_height() / 2
        bottom += height

    plt.gca().set_axis_off()
    plt.margins(0,0)
    plt.savefig(output_dir+"_bar.png", bbox_inches="tight", transparent=True, format="png")




# Needs colors and labels (ie, [['SBS1', '45.4%'], ['SBS2', '44.6%'] and ['red', 'blue'])
# Two lists, one with the decomposition, second the color to associate with it
def gen_legend_png(cosmic_contribs, output_dir):
    plt.clf()
    plt.figure(figsize=(3,3))
    patch_arr = []

    for i in range(0, len(cosmic_contribs)):
        tmp_lab = cosmic_contribs[i][0] + " (" + cosmic_contribs[i][1] + ")"
        tmp_patch = mpatches.Patch(color=colors[i], label=tmp_lab)
        patch_arr.append(tmp_patch)

    font = font_manager.FontProperties(family='Arial', style='normal', size='96', weight='normal')
    plt.legend(handles=patch_arr, loc='center', prop=font)
    plt.gca().set_axis_off()
    plt.savefig(output_dir + "_legend.png", bbox_inches="tight", transparent=True, format="png")


def gen_figures(cosmic_contribs, output_dir):
    gen_bar_png(cosmic_contribs,output_dir)
    gen_legend_png(cosmic_contribs, output_dir)

#example
#gen_figures([['SBS1', '11.28%'], ['SBS5', '42.58%'], ['SBS15', '16.90%'], ['SBS20', '29.24%']])
