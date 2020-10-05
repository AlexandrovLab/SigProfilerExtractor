#!/usr/bin/env python3

#Author: Erik Bergstrom

# Modifications:
# SBS-96, SBS-1536, DBS-78, and ID-83 output to .png extension
# Input parameter matrix_path has been updated to be a pandas dataframe
# modifications made by Mark Barnes

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mplpatches
import re
import os
import sys
import argparse
from collections import OrderedDict
import pandas as pd
import numpy as np
# imports for saving plots to memory
import io
from PIL import Image

def plotSBS(matrix_path, output_path, project, plot_type, percentage=False, custom_text_upper=None, custom_text_middle=None, custom_text_bottom=None):
	plot_custom_text = False
	sig_probs = False
	pcawg = False
	buff_list = dict()
	if plot_type == '96':
		first_line=matrix_path.iloc[0,:]
		if first_line[0][1] == ">":
			pcawg = True
		if first_line[0][5] != "]" and first_line[0][1] != ">":
			sys.exit("The matrix does not match the correct SBS96 format. Please check you formatting and rerun this plotting function.")

		# pp = PdfPages(output_path + 'SBS_96_plots_' + project + '.pdf')

		mutations = OrderedDict()
		total_count = []

		if pcawg:
			samples=matrix_path.columns[:]
			samples = samples[2:]
		else:
			#samples = first_line.strip().split("\t")
			samples=matrix_path.columns[:]
			samples = samples[1:]

		for sample in samples:
			mutations[sample] = OrderedDict()
			mutations[sample]['C>A'] = OrderedDict()
			mutations[sample]['C>G'] = OrderedDict()
			mutations[sample]['C>T'] = OrderedDict()
			mutations[sample]['T>A'] = OrderedDict()
			mutations[sample]['T>C'] = OrderedDict()
			mutations[sample]['T>G'] = OrderedDict()

		for lines_tmp in range(0,matrix_path.shape[0]):
			if pcawg:
				line = matrix_path.iloc[lines_tmp,:]
				#line = lines.strip().split(",")
				mut_type = line[0]
				nuc = line[1][0] + "[" + mut_type + "]" + line[1][2]
				sample_index = 2
			else:
				line = matrix_path.iloc[lines_tmp,:]
				#line = lines.strip().split()
				nuc = line[0]
				mut_type = line[0][2:5]
				sample_index = 1

			for sample in samples:
				if percentage:
					mutCount = float(line[sample_index])
					if mutCount < 1 and mutCount > 0:
						sig_probs = True
				else:
					mutCount = int(line[sample_index])
				mutations[sample][mut_type][nuc] = mutCount
				sample_index += 1

		sample_count = 0
		for sample in mutations.keys():
			total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
			plt.rcParams['axes.linewidth'] = 2
			plot1 = plt.figure(figsize=(43.93,9.92))
			plt.rc('axes', edgecolor='lightgray')
			panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
			xlabels = []

			x = 0.4
			ymax = 0
			colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
			i = 0
			for key in mutations[sample]:
				for seq in mutations[sample][key]:
					xlabels.append(seq[0]+seq[2]+seq[6])
					if percentage:
						if total_count > 0:
							plt.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
							if mutations[sample][key][seq]/total_count*100 > ymax:
								ymax = mutations[sample][key][seq]/total_count*100
					else:
						plt.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
						if mutations[sample][key][seq] > ymax:
								ymax = mutations[sample][key][seq]
					x += 1
				i += 1

			x = .043
			y3 = .87
			y = int(ymax*1.25)
			y2 = y+2
			for i in range(0, 6, 1):
				panel1.add_patch(plt.Rectangle((x,y3), .15, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
				x += .159

			yText = y3 + .06
			plt.text(.1, yText, 'C>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.255, yText, 'C>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.415, yText, 'C>T', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.575, yText, 'T>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.735, yText, 'T>C', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.89, yText, 'T>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

			if y <= 4:
				y += 4

			while y%4 != 0:
				y += 1
			ytick_offest = int(y/4)


			if percentage:
				ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
				ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%",
						  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
			else:
				ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
				ylabels= [0, ytick_offest, ytick_offest*2,
						  ytick_offest*3, ytick_offest*4]

			labs = np.arange(0.375,96.375,1)

			if not percentage:
				ylabels = ['{:,}'.format(int(x)) for x in ylabels]
				if len(ylabels[-1]) > 3:
					ylabels_temp = []
					if len(ylabels[-1]) > 7:
						for label in ylabels:
							if len(label) > 7:
								ylabels_temp.append(label[0:-8] + "m")
							elif len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)

					else:
						for label in ylabels:
							if len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)
					ylabels = ylabels_temp

			panel1.set_xlim([0, 96])
			panel1.set_ylim([0, y])
			panel1.set_xticks(labs)
			panel1.set_yticks(ylabs)
			count = 0
			m = 0
			for i in range (0, 96, 1):
				plt.text(i/101 + .0415, .02, xlabels[i][0], fontsize=30, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				plt.text(i/101 + .0415, .044, xlabels[i][1], fontsize=30, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
				plt.text(i/101 + .0415, .071, xlabels[i][2], fontsize=30, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				count += 1
				if count == 16:
					count = 0
					m += 1

			if sig_probs:
				plt.text(0.045, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
			else:
				plt.text(0.045, 0.75, sample + ": " + "{:,}".format(int(total_count)) + " subs", fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)



			custom_text_upper_plot = ''
			try:
				custom_text_upper[sample_count]
			except:
				custom_text_upper = False
			try:
				custom_text_middle[sample_count]
			except:
				custom_text_middle = False
			try:
				custom_text_bottom[sample_count]
			except:
				custom_text_bottom = False

			if custom_text_upper:
				plot_custom_text = True
				if len(custom_text_upper[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False
			if custom_text_middle:
				if len(custom_text_middle[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False

			if plot_custom_text:
				x_pos_custom = 0.98
				if custom_text_upper and custom_text_middle:
					custom_text_upper_plot = custom_text_upper[sample_count] + "\n" + custom_text_middle[sample_count]
					if custom_text_bottom:
						custom_text_upper_plot += "\n" + custom_text_bottom[sample_count]

				if custom_text_upper and not custom_text_middle:
					custom_text_upper_plot = custom_text_upper[sample_count]
					panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

				elif custom_text_upper and custom_text_middle:
					if not custom_text_bottom:
						panel1.text(x_pos_custom, 0.72, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')
					else:
						panel1.text(x_pos_custom, 0.68, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

				elif not custom_text_upper and custom_text_middle:
					custom_text_upper_plot = custom_text_middle[sample_count]
					panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')



			panel1.set_yticklabels(ylabels, fontsize=30)
			plt.gca().yaxis.grid(True)
			plt.gca().grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
			panel1.set_xlabel('')
			panel1.set_ylabel('')

			if percentage:
				plt.ylabel("Percentage of Single Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')
			else:
				plt.ylabel("Number of Single Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')



			panel1.tick_params(axis='both',which='both',\
							   bottom=False, labelbottom=False,\
							   left=True, labelleft=True,\
							   right=True, labelright=False,\
							   top=False, labeltop=False,\
							   direction='in', length=25, colors='lightgray', width=2)


			[i.set_color("black") for i in plt.gca().get_yticklabels()]
			# Save matplot lib as a BytesIO, add to list, and return this list
			buffer = io.BytesIO()
			plt.savefig(buffer, format="png")
			buff_list[sample]=buffer
			plt.close()
			sample_count += 1
		return buff_list

	elif plot_type == '192' or plot_type == '96SB' or plot_type == '384':
		first_line=matrix_path.iloc[0,:]
		if first_line[0][6] == ">" or first_line[0][3] == ">":
			pcawg = True
		if first_line[0][7] != "]" and first_line[0][6] != ">" and first_line[0][3] != ">":
			sys.exit("The matrix does not match the correct SBS192 format. Please check you formatting and rerun this plotting function.")

		mutations = OrderedDict()

		if pcawg:
			samples=matrix_path.columns[:]
			samples = samples[3:]
			samples = [x.replace('"','') for x in samples]
		else:
			samples=matrix_path.columns[:]
			samples = samples[1:]

		for sample in samples:
			mutations[sample] = OrderedDict()
			mutations[sample]['C>A'] = OrderedDict()
			mutations[sample]['C>G'] = OrderedDict()
			mutations[sample]['C>T'] = OrderedDict()
			mutations[sample]['T>A'] = OrderedDict()
			mutations[sample]['T>C'] = OrderedDict()
			mutations[sample]['T>G'] = OrderedDict()


		for lines_tmp in range(0,matrix_path.shape[0]):
			if pcawg:
				line = matrix_path.iloc[lines_tmp,:]
				line = [x.replace('"','') for x in line]
				nuc = line[2][0] + "[" + line[1] + "]" + line[2][2]
				bias = line[0][0]
			else:
				line = matrix_path.iloc[lines_tmp,:]
				nuc = line[0][2:]
				bias = line[0][0]
			if bias == 'N' or bias == 'B':
				continue
			else:
				if pcawg:
					mut_type = line[1]
					sample_index = 3
				else:
					mut_type = line[0][4:7]
					sample_index = 1

				for sample in samples:
					if percentage:
						mutCount = float(line[sample_index])
						if mutCount < 1 and mutCount > 0:
							sig_probs = True
					else:
						mutCount = int(line[sample_index])
					if nuc not in mutations[sample][mut_type].keys():
						mutations[sample][mut_type][nuc] = [0,0]
					if bias == 'T':
						mutations[sample][mut_type][nuc][0] = mutCount
					else:
						mutations[sample][mut_type][nuc][1] = mutCount
					sample_index += 1

		sample_count = 0
		for sample in mutations.keys():
			total_count = sum(sum(sum(tsb) for tsb in nuc.values()) for nuc in mutations[sample].values())
			plt.rcParams['axes.linewidth'] = 2
			plot1 = plt.figure(figsize=(43.93,9.92))
			plt.rc('axes', edgecolor='lightgray')
			panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
			xlabels = []

			x = 0.7
			ymax = 0
			colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
			i = 0
			for key in mutations[sample]:
				for seq in mutations[sample][key]:
					xlabels.append(seq[0]+seq[2]+seq[6])
					if percentage:
						if total_count > 0:
							trans = plt.bar(x, mutations[sample][key][seq][0]/total_count*100,width=0.75,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
							x += 0.75
							untrans = plt.bar(x, mutations[sample][key][seq][1]/total_count*100,width=0.75,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
							x += .2475
							if mutations[sample][key][seq][0]/total_count*100 > ymax:
									ymax = mutations[sample][key][seq][0]/total_count*100
							if mutations[sample][key][seq][1]/total_count*100 > ymax:
									ymax = mutations[sample][key][seq][1]/total_count*100

					else:
						trans = plt.bar(x, mutations[sample][key][seq][0],width=0.75,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
						x += 0.75
						untrans = plt.bar(x, mutations[sample][key][seq][1],width=0.75,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
						x += .2475
						if mutations[sample][key][seq][0] > ymax:
								ymax = mutations[sample][key][seq][0]
						if mutations[sample][key][seq][1] > ymax:
								ymax = mutations[sample][key][seq][1]
					x += 1
				i += 1

			x = .0415
			y3 = .87
			y = int(ymax*1.25)
			x_plot = 0


			yText = y3 + .06
			plt.text(.1, yText, 'C>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.255, yText, 'C>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.415, yText, 'C>T', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.575, yText, 'T>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.735, yText, 'T>C', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.89, yText, 'T>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

			if y <= 4:
				y += 4

			while y%4 != 0:
				y += 1

			ytick_offest = int(y/4)

			for i in range(0, 6, 1):
				panel1.add_patch(plt.Rectangle((x,y3), .155, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
				panel1.add_patch(plt.Rectangle((x_plot,0), 32, y, facecolor=colors[i], zorder=0, alpha = 0.25, edgecolor='grey'))
				x += .1585
				x_plot += 32

			if percentage:
				ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
				ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%",
						  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
			else:
				ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
				ylabels= [0, ytick_offest, ytick_offest*2,
						  ytick_offest*3, ytick_offest*4]

			labs = np.arange(0.750,192.750,1)



			if not percentage:
				ylabels = ['{:,}'.format(int(x)) for x in ylabels]
				if len(ylabels[-1]) > 3:
					ylabels_temp = []
					if len(ylabels[-1]) > 7:
						for label in ylabels:
							if len(label) > 7:
								ylabels_temp.append(label[0:-8] + "m")
							elif len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)

					else:
						for label in ylabels:
							if len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)
					ylabels = ylabels_temp

			panel1.set_xlim([0, 96])
			panel1.set_ylim([0, y])
			panel1.set_xticks(labs)
			panel1.set_yticks(ylabs)
			count = 0
			m = 0
			for i in range (0, 96, 1):
				plt.text(i/101 + .0415, .02, xlabels[i][0], fontsize=30, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				plt.text(i/101 + .0415, .044, xlabels[i][1], fontsize=30, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
				plt.text(i/101 + .0415, .071, xlabels[i][2], fontsize=30, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				count += 1
				if count == 16:
					count = 0
					m += 1

			if sig_probs:
				plt.text(0.045, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
			else:
				plt.text(0.045, 0.75, sample + ": " + "{:,}".format(int(total_count)) + " transcribed subs", fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)


			custom_text_upper_plot = ''
			try:
				custom_text_upper[sample_count]
			except:
				custom_text_upper = False
			try:
				custom_text_bottom[sample_count]
			except:
				custom_text_bottom = False

			if custom_text_upper:
				plot_custom_text = True
				if len(custom_text_upper[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False
			if custom_text_bottom:
				if len(custom_text_bottom[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False

			if plot_custom_text:
				x_pos_custom = 0.84
				if custom_text_upper and custom_text_bottom:
					custom_text_upper_plot = custom_text_upper[sample_count] + "\n" + custom_text_bottom[sample_count]

				if custom_text_upper and not custom_text_bottom:
					custom_text_upper_plot = custom_text_upper[sample_count]
					panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

				elif custom_text_upper and custom_text_bottom:
					panel1.text(x_pos_custom, 0.72, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

				elif not custom_text_upper and custom_text_bottom:
					custom_text_upper_plot = custom_text_bottom[sample_count]
					panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')




			panel1.set_yticklabels(ylabels, fontsize=30)
			plt.gca().yaxis.grid(True)
			plt.gca().grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
			panel1.set_xlabel('')
			panel1.set_ylabel('')
			plt.legend(handles=[trans, untrans], prop={'size':30})
			if percentage:
				plt.ylabel("Percentage of Single Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')
			else:
				plt.ylabel("Number of Single Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')

			panel1.tick_params(axis='both',which='both',\
							   bottom=False, labelbottom=False,\
							   left=False, labelleft=True,\
							   right=False, labelright=False,\
							   top=False, labeltop=False,\
							   direction='in', length=25, colors=[0.6, 0.6, 0.6])

			[i.set_color("black") for i in plt.gca().get_yticklabels()]

			buffer = io.BytesIO()
			plt.savefig(buffer, format="png")
			buff_list[sample]=buffer
			plt.close()
			sample_count += 1
		return buff_list
	# IMPORTANT: SBS6 has yet to be updated
	elif plot_type == '6':
		with open(matrix_path) as f:
			next(f)
			first_line = f.readline()
			first_line = first_line.strip().split()
			if len(first_line[0]) > 3:
				sys.exit("The matrix does not match the correct SBS6 format. Please check you formatting and rerun this plotting function.")

		pp = PdfPages(output_path + 'SBS_6_plots_' + project + '.pdf')

		mutations = OrderedDict()
		total_count = []
		try:
			with open (matrix_path) as f:
				first_line = f.readline()
				samples = first_line.strip().split("\t")
				samples = samples[1:]
				for sample in samples:
					mutations[sample] = OrderedDict()
					mutations[sample]['C>A'] = 0
					mutations[sample]['C>G'] = 0
					mutations[sample]['C>T'] = 0
					mutations[sample]['T>A'] = 0
					mutations[sample]['T>C'] = 0
					mutations[sample]['T>G'] = 0

				for lines in f:
					line = lines.strip().split()
					nuc = line[0]
					mut_type = line[0]
					sample_index = 1

					for sample in samples:
						if percentage:
							mutCount = float(line[sample_index])
							if mutCount < 1 and mutCount > 0:
								sig_probs = True
						else:
							mutCount = int(line[sample_index])
						mutations[sample][mut_type] = mutCount
						sample_index += 1

			for sample in mutations:
				total_count = sum(mutations[sample].values())
				plt.rcParams['axes.linewidth'] = 2
				plot1 = plt.figure(figsize=(15,10))
				plt.rc('axes', edgecolor='lightgray')
				panel1 = plt.axes([0.12, 0.12, 0.8, 0.77])
				xlabels = []

				y = -0.5
				xmax = 0
				colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
				i = 0
				for key in mutations[sample]:
					xlabels.append(key)
					if percentage:
						if total_count > 0:
							plt.barh(y, mutations[sample][key]/total_count*100,height=0.7,color=colors[i],align='center', zorder=1000)
							if mutations[sample][key]/total_count*100 > xmax:
								xmax = mutations[sample][key]/total_count*100
					else:
						plt.barh(y, mutations[sample][key],height=0.7,color=colors[i],align='center', zorder=1000)
						if mutations[sample][key] > xmax:
								xmax = mutations[sample][key]
					y -= 1
					i += 1

				y = .043
				x3 = .87
				x = int(xmax*1.1)


				while x%4 != 0:
					x += 1
				xtick_offest = int(x/4)

				if percentage:
					xlabs = [0, round(xtick_offest, 1), round(xtick_offest*2, 1), round(xtick_offest*3, 1), round(xtick_offest*4, 1)]
					xlabels= [str(0), str(round(xtick_offest, 1)) + "%", str(round(xtick_offest*2, 1)) + "%",
							  str(round(xtick_offest*3, 1)) + "%", str(round(xtick_offest*4, 1)) + "%"]
				else:
					xlabs = [0, xtick_offest, xtick_offest*2, xtick_offest*3, xtick_offest*4]
					xlabels= [0, xtick_offest, xtick_offest*2,
							  xtick_offest*3, xtick_offest*4]

				# if not percentage:
				# 	xlabels = ['{:,}'.format(int(x)) for x in xlabels]

				if not percentage:
					xlabels = ['{:,}'.format(int(x)) for x in xlabels]
					if len(xlabels[-1]) > 3:
						xlabels_temp = []
						if len(xlabels[-1]) > 7:
							for label in xlabels:
								if len(label) > 7:
									xlabels_temp.append(label[0:-8] + "m")
								elif len(label) > 3:
									xlabels_temp.append(label[0:-4] + "k")
								else:
									xlabels_temp.append(label)

						else:
							for label in xlabels:
								if len(label) > 3:
									xlabels_temp.append(label[0:-4] + "k")
								else:
									xlabels_temp.append(label)
						xlabels = xlabels_temp

				ylabs = np.arange(-5.5, 0.5, 1)
				ylabels = (['T>G','T>C','T>A','C>T','C>G','C>A'])
				panel1.set_xlim([0, x])
				panel1.set_ylim([-6, 0])
				panel1.set_xticks(xlabs)
				panel1.set_yticks(ylabs)
				panel1.set_xticklabels(xlabels, fontsize=30)
				panel1.set_yticklabels(ylabels, fontsize=30)
				panel1.spines['right'].set_visible(False)
				panel1.spines['top'].set_visible(False)

				if sig_probs:
					plt.text(.125, .9, sample, fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
				else:
					plt.text(.125, .9, sample + ": " + "{:,}".format(int(total_count)) + " subs", fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)


				if percentage:
					plt.xlabel("Percentage of Single Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')
				else:
					plt.xlabel("Number of Single Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')

				panel1.set_ylabel('')

				panel1.tick_params(axis='both',which='both',\
								   bottom=True, labelbottom=True,\
								   left=False, labelleft=True,\
								   right=False, labelright=False,\
								   top=False, labeltop=False,\
								   width=2)



				pp.savefig(plot1)
				plt.close()
			pp.close()

		except:
			print("There may be an issue with the formatting of you matrix file.")
			os.remove(output_path + 'SBS_6_plots_' + project + '.pdf')
	# IMPORTANT: SBS24 has yet to be updated
	elif plot_type == '12' or plot_type == '6SB' or plot_type == '24':
		with open(matrix_path) as f:
			next(f)
			first_line = f.readline()
			first_line = first_line.strip().split()
			if first_line[0][1] != ":" or len(first_line[0]) != 5:
				sys.exit("The matrix does not match the correct SBS192 format. Please check you formatting and rerun this plotting function.")

		pp = PdfPages(output_path + 'SBS_24_plots_' + project + '.pdf')

		mutations = OrderedDict()

		try:
			with open (matrix_path) as f:
				first_line = f.readline()
				samples = first_line.strip().split("\t")
				samples = samples[1:]
				for sample in samples:
					mutations[sample] = OrderedDict()
					mutations[sample]['C>A'] = [0,0]
					mutations[sample]['C>G'] = [0,0]
					mutations[sample]['C>T'] = [0,0]
					mutations[sample]['T>A'] = [0,0]
					mutations[sample]['T>C'] = [0,0]
					mutations[sample]['T>G'] = [0,0]

				for lines in f:
					line = lines.strip().split()
					nuc = line[0][2:]
					bias = line[0][0]
					if bias == 'N' or bias == 'B':
						continue
					else:
						sample_index = 1
						for sample in samples:
							if percentage:
								mutCount = float(line[sample_index])
								if mutCount < 1 and mutCount > 0:
									sig_probs = True
							else:
								mutCount = int(line[sample_index])
							if bias == 'T':
								mutations[sample][nuc][0] = mutCount
							else:
								mutations[sample][nuc][1] = mutCount
							sample_index += 1
			for sample in mutations:
				total_count = sum(sum(tsb) for tsb in mutations[sample].values())
				plt.rcParams['axes.linewidth'] = 2
				plot1 = plt.figure(figsize=(15,10))
				plt.rc('axes', edgecolor='lightgray')
				panel1 = plt.axes([0.12, 0.12, 0.8, 0.77])

				y = 12.485
				xmax = 0
				for key in mutations[sample]:
					if percentage:
						if total_count > 0:
							trans = plt.barh(y, mutations[sample][key][0]/total_count*100,height=0.75,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
							y -= 0.75
							untrans = plt.barh(y, mutations[sample][key][1]/total_count*100,height=0.75,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
							y -= .2475
							if mutations[sample][key][0]/total_count*100 > xmax:
									xmax = mutations[sample][key][0]/total_count*100
							if mutations[sample][key][1]/total_count*100 > xmax:
									xmax = mutations[sample][key][1]/total_count*100

					else:
						trans = plt.barh(y, mutations[sample][key][0],height=0.75,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
						y -= 0.75
						untrans = plt.barh(y, mutations[sample][key][1],height=0.75,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
						y -= .2475
						if mutations[sample][key][0] > xmax:
								xmax = mutations[sample][key][0]
						if mutations[sample][key][1] > xmax:
								xmax = mutations[sample][key][1]
					y -= 1

				x = int(xmax*1.1)

				while x%4 != 0:
					x += 1

				xtick_offest = int(x/4)


				if percentage:
					xlabs = [0, round(xtick_offest, 1), round(xtick_offest*2, 1), round(xtick_offest*3, 1), round(xtick_offest*4, 1)]
					xlabels= [str(0), str(round(xtick_offest, 1)) + "%", str(round(xtick_offest*2, 1)) + "%",
							  str(round(xtick_offest*3, 1)) + "%", str(round(xtick_offest*4, 1)) + "%"]
				else:
					xlabs = [0, xtick_offest, xtick_offest*2, xtick_offest*3, xtick_offest*4]
					xlabels= [0, xtick_offest, xtick_offest*2,
							  xtick_offest*3, xtick_offest*4]
				if not percentage:
					xlabels = ['{:,}'.format(int(x)) for x in xlabels]
					if len(xlabels[-1]) > 3:
						xlabels_temp = []
						if len(xlabels[-1]) > 7:
							for label in xlabels:
								if len(label) > 7:
									xlabels_temp.append(label[0:-8] + "m")
								elif len(label) > 3:
									xlabels_temp.append(label[0:-4] + "k")
								else:
									xlabels_temp.append(label)

						else:
							for label in xlabels:
								if len(label) > 3:
									xlabels_temp.append(label[0:-4] + "k")
								else:
									xlabels_temp.append(label)
						xlabels = xlabels_temp
				# if not percentage:
				# 	xlabels = ['{:,}'.format(int(x)) for x in xlabels]
				ylabs = np.arange(2.15, 13, 2)
				ylabels = (['T>G','T>C','T>A','C>T','C>G','C>A'])
				#labs = np.arange(0,12.485,)
				panel1.set_xlim([0, x])
				panel1.set_ylim([1.2524, 13.235])
				panel1.set_yticks(ylabs)
				panel1.set_xticks(xlabs)
				panel1.set_xticklabels(xlabels, fontsize=30)
				panel1.set_yticklabels(ylabels, fontsize=30)
				panel1.set_xlabel('')
				panel1.set_ylabel('')
				panel1.spines['right'].set_visible(False)
				panel1.spines['top'].set_visible(False)

				if sig_probs:
					plt.text(.125, .9, sample, fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
				else:
					plt.text(.125, .9, sample + ": " + "{:,}".format(int(total_count)) + " transcribed subs", fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

				if percentage:
					plt.xlabel("Percentage of Single Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')
				else:
					plt.xlabel("Number of Single Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')

				panel1.tick_params(axis='both',which='both',\
								   bottom=True, labelbottom=True,\
								   left=False, labelleft=True,\
								   right=False, labelright=False,\
								   top=False, labeltop=False,\
								   width=2)

				plt.legend(handles=[trans, untrans], prop={'size':25})

				pp.savefig(plot1)
				plt.close()
			pp.close()

		except:
			print("There may be an issue with the formatting of you matrix file.")
			os.remove(output_path + 'SBS_24_plots_' + project + '.pdf')

	elif plot_type == '288':
		first_line=matrix_path.iloc[0,:]
		if first_line[0][6] == ">" or first_line[0][3] == ">":
			pcawg = True
		if first_line[0][7] != "]" and first_line[0][6] != ">" and first_line[0][3] != ">":
			sys.exit("The matrix does not match the correct SBS288 format. Please check you formatting and rerun this plotting function.")

		#pp = PdfPages(output_path + 'SBS_288_plots_' + project + '.pdf')

		mutations = OrderedDict()
		mutations_TSB = OrderedDict()
		total_count = []

		if pcawg:
			samples=matrix_path.columns[:]
			samples = samples[2:]
		else:
			samples=matrix_path.columns[:]
			samples = samples[1:]
			
		for sample in samples:
			mutations[sample] = OrderedDict()
			mutations[sample]['C>A'] = OrderedDict()
			mutations[sample]['C>G'] = OrderedDict()
			mutations[sample]['C>T'] = OrderedDict()
			mutations[sample]['T>A'] = OrderedDict()
			mutations[sample]['T>C'] = OrderedDict()
			mutations[sample]['T>G'] = OrderedDict()

			mutations_TSB[sample] = OrderedDict()
			mutations_TSB[sample]['All'] = OrderedDict({'T':0, 'U':0,'N':0})
			mutations_TSB[sample]['C>A'] = OrderedDict({'T':0, 'U':0,'N':0})
			mutations_TSB[sample]['C>G'] = OrderedDict({'T':0, 'U':0,'N':0})
			mutations_TSB[sample]['C>T'] = OrderedDict({'T':0, 'U':0,'N':0})
			mutations_TSB[sample]['T>A'] = OrderedDict({'T':0, 'U':0,'N':0})
			mutations_TSB[sample]['T>C'] = OrderedDict({'T':0, 'U':0,'N':0})
			mutations_TSB[sample]['T>G'] = OrderedDict({'T':0, 'U':0,'N':0})

		for lines_tmp in range(0,matrix_path.shape[0]):
			if pcawg:
				line = matrix_path.iloc[lines_tmp,:]
				mut_type = line[0]
				nuc = line[1][0] + "[" + mut_type + "]" + line[1][2]
				sample_index = 2
			else:
				line = matrix_path.iloc[lines_tmp,:]
				nuc = line[0]
				mut_type = line[0][4:7]
				sample_index = 1
				tsb = nuc[0]


			for sample in samples:
				if percentage:
					mutCount = float(line[sample_index])
					if mutCount < 1 and mutCount > 0:
						sig_probs = True
				else:
					try:
						mutCount = int(line[sample_index])
					except:
						print("It appears that the provided matrix does not contain mutation counts.\n\tIf you have provided a signature activity matrix, please change the percentage parameter to True.\n\tOtherwise, ", end='')
				if nuc[2:] not in mutations[sample][mut_type]:
					mutations[sample][mut_type][nuc[2:]] = mutCount
				else:
					mutations[sample][mut_type][nuc[2:]] += mutCount
				mutations_TSB[sample][mut_type][tsb] += mutCount
				mutations_TSB[sample]['All'][tsb] += mutCount
				sample_index += 1

		sample_count = 0

		for sample in mutations.keys():
			total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
			plt.rcParams['axes.linewidth'] = 2
			plot1 = plt.figure(figsize=(43.93,9.92))
			plt.rc('axes', edgecolor='lightgray')
			panel1 = plt.axes([0.04, 0.09, 0.7, 0.77])
			panel2 = plt.axes([0.77, 0.09, 0.21, 0.77])
			xlabels = []
			
			x = 0.4
			ymax = 0
			colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
			i = 0
			for key in mutations[sample]:
				for seq in mutations[sample][key]:
					xlabels.append(seq[0]+seq[2]+seq[6])
					if percentage:
						if total_count > 0: 
							panel1.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
							if mutations[sample][key][seq]/total_count*100 > ymax:
								ymax = mutations[sample][key][seq]/total_count*100
					else:
						panel1.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
						if mutations[sample][key][seq] > ymax:
								ymax = mutations[sample][key][seq]
					x += 1
				i += 1

			x = .043
			y3 = .87
			y = int(ymax*1.25)
			y2 = y+2
			for i in range(0, 6, 1):
				panel1.add_patch(plt.Rectangle((x,y3), .11, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
				x += .117

			yText = y3 + .06
			

			plt.text(.082, yText, 'C>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.1975, yText, 'C>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.315, yText, 'C>T', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.43, yText, 'T>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.55, yText, 'T>C', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.665, yText, 'T>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

			if y <= 4:
				y += 4

			while y%4 != 0:
				y += 1
			y = ymax/1.025
			ytick_offest = float(y/3)
			# ytick_offest = int(y/4)


			if percentage:
				ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
				ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
						  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
				# ylabels= [str(0), str(round(ytick_offest)) + "%", str(round(ytick_offest*2)) + "%", 
				# 		  str(round(ytick_offest*3)) + "%", str(round(ytick_offest*4)) + "%"]
			else:
				ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
				ylabels= [0, ytick_offest, ytick_offest*2, 
						  ytick_offest*3, ytick_offest*4]       

			font_label_size = 30
			if not percentage:
				if int(ylabels[3]) >= 1000:
					font_label_size = 20

			if percentage:
				if len(ylabels) > 2:
					font_label_size = 20
			
			labs = np.arange(0.375,96.375,1)

			if not percentage:
				ylabels = ['{:,}'.format(int(x)) for x in ylabels]
				if len(ylabels[-1]) > 3:
					ylabels_temp = []
					if len(ylabels[-1]) > 7:
						for label in ylabels:
							if len(label) > 7:
								ylabels_temp.append(label[0:-8] + "m")
							elif len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)

					else:
						for label in ylabels:
							if len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)
					ylabels = ylabels_temp

			panel1.set_xlim([0, 96])
			panel1.set_ylim([0, y])
			panel1.set_xticks(labs)
			panel1.set_yticks(ylabs)
			count = 0
			m = 0
			for i in range (0, 96, 1):
				plt.text(i/137 + .04, .02, xlabels[i][0], fontsize=25, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				plt.text(i/137 + .04, .044, xlabels[i][1], fontsize=25, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
				plt.text(i/137 + .04, .071, xlabels[i][2], fontsize=25, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				count += 1
				if count == 16:
					count = 0
					m += 1  

			if sig_probs:
				plt.text(0.045, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
			else:
				plt.text(0.045, 0.75, sample + ": " + "{:,}".format(int(total_count)) + " subs", fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)


			
			custom_text_upper_plot = ''
			try:
				custom_text_upper[sample_count]
			except:
				custom_text_upper = False
			try:
				custom_text_middle[sample_count]
			except:
				custom_text_middle = False
			try:
				custom_text_bottom[sample_count]
			except:
				custom_text_bottom = False
									
			if custom_text_upper: 
				plot_custom_text = True
				if len(custom_text_upper[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False
			if custom_text_middle:
				if len(custom_text_middle[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False

			if plot_custom_text:
				x_pos_custom = 0.73
				if custom_text_upper and custom_text_middle:
					custom_text_upper_plot = custom_text_upper[sample_count] + "\n" + custom_text_middle[sample_count]
					if custom_text_bottom:
						custom_text_upper_plot += "\n" + custom_text_bottom[sample_count]

				if custom_text_upper and not custom_text_middle:
					custom_text_upper_plot = custom_text_upper[sample_count]
					panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right',zorder=1)
					
				elif custom_text_upper and custom_text_middle:
					if not custom_text_bottom:
						panel1.text(x_pos_custom, 0.72, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')
					else:
						panel1.text(x_pos_custom, 0.68, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

				elif not custom_text_upper and custom_text_middle:
					custom_text_upper_plot = custom_text_middle[sample_count]
					panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')
			


			panel1.set_yticklabels(ylabels, fontsize=font_label_size)
			panel1.yaxis.grid(True)
			panel1.grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
			panel1.set_xlabel('')
			panel1.set_ylabel('')

			if percentage:
				panel1.set_ylabel("Percentage of Single Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')
			else:
				panel1.set_ylabel("Number of Single Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')



			panel1.tick_params(axis='both',which='both',\
							   bottom=False, labelbottom=False,\
							   left=True, labelleft=True,\
							   right=True, labelright=False,\
							   top=False, labeltop=False,\
							   direction='in', length=25, colors='lightgray', width=2)


			[i.set_color("black") for i in panel1.get_yticklabels()]


			yp2 = 28
			labels = []
			y2max = 0
			tsbColors = [[1/256,70/256,102/256], [228/256,41/256,38/256], 'green']
			for mut in mutations_TSB[sample]:
				labels.append(mut)
				i = 0
				for tsb in mutations_TSB[sample][mut]:
					if tsb == "T":
						label = "Transcribed"
					elif tsb == "U":
						label = "Untranscribed"
					else:
						label = "Nontranscribed"
					if percentage:
						if total_count > 0: 
							panel2.barh(yp2, mutations_TSB[sample][mut][tsb]/total_count*100,color=tsbColors[i], label=label)
							if mutations_TSB[sample][mut][tsb]/total_count*100 > y2max:
								y2max = mutations_TSB[sample][mut][tsb]/total_count*100
					else:
						if mutations_TSB[sample][mut][tsb] > y2max:
							y2max = mutations_TSB[sample][mut][tsb]

						panel2.barh(yp2, mutations_TSB[sample][mut][tsb], color=tsbColors[i], label=label)
					yp2 -= 1
					i += 1
				yp2 -=1 

			y = int(y2max*1.1)
			if y <= 4:
				y += 4

			while y%4 != 0:
				y += 1
			ytick_offest = int(y/4)

			if percentage:
				xlabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
				xlabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
						  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
			else:
				xlabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
				xlabels= [0, ytick_offest, ytick_offest*2, 
						  ytick_offest*3, ytick_offest*4]


			if not percentage:
				xlabels = ['{:,}'.format(int(x)) for x in xlabels]
				if len(xlabels[-1]) > 3:
					xlabels_temp = []
					if len(xlabels[-1]) > 7:
						for label in xlabels:
							if len(label) > 7:
								xlabels_temp.append(label[0:-8] + "m")
							elif len(label) > 3:
								xlabels_temp.append(label[0:-4] + "k")
							else:
								xlabels_temp.append(label)

					else:
						for label in xlabels:
							if len(label) > 3:
								xlabels_temp.append(label[0:-4] + "k")
							else:
								xlabels_temp.append(label)
					xlabels = xlabels_temp
			panel2.spines['right'].set_visible(False)
			panel2.spines['top'].set_visible(False)
			labels.reverse()
			panel2.set_yticks([3, 7, 11, 15, 19, 23, 27])
			panel2.set_yticklabels(labels, fontsize=30,fontname="Arial", weight = 'bold')
			panel2.set_xticklabels(xlabels, fontsize=30)
			panel2.set_xticks(xlabs)
			handles, labels = panel2.get_legend_handles_labels()
			panel2.legend(handles[:3], labels[:3], loc='best', prop={'size':30})
			buffer = io.BytesIO()
			plt.savefig(buffer, format="png")
			buff_list[sample]=buffer
			plt.close()
			sample_count += 1
		return buff_list

###########################################################################################################################
	elif plot_type == '1536':
		first_line=matrix_path.iloc[0,:]
		if first_line[0][1] == '>':
			pcawg = True
		if first_line[0][6] != "]" and first_line[0][1] != '>':
			sys.exit("The matrix does not match the correct SBS1536 format. Please check you formatting and rerun this plotting function.")

		mutations_96 = OrderedDict()
		mutations = OrderedDict()
		mutations_5 = OrderedDict()
		mutations_3 = OrderedDict()
		max_count = {}
		max_all = {}
		max_5 = {}
		max_3 = {}
		total_count = []
		total_counts = {'TT':0, 'TG':0,'TC':0,'TA':0,
						'GT':0,'GG':0,'GC':0,'GA':0,
						'CT':0,'CG':0,'CC':0,'CA':0,
						'AT':0,'AG':0,'AC':0,'AA':0,}
		total_counts_5 = {'T':0, 'G':0,'C':0,'A':0}
		total_counts_3 = {'T':0, 'G':0,'C':0,'A':0}

		first_line=matrix_path.iloc[0,:]
		if pcawg:
			samples=matrix_path.columns[:]
			samples = samples[2:]
			samples = [x.replace('"','') for x in samples]
		else:
			samples=matrix_path.columns[:]
			samples = samples[1:]

		for sample in samples:
			max_all[sample] = 0
			max_5[sample] = 0
			max_3[sample] = 0
			total_counts[sample]= {'TT':0, 'TG':0,'TC':0,'TA':0,
								   'GT':0,'GG':0,'GC':0,'GA':0,
								   'CT':0,'CG':0,'CC':0,'CA':0,
								   'AT':0,'AG':0,'AC':0,'AA':0,}
			total_counts_5[sample]= {'T':0, 'G':0,'C':0,'A':0}
			total_counts_3[sample]= {'T':0, 'G':0,'C':0,'A':0}

			mutations_96[sample] = OrderedDict()
			mutations_96[sample]['C>A'] = OrderedDict()
			mutations_96[sample]['C>G'] = OrderedDict()
			mutations_96[sample]['C>T'] = OrderedDict()
			mutations_96[sample]['T>A'] = OrderedDict()
			mutations_96[sample]['T>C'] = OrderedDict()
			mutations_96[sample]['T>G'] = OrderedDict()

			max_count[sample] = 0
			mutations[sample] = OrderedDict()
			mutations_5[sample] = OrderedDict()
			mutations_3[sample] = OrderedDict()

			mutations[sample]['C>A'] = OrderedDict()
			mutations_5[sample]['C>A'] = OrderedDict()
			mutations_3[sample]['C>A'] = OrderedDict()
			mutations[sample]['C>A'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
										'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
										'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(),
										'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
			mutations_5[sample]['C>A'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
			mutations_3[sample]['C>A'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}

			mutations[sample]['C>G'] = OrderedDict()
			mutations_5[sample]['C>G'] = OrderedDict()
			mutations_3[sample]['C>G'] = OrderedDict()
			mutations[sample]['C>G'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
										'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
										'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(),
										'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
			mutations_5[sample]['C>G'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
			mutations_3[sample]['C>G'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}

			mutations[sample]['C>T'] = OrderedDict()
			mutations_5[sample]['C>T'] = OrderedDict()
			mutations_3[sample]['C>T'] = OrderedDict()
			mutations[sample]['C>T'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
										'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
										'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(),
										'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
			mutations_5[sample]['C>T'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
			mutations_3[sample]['C>T'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}

			mutations[sample]['T>A'] = OrderedDict()
			mutations_5[sample]['T>A'] = OrderedDict()
			mutations_3[sample]['T>A'] = OrderedDict()
			mutations[sample]['T>A'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
										'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
										'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(),
										'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
			mutations_5[sample]['T>A'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
			mutations_3[sample]['T>A'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}

			mutations[sample]['T>C'] = OrderedDict()
			mutations_5[sample]['T>C'] = OrderedDict()
			mutations_3[sample]['T>C'] = OrderedDict()
			mutations[sample]['T>C'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
										'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
										'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(),
										'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
			mutations_5[sample]['T>C'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
			mutations_3[sample]['T>C'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}

			mutations[sample]['T>G'] = OrderedDict()
			mutations_5[sample]['T>G'] = OrderedDict()
			mutations_3[sample]['T>G'] = OrderedDict()
			mutations[sample]['T>G'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
										'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
										'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(),
										'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
			mutations_5[sample]['T>G'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
			mutations_3[sample]['T>G'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}


		for lines_tmp in range(0,matrix_path.shape[0]):
			if pcawg:
				line = matrix_path.iloc[lines_tmp,:]
				line = [x.replace('"','') for x in line]
				nuc = line[1][0:2] + "[" + line[0] + "]" + line[1][3:]
				mut_type = line[0]
				penta_key = line[1][0] + line[1][-1]
				tri_key = line[1][1] + line[1][-2]
				sample_index = 2
			else:
				line = matrix_path.iloc[lines_tmp,:]
				nuc = line[0]
				mut_type = line[0][3:6]
				penta_key = line[0][0] + line[0][-1]
				tri_key = line[0][1] + line[0][-2]
				sample_index = 1

				tri = line[0][1:8]

			for sample in samples:

				if tri not in mutations_96[sample][mut_type]:
					mutations_96[sample][mut_type][tri] = 0
				if percentage:
					mutCount = float(line[sample_index])
					if mutCount < 1 and mutCount > 0:
						sig_probs = True
				else:
					mutCount = int(line[sample_index])

				if pcawg:
					sample_ref = sample_index - 2
				else:
					sample_ref = sample_index - 1
				if mutCount > max_count[samples[sample_ref]]:
					max_count[samples[sample_ref]] = mutCount

				if mutCount > max_all[sample]:
					max_all[sample] = mutCount

				mutations[sample][mut_type][penta_key][tri_key] = mutCount
				total_counts[sample][penta_key] += mutCount
				total_counts_5[sample][penta_key[0]] += mutCount
				total_counts_3[sample][penta_key[1]] += mutCount
				penta_key_short = penta_key[0]
				mutations_5[sample][mut_type][penta_key_short][tri_key] = 0
				mutations_3[sample][mut_type][penta_key_short][tri_key] = 0
				mutations_96[sample][mut_type][tri] += mutCount
				sample_index += 1


		sample_count = 0
		for sample in mutations.keys():
			total_count_sample = sum(sum(nuc.values()) for nuc in mutations_96[sample].values())
			total_count = max_all[sample]*1.1
			ratio = total_count/total_count_sample
			plt.rcParams['axes.linewidth'] = 2
			plot1 = plt.figure(figsize=(43.93,22.5))
			plt.rc('axes', edgecolor='lightgray')
			panel1 = plt.axes([0.03, 0.0677, 0.92, 0.267])
			panel2 = plt.axes([0.03, 0.65, 0.92, 0.267])
			panel3 = plt.axes([0.03, 0.35, 0.92, 0.1335]) # 3' context
			panel4 = plt.axes([0.03, 0.5, 0.92, 0.1335]) # 5' context
			xlabels = []
			ylabels = []
			ylabels_5 = []
			ylabels_3 = []


			colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
			colors_heat = [np.linspace(56/255,255/255,5), np.linspace(66/255,225/255,5), np.linspace(157/255,40/255,5)]
			colors_heat_compact = [np.linspace(56/255,255/255,5), np.linspace(66/255,225/255,5), np.linspace(157/255,40/255,5)]


			i = 0
			x_pos = 0
			x_inter = 0
			for key in mutations[sample]:
				y_pos = 15
				for penta in mutations[sample][key]:
					#total_count = total_counts[sample][penta]
					key_5 = penta[0]
					key_3 = penta[1]
					ylabels.append(penta[0] + "---" + penta[1])
					for tri in mutations[sample][key][penta]:
						tri_nuc = tri[0] + "[" + key + "]" + tri[1]
						normalized = mutations_96[sample][key][tri_nuc]
						try:
							mut_count = int(int(20 * round(float(mutations[sample][key][penta][tri]/total_count_sample/ratio* 100))/20)/20)
							mutations_5[sample][key][key_5][tri] += float(mutations[sample][key][penta][tri])
							mutations_3[sample][key][key_3][tri] += float(mutations[sample][key][penta][tri])
							if mutations_5[sample][key][key_5][tri] > max_5[sample]:
								max_5[sample] = mutations_5[sample][key][key_5][tri]
							if mutations_3[sample][key][key_3][tri] > max_3[sample]:
								max_3[sample] = mutations_3[sample][key][key_3][tri]

						except:
							mut_count = 0
						xlabels.append(tri[0]+"-"+tri[1])
						rectangle=mplpatches.Rectangle((x_pos, y_pos), 1, 1,\
														linewidth=1,\
														facecolor=(colors_heat[0][mut_count], colors_heat[1][mut_count], colors_heat[2][mut_count]))

						panel1.add_patch(rectangle)
						x_pos += 1
					y_pos -= 1
					x_pos = x_inter

				x_inter += 17
				x_pos = x_inter
				i += 1



			x_pos = 0
			x_inter = 0
			total_count_5 = max_5[sample]*1.1
			total_count_3 = max_3[sample]*1.1
			ratio_5 = total_count_5/total_count_sample
			ratio_3 = total_count_3/total_count_sample
			ratio_total = max(ratio_5, ratio_3)

			for key in mutations_5[sample]:
				y_pos = 3
				for penta in mutations_5[sample][key]:
					# total_count_5 = total_counts_5[sample][penta]
					# total_count_3 = total_counts_3[sample][penta]

					ylabels_5.append(penta + "---N" )
					ylabels_3.append("N---" + penta)
					for tri in mutations_5[sample][key][penta]:
						mut_count = int(int(20 * round(float(mutations_5[sample][key][penta][tri])/total_count_sample/ratio_total*100)/20)/20)
						mut_count_3 = int(int(20 * round(float(mutations_3[sample][key][penta][tri])/total_count_sample/ratio_total*100)/20)/20)
						rectangle=mplpatches.Rectangle((x_pos, y_pos), 1, 1,\
														linewidth=1,\
														facecolor=(colors_heat_compact[0][mut_count], colors_heat_compact[1][mut_count], colors_heat_compact[2][mut_count]))

						panel4.add_patch(rectangle)
						rectangle=mplpatches.Rectangle((x_pos, y_pos), 1, 1,\
														linewidth=1,\
														facecolor=(colors_heat[0][mut_count_3], colors_heat[1][mut_count_3], colors_heat[2][mut_count_3]))

						panel3.add_patch(rectangle)
						x_pos += 1
					y_pos -= 1
					x_pos = x_inter
				x_inter += 17
				x_pos = x_inter
				i += 1

			x = 0.5
			ymax = 0
			colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
			i = 0

			for key in mutations_96[sample]:
				for seq in mutations_96[sample][key]:
					xlabels.append(seq[0]+seq[2]+seq[6])
					if percentage:
						if total_count_sample > 0:
							panel2.bar(x, mutations_96[sample][key][seq]/total_count_sample*100,width=0.5,color=colors[i],align='center', zorder=1000)
							if mutations_96[sample][key][seq]/total_count_sample*100 > ymax:
								ymax = mutations_96[sample][key][seq]/total_count_sample*100
					else:
						panel2.bar(x, mutations_96[sample][key][seq],width=0.5,color=colors[i],align='center', zorder=1000)
						if mutations_96[sample][key][seq] > ymax:
								ymax = mutations_96[sample][key][seq]
					x += 1
				x += 1
				i += 1


			# scale bar for bottom 1536 plot
			y_grad = .267/len(colors_heat[0])
			y_start = 0.0677
			for l in range (0, len(colors_heat[0]), 1):
				rectangle=mplpatches.Rectangle((.96, y_start), .02, y_grad,\
												linewidth=1,\
												facecolor=(colors_heat[0][l], colors_heat[1][l], colors_heat[2][l]),\
												transform=plt.gcf().transFigure, clip_on=False,\
												edgecolor=(colors_heat[0][l], colors_heat[1][l], colors_heat[2][l]))
				panel1.add_patch(rectangle)
				y_start += y_grad

			# scale bar for top 1536 plot
			y_grad = .1335/len(colors_heat_compact[0])
			y_start = 0.5
			for l in range (0, len(colors_heat_compact[0]), 1):
				rectangle=mplpatches.Rectangle((.96, y_start), .02, y_grad,\
												linewidth=1,\
												facecolor=(colors_heat_compact[0][l], colors_heat_compact[1][l], colors_heat_compact[2][l]),\
												transform=plt.gcf().transFigure, clip_on=False,\
												edgecolor=(colors_heat_compact[0][l], colors_heat_compact[1][l], colors_heat_compact[2][l]))
				panel1.add_patch(rectangle)
				y_start += y_grad

			# scale bar for middle 1536 plot
			y_grad = .1335/len(colors_heat_compact[0])
			y_start = 0.35
			for l in range (0, len(colors_heat_compact[0]), 1):
				rectangle=mplpatches.Rectangle((.96, y_start), .02, y_grad,\
												linewidth=1,\
												facecolor=(colors_heat_compact[0][l], colors_heat_compact[1][l], colors_heat_compact[2][l]),\
												transform=plt.gcf().transFigure, clip_on=False,\
												edgecolor=(colors_heat_compact[0][l], colors_heat_compact[1][l], colors_heat_compact[2][l]))
				panel1.add_patch(rectangle)
				y_start += y_grad
			y_tick_grad = max_count[sample]/2

			# scale numbers for bottom 1536 plot
			plt.text(.9825, .0677, '0', fontsize=15, fontweight='bold', transform=plt.gcf().transFigure)
			plt.text(.9825, .2012, str(ratio/2)[:5], fontsize=15, fontweight='bold', transform=plt.gcf().transFigure)
			plt.text(.9825, .325, str(ratio)[:5], fontsize=15, fontweight='bold', transform=plt.gcf().transFigure)
			y = int(ymax*1.25)

			# scale numbers for top 1536 plot
			plt.text(.9825, .5, '0', fontsize=20, fontweight='bold', transform=plt.gcf().transFigure)
			plt.text(.9825, .56675, str(ratio_total/2)[:5], fontsize=15, fontweight='bold', transform=plt.gcf().transFigure)
			plt.text(.9825, .625, str(ratio_total)[:5], fontsize=15, fontweight='bold', transform=plt.gcf().transFigure)
			y = int(ymax*1.25)

			# scale numbers for middle 1536 plot
			plt.text(.9825, .35, '0', fontsize=20, fontweight='bold', transform=plt.gcf().transFigure)
			plt.text(.9825, .41675, str(ratio_total/2)[:5], fontsize=15, fontweight='bold', transform=plt.gcf().transFigure)
			plt.text(.9825, .475, str(ratio_total)[:5], fontsize=15, fontweight='bold', transform=plt.gcf().transFigure)
			y = int(ymax*1.25)

			x = .033
			y3 = .92
			for i in range(0, 6, 1):
				panel1.add_patch(plt.Rectangle((x,y3), .143, .03, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
				x += .154


			y3 = .9
			yText = y3 + 0.06
			plt.text(.085, yText, 'C>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.24, yText, 'C>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.395, yText, 'C>T', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.552, yText, 'T>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.705, yText, 'T>C', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.86, yText, 'T>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)



			if y <= 4:
				y += 4

			while y%4 != 0:
				y += 1
			ytick_offest = int(y/4)

			ylabels_96 = []
			if percentage:
				ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
				ylabels_96= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%",
						  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
			else:
				ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
				ylabels_96= [0, ytick_offest, ytick_offest*2,
						  ytick_offest*3, ytick_offest*4]

			font_label_size = 25
			if not percentage:
				if int(ylabels_96[3]) >= 1000:
					font_label_size = 20


			if not percentage:
				ylabels_96 = ['{:,}'.format(int(x)) for x in ylabels_96]
				if len(ylabels_96[-1]) > 3:
					ylabels_temp = []
					if len(ylabels_96[-1]) > 7:
						for label in ylabels_96:
							if len(label) > 7:
								ylabels_temp.append(label[0:-8] + "m")
							elif len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)

					else:
						for label in ylabels_96:
							if len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)
					ylabels_96 = ylabels_temp


			panel1.set_xlim([0, 101])
			panel1.set_ylim([0, 16])
			panel2.set_xlim([0, 101])
			panel2.set_yticks(ylabs)
			panel1.set_yticks([])
			panel1.set_xticks([])
			panel2.set_xticks([])
			panel4.set_xlim([0, 101])
			panel4.set_ylim([0, 4])
			panel3.set_xlim([0, 101])
			panel3.set_ylim([0, 4])
			panel4.set_yticks([])
			panel3.set_yticks([])
			panel3.set_xticks([])
			panel4.set_xticks([])



			# x-axis 1536 bottom plot
			m = 0
			count = 0
			x_letter = 0
			for i in range (0, 96, 1):
				plt.text(x_letter/101 + .032, .04, xlabels[i][0], fontsize=25, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				plt.text(x_letter/101 + .032, .05, xlabels[i][1], fontsize=25, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
				plt.text(x_letter/101 + .032, .06, xlabels[i][2], fontsize=25, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				count += 1
				x_letter += .92
				if (i+1)%16 == 0 and i != 0:
					x_letter += .92
				if count == 16:
					count = 0
					m += 1

			# y-axis 1536 bottom plot
			m = 0
			count = 0
			y_letter = 5.2
			for i in range (0, 16, 1):
				plt.text(.003, y_letter/16 , ylabels[i][0], fontsize=25, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				plt.text(.008, y_letter/16 + 0, ylabels[i][1:4], fontsize=25, color='black', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
				plt.text(.022, y_letter/16 + 0, ylabels[i][4], fontsize=25, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				count += 1
				y_letter -= .2675

				if count == 16:
					count = 0
					m += 1

			# y-axis 1536 top matrix plot
			m = 0
			count = 0
			y_letter = 9.85
			for i in range (0, 4, 1):
				plt.text(.003, y_letter/16 , ylabels_5[i][0], fontsize=25, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				plt.text(.008, y_letter/16 + 0, ylabels_5[i][1:4], fontsize=25, color='black', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
				plt.text(.022, y_letter/16 + 0, ylabels_5[i][4], fontsize=25, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				count += 1
				y_letter -= .53

				if count == 16:
					count = 0
					m += 1

			# y-axis 1536 middle matrix plot
			m = 0
			count = 0
			y_letter = 7.45
			for i in range (0, 4, 1):
				plt.text(.003, y_letter/16 , ylabels_3[i][0], fontsize=25, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				plt.text(.008, y_letter/16 + 0, ylabels_3[i][1:4], fontsize=25, color='black', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
				plt.text(.022, y_letter/16 + 0, ylabels_3[i][4], fontsize=25, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
				count += 1
				y_letter -= .53

				if count == 16:
					count = 0
					m += 1


			custom_text_upper_plot = ''
			try:
				custom_text_upper[sample_count]
			except:
				custom_text_upper = False
			try:
				custom_text_middle[sample_count]
			except:
				custom_text_middle = False
			try:
				custom_text_bottom[sample_count]
			except:
				custom_text_bottom = False

			if custom_text_upper:
				plot_custom_text = True
				if len(custom_text_upper[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False
			if custom_text_middle:
				if len(custom_text_middle[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False

			if plot_custom_text:
				x_pos_custom = 0.94
				if custom_text_upper and custom_text_middle:
					custom_text_upper_plot = custom_text_upper[sample_count] + "\n" + custom_text_middle[sample_count]
					if custom_text_bottom:
						custom_text_upper_plot += "\n" + custom_text_bottom[sample_count]

				if custom_text_upper and not custom_text_middle:
					custom_text_upper_plot = custom_text_upper[sample_count]
					panel2.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

				elif custom_text_upper and custom_text_middle:
					if not custom_text_bottom:
						panel2.text(x_pos_custom, 0.86, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')
					else:
						panel2.text(x_pos_custom, 0.835, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

				elif not custom_text_upper and custom_text_middle:
					custom_text_upper_plot = custom_text_middle[sample_count]
					panel2.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')


			if sig_probs:
				panel2.text(0.04, 0.875, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
			else:
				panel2.text(0.04, 0.875, sample + ": " + "{:,}".format(int(total_count_sample)) + " subs", fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)


			panel2.grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
			panel2.set_yticklabels(ylabels_96, fontsize=font_label_size)
			panel1.set_xlabel('')
			panel1.set_ylabel('')
			panel1.set_yticklabels([])
			panel1.set_xticklabels([])
			panel2.set_xticklabels([])

			if percentage:
				panel2.set_ylabel("Percentage of Single Base Substitutions", fontsize=25, fontname="Times New Roman", weight = 'bold')
			else:
				panel2.set_ylabel("Number of Single Base Substitutions", fontsize=25, fontname="Times New Roman", weight = 'bold')


			panel1.axis('off')
			panel1.tick_params(axis='both',which='both',\
							   bottom=False, labelbottom=False,\
							   left=False, labelleft=False,\
							   right=False, labelright=False,\
							   top=False, labeltop=False,\
							   direction='in', length=25, colors='white', width=2)
			panel2.tick_params(axis='both',which='both',\
							   bottom=False, labelbottom=False,\
							   left=False, labelleft=False,\
							   right=True, labelright=False,\
							   top=False, labeltop=False,\
							   direction='in', length=25, colors='white', width=2)
			panel3.axis('off')
			panel3.tick_params(axis='both',which='both',\
							   bottom=False, labelbottom=False,\
							   left=False, labelleft=False,\
							   right=True, labelright=False,\
							   top=False, labeltop=False,\
							   direction='in', length=25, colors='white', width=2)
			panel4.axis('off')
			panel4.tick_params(axis='both',which='both',\
							   bottom=False, labelbottom=False,\
							   left=False, labelleft=False,\
							   right=True, labelright=False,\
							   top=False, labeltop=False,\
							   direction='in', length=25, colors='white', width=2)

			buffer = io.BytesIO()
			plt.savefig(buffer, format="png")
			buff_list[sample]=buffer
			plt.close()
			sample_count += 1
		return buff_list

	else:
		print("The provided plot_type: ", plot_type, " is not supported by this plotting function")


def plotID(matrix_path, output_path, project, plot_type, percentage=False, custom_text_upper=None, custom_text_middle=None, custom_text_bottom=None):
	# if 'roman' in matplotlib.font_manager.weight_dict:
	# 	del matplotlib.font_manager.weight_dict['roman']
	# 	matplotlib.font_manager._rebuild()
	plot_custom_text = False
	sig_probs = False
	pcawg = False
	buff_list = dict()
	if plot_type == '94' or plot_type == 'ID94' or plot_type == '94ID' or plot_type == '83':
		first_line=matrix_path.iloc[0,:]
		if first_line[0][1] == 'D' or first_line[0][0] == 'D':
			pcawg = True
		mutation_type = first_line[0]
		mutation_type_list = mutation_type.split(":")
		if len(mutation_type_list) != 4 and first_line[0][1] != 'D' and first_line[0][0] != 'D':
			sys.exit("The matrix does not match the correct ID96 format. Please check you formatting and rerun this plotting function.")
		#pp = PdfPages(output_path + 'ID_83_plots_' + project + '.pdf')

		indel_types = ['1:Del:C:0', '1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5',
					   '1:Del:T:0', '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5',
					   '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
					   '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5',
							# >1bp INDELS
					   '2:Del:R:0', '2:Del:R:1', '2:Del:R:2', '2:Del:R:3', '2:Del:R:4', '2:Del:R:5',
					   '3:Del:R:0', '3:Del:R:1', '3:Del:R:2', '3:Del:R:3', '3:Del:R:4', '3:Del:R:5',
					   '4:Del:R:0', '4:Del:R:1', '4:Del:R:2', '4:Del:R:3', '4:Del:R:4', '4:Del:R:5',
					   '5:Del:R:0', '5:Del:R:1', '5:Del:R:2', '5:Del:R:3', '5:Del:R:4', '5:Del:R:5',
					   '2:Ins:R:0', '2:Ins:R:1', '2:Ins:R:2', '2:Ins:R:3', '2:Ins:R:4', '2:Ins:R:5',
					   '3:Ins:R:0', '3:Ins:R:1', '3:Ins:R:2', '3:Ins:R:3', '3:Ins:R:4', '3:Ins:R:5',
					   '4:Ins:R:0', '4:Ins:R:1', '4:Ins:R:2', '4:Ins:R:3', '4:Ins:R:4', '4:Ins:R:5',
					   '5:Ins:R:0', '5:Ins:R:1', '5:Ins:R:2', '5:Ins:R:3', '5:Ins:R:4', '5:Ins:R:5',
							#MicroHomology INDELS
					   '2:Del:M:1', '3:Del:M:1', '3:Del:M:2', '4:Del:M:1', '4:Del:M:2', '4:Del:M:3',
					   '5:Del:M:1', '5:Del:M:2', '5:Del:M:3', '5:Del:M:4', '5:Del:M:5']
		mutations = OrderedDict()

		first_line=matrix_path.iloc[0,:]
		if pcawg:
			samples=matrix_path.columns[:]
			samples = samples[4:]
			samples = [x.replace('"','') for x in samples]
		else:
			samples=matrix_path.columns[:]
			samples = samples[1:]
		for sample in samples:
			mutations[sample] = OrderedDict()
			mutations[sample]['1DelC'] = [0,0,0,0,0,0]
			mutations[sample]['1DelT'] = [0,0,0,0,0,0]
			mutations[sample]['1InsC'] = [0,0,0,0,0,0]
			mutations[sample]['1InsT'] = [0,0,0,0,0,0]
			mutations[sample]['2DelR'] = [0,0,0,0,0,0]
			mutations[sample]['3DelR'] = [0,0,0,0,0,0]
			mutations[sample]['4DelR'] = [0,0,0,0,0,0]
			mutations[sample]['5DelR'] = [0,0,0,0,0,0]
			mutations[sample]['2InsR'] = [0,0,0,0,0,0]
			mutations[sample]['3InsR'] = [0,0,0,0,0,0]
			mutations[sample]['3InsR'] = [0,0,0,0,0,0]
			mutations[sample]['4InsR'] = [0,0,0,0,0,0]
			mutations[sample]['5InsR'] = [0,0,0,0,0,0]
			mutations[sample]['2DelM'] = [0]
			mutations[sample]['3DelM'] = [0,0]
			mutations[sample]['4DelM'] = [0,0,0]
			mutations[sample]['5DelM'] = [0,0,0,0,0]

		for lines_tmp in range(0,matrix_path.shape[0]):
			if pcawg:
				line = matrix_path.iloc[lines_tmp,:]
				line = [x.replace('"','') for x in line]
				if line[1] == 'repeats':
					mut_type = line[2][0] + line[0][0] + line[0][1].lower() + line[0][2].lower() + "R"
				else:
					mut_type = line[2][0] + line[0][0] + line[0][1].lower() + line[0][2].lower() + line[1][0]
				try:
					repeat_size = int(line[3])
				except:
					repeat_size = int(line[3][0])
				if line[1] == 'MH':
					repeat_size -= 1
				sample_index = 4
			else:
				line = matrix_path.iloc[lines_tmp,:]
				if line[0] not in indel_types:
					continue
				categories = line[0].split(":")
				mut_type = categories[0] + categories[1] + categories[2]
				repeat_size = int(categories[3])
				if categories[2] == 'M':
					repeat_size -= 1
				sample_index = 1

			for sample in samples:
				if mut_type in mutations[sample].keys():
					if percentage:
						mutCount = float(line[sample_index])
						if mutCount < 1 and mutCount > 0:
							sig_probs = True
					else:
						mutCount = int(line[sample_index])
					mutations[sample][mut_type][repeat_size] = mutCount
				else:
					continue
				sample_index += 1

		sample_count = 0
		for sample in mutations.keys():
			total_count = sum(sum(nuc) for nuc in mutations[sample].values())
			plt.rcParams['axes.linewidth'] = 2
			plot1 = plt.figure(figsize=(43.93,12))
			plt.rc('axes', edgecolor='black')
			panel1 = plt.axes([0.045, 0.17, 0.92, 0.65])
			xlabels = []

			x = 0.4
			ymax = 0
			colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256],
					  [253/256,202/256,181/256], [252/256,138/256,106/256], [241/256,68/256,50/256], [188/256,25/256,26/256],
					  [208/256,225/256,242/256], [148/256,196/256,223/256], [74/256,152/256,201/256], [23/256,100/256,171/256],
					  [226/256,226/256,239/256], [182/256,182/256,216/256], [134/256,131/256,189/256], [98/256,64/256,155/256]]

			i = 0
			for key in mutations[sample]:
				l = 1
				for seq in mutations[sample][key]:
					xlabels.append(l)
					if percentage:
						if total_count > 0:
							plt.bar(x, seq/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
							if seq/total_count*100 > ymax:
								ymax = seq/total_count*100
					else:
						plt.bar(x, seq,width=0.4,color=colors[i],align='center', zorder=1000)
						if seq > ymax:
								ymax = seq
					x += 1
					l += 1
				i += 1

			x = .0475
			y_top = .827
			y_bottom = .114
			y = int(ymax*1.25)
			y2 = y+2
			for i in range(0, 12, 1):
				panel1.add_patch(plt.Rectangle((x,y_top), .0595, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
				panel1.add_patch(plt.Rectangle((x,y_bottom), .0595, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
				x += .0665

			panel1.add_patch(plt.Rectangle((x-.001,y_top), .006, .05, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((x-.001,y_bottom), .006, .05, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
			x +=.011
			panel1.add_patch(plt.Rectangle((x,y_top), .0155, .05, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((x,y_bottom), .0155, .05, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
			x += .022
			panel1.add_patch(plt.Rectangle((x,y_top), .027, .05, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((x,y_bottom), .027, .05, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
			x += .0335
			panel1.add_patch(plt.Rectangle((x,y_top), .049, .05, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((x,y_bottom), .049, .05, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))




			yText = y_top + .01
			plt.text(.072, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
			plt.text(.1385, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
			plt.text(.205, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
			plt.text(.2715, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
			plt.text(.338, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.4045, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.471, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.5375, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			plt.text(.604, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.6705, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.737, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.8035, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			plt.text(.844, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.861, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.888, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.93, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

			yText_labels_top = yText + .075
			yText_labels_bottom = y_bottom - .03
			yText_labels_bottom_sec = yText_labels_bottom - .045

			plt.text(.08, yText_labels_top, '1bp Deletion', fontsize=40, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
			plt.text(.21, yText_labels_top, '1bp Insertion', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.375, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.64, yText_labels_top, '>1bp Insertions at Repeats\n       (Insertion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.85, yText_labels_top, ' Mircohomology\n(Deletion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

			plt.text(.058, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=35, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
			plt.text(.19, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.39, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.65, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			plt.text(.85, yText_labels_bottom_sec, 'Mircohomology Length', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

			x = .0477
			for i in range (0, 8, 1):
				if i != 2 and i != 3:
					plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				else:
					plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

				x += .0665

			for i in range (0, 4, 1):
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				x += .0665

			plt.text(x, yText_labels_bottom, '1', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			x += .011
			plt.text(x, yText_labels_bottom, '1  2', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			x += .022
			plt.text(x, yText_labels_bottom, '1  2  3', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			x += .0335
			plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)


			if y <= 4:
				y += 4

			while y%4 != 0:
				y += 1
			ytick_offest = int(y/4)

			if percentage:
				ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
				ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%",
						  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
			else:
				ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
				ylabels= [0, ytick_offest, ytick_offest*2,
						  ytick_offest*3, ytick_offest*4]

			labs = np.arange(0.375,83.375,1)

			if not percentage:
				ylabels = ['{:,}'.format(int(x)) for x in ylabels]
				if len(ylabels[-1]) > 3:
					ylabels_temp = []
					if len(ylabels[-1]) > 7:
						for label in ylabels:
							if len(label) > 7:
								ylabels_temp.append(label[0:-8] + "m")
							elif len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)

					else:
						for label in ylabels:
							if len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)
					ylabels = ylabels_temp

			panel1.set_xlim([0, 83])
			panel1.set_ylim([0, y])
			panel1.set_xticks(labs)
			panel1.set_yticks(ylabs)

			if sig_probs:
				plt.text(0.0475, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
			else:
				plt.text(0.0475, 0.75, sample + ": " + "{:,}".format(int(total_count)) + " indels", fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

			custom_text_upper_plot = ''
			try:
				custom_text_upper[sample_count]
			except:
				custom_text_upper = False
			try:
				custom_text_middle[sample_count]
			except:
				custom_text_middle = False
			try:
				custom_text_bottom[sample_count]
			except:
				custom_text_bottom = False

			if custom_text_upper:
				plot_custom_text = True
				if len(custom_text_upper[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False
			if custom_text_middle:
				if len(custom_text_middle[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False

			if plot_custom_text:
				x_pos_custom = 0.95
				if custom_text_upper and custom_text_middle:
					custom_text_upper_plot = custom_text_upper[sample_count] + "\n" + custom_text_middle[sample_count]
					if custom_text_bottom:
						custom_text_upper_plot += "\n" + custom_text_bottom[sample_count]

				if custom_text_upper and not custom_text_middle:
					custom_text_upper_plot = custom_text_upper[sample_count]
					panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

				elif custom_text_upper and custom_text_middle:
					if not custom_text_bottom:
						panel1.text(x_pos_custom, 0.72, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')
					else:
						panel1.text(x_pos_custom, 0.68, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

				elif not custom_text_upper and custom_text_middle:
					custom_text_upper_plot = custom_text_middle[sample_count]
					panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')




			panel1.set_yticklabels(ylabels, fontsize=30)
			plt.gca().yaxis.grid(True)
			plt.gca().grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
			panel1.set_xlabel('')
			panel1.set_ylabel('')

			if percentage:
				plt.ylabel("Percentage of Indels", fontsize=35, fontname="Times New Roman", weight = 'bold')
			else:
				plt.ylabel("Number of Indels", fontsize=35, fontname="Times New Roman", weight = 'bold')

			panel1.tick_params(axis='both',which='both',\
							   bottom=False, labelbottom=False,\
							   left=False, labelleft=True,\
							   right=False, labelright=False,\
							   top=False, labeltop=False,\
							   direction='in', length=25, colors='gray', width=2)

			[i.set_color("black") for i in plt.gca().get_yticklabels()]

			buffer = io.BytesIO()
			plt.savefig(buffer, format="png")
			buff_list[sample]=buffer
			plt.close()
			sample_count += 1
		return buff_list

	# =======================================OLD INDEL-TSB PLOT==================================================================================
	# elif plot_type == '96' or plot_type == 'ID96' or plot_type == '96ID' or plot_type == 'IDSB':
	# 	with open(matrix_path) as f:
	# 		next(f)
	# 		first_line = f.readline()
	# 		first_line = first_line.strip().split("\t")
	# 		mutation_type = first_line[0]
	# 		mutation_type_list = mutation_type.split(":")
	# 		if len(mutation_type_list) != 5:
	# 			print(mutation_type_list)
	# 			sys.exit("The matrix does not match the correct ID-96 format. Please check you formatting and rerun this plotting function.")

	# 	pp = PdfPages(output_path + 'ID_96_plots_' + project + '.pdf')

	# 	indel_types = ['1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5', '1:Del:C:6'
	# 				   '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5', '1:Del:T:6'
	# 				   '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
	# 				   '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5']

	# 	sig_probs = False
	# 	mutations = OrderedDict()
	# 	with open (matrix_path) as f:
	# 		first_line = f.readline()
	# 		samples = first_line.strip().split("\t")
	# 		samples = samples[1:]
	# 		for sample in samples:
	# 			mutations[sample] = OrderedDict()
	# 			mutations[sample]['1DelC'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
	# 			mutations[sample]['1DelT'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
	# 			mutations[sample]['1InsC'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
	# 			mutations[sample]['1InsT'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]

	# 		for lines in f:
	# 			line = lines.strip().split()
	# 			categories = line[0].split(":")
	# 			if categories[1] != '1':
	# 				break
	# 			else:
	# 				mut_type = categories[1] + categories[2] + categories[3]
	# 				bias = categories[0]
	# 				if bias != 'T' and bias != 'U':
	# 					continue
	# 				else:
	# 					repeat_size = int(categories[4])
	# 					sample_index = 1
	# 					for sample in samples:
	# 						if mut_type in mutations[sample]:
	# 							if percentage:
	# 								mutCount = float(line[sample_index])
	# 								if mutCount < 1 and mutCount > 0:
	# 									sig_probs = True
	# 							else:
	# 								mutCount = int(line[sample_index])
	# 							if bias == 'T':
	# 								mutations[sample][mut_type][repeat_size][0] = mutCount
	# 							else:
	# 								mutations[sample][mut_type][repeat_size][1] = mutCount
	# 						else:
	# 							continue
	# 						sample_index += 1

	# 	for sample in mutations:
	# 		total_count = sum(sum(sum(tsb) for tsb in nuc) for nuc in mutations[sample].values())
	# 		plt.rcParams['axes.linewidth'] = 2
	# 		plot1 = plt.figure(figsize=(15,13))
	# 		plt.rc('axes', edgecolor='black')
	# 		panel1 = plt.axes([0.12, 0.12, 0.8, 0.77])
	# 		xlabels = []

	# 		x = 0.3
	# 		ymax = 0
	# 		colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256]]

	# 		i = 0
	# 		for key in mutations[sample]:
	# 			l = 1
	# 			for seq in mutations[sample][key]:
	# 				xlabels.append(l)
	# 				if percentage:
	# 					trans = plt.bar(x, seq[0]/total_count*100,width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
	# 					x += 0.2
	# 					untrans = plt.bar(x, seq[1]/total_count*100,width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
	# 					x += 0.8
	# 					if seq[0]/total_count*100 > ymax:
	# 						ymax = seq[0]/total_count*100
	# 					if seq[1]/total_count*100 > ymax:
	# 						ymax = seq[1]/total_count*100
	# 				else:
	# 					trans = plt.bar(x, seq[0],width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
	# 					x += 0.2
	# 					untrans = plt.bar(x, seq[1],width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
	# 					x += 0.8
	# 					if seq[0] > ymax:
	# 							ymax = seq[0]
	# 					if seq[1] > ymax:
	# 							ymax = seq[1]
	# 				l += 1
	# 			i += 1

	# 		x = .125
	# 		y_top = .8975
	# 		#y_bottom = .06525
	# 		y_bottom = .075
	# 		y = int(ymax*1.25)
	# 		y2 = y+2





	# 		for i in range(0, 4, 1):
	# 			panel1.add_patch(plt.Rectangle((x,y_top), .185, .037, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
	# 			panel1.add_patch(plt.Rectangle((x,y_bottom), .185, .037, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
	# 			x += .202


	# 		yText = y_top + .005
	# 		plt.text(.205, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	# 		plt.text(.407, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	# 		plt.text(.609, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	# 		plt.text(.811, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)

	# 		yText_labels_top = yText + .05
	# 		yText_labels_bottom = y_bottom - .03
	# 		yText_labels_bottom_sec = yText_labels_bottom -.025

	# 		plt.text(.23, yText_labels_top, '1bp Deletion', fontsize=35, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	# 		plt.text(.634, yText_labels_top, '1bp Insertion', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	# 		plt.text(.18, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=30, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	# 		plt.text(.58, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=30, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	# 		x = .127
	# 		yText_labels_bottom = y_bottom - 0.025

	# 		for l in range (0, 4, 1):
	# 			if l < 2:
	# 				for i in range(1, 6, 1):
	# 					plt.text(x, yText_labels_bottom, str(i), fontsize=25, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	# 					x += 0.034
	# 				x -= 0.008
	# 				plt.text(x, yText_labels_bottom, '6+', fontsize=25, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	# 				x += 0.041
	# 			else:
	# 				for i in range(0, 5, 1):
	# 					plt.text(x, yText_labels_bottom, str(i), fontsize=25, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	# 					x += 0.033
	# 				x -= 0.005
	# 				plt.text(x, yText_labels_bottom, '5+', fontsize=25, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	# 				x += 0.044

	# 		if y <= 4:
	# 			y += 4

	# 		while y%4 != 0:
	# 			y += 1
	# 		ytick_offest = int(y/4)


	# 		x_shaded = 0
	# 		for i in range(0, 4, 1):
	# 			panel1.add_patch(plt.Rectangle((x_shaded,0), 6, y, facecolor=colors[i], zorder=0, alpha = 0.25, edgecolor='grey'))
	# 			x_shaded += 6

	# 		if percentage:
	# 			ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
	# 			ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%",
	# 					  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
	# 		else:
	# 			ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
	# 			ylabels= [0, ytick_offest, ytick_offest*2,
	# 					  ytick_offest*3, ytick_offest*4]


	# 		if not percentage:
	# 			ylabels = ['{:,}'.format(int(x)) for x in ylabels]
	# 			if len(ylabels[-1]) > 3:
	# 				ylabels_temp = []
	# 				if len(ylabels[-1]) > 7:
	# 					for label in ylabels:
	# 						if len(label) > 7:
	# 							ylabels_temp.append(label[0:-8] + "m")
	# 						elif len(label) > 3:
	# 							ylabels_temp.append(label[0:-4] + "k")
	# 						else:
	# 							ylabels_temp.append(label)

	# 				else:
	# 					for label in ylabels:
	# 						if len(label) > 3:
	# 							ylabels_temp.append(label[0:-4] + "k")
	# 						else:
	# 							ylabels_temp.append(label)
	# 				ylabels = ylabels_temp

	# 		panel1.set_xlim([0, 23.8])
	# 		panel1.set_ylim([0, y])
	# 		panel1.set_yticks(ylabs)

	# 		if sig_probs:
	# 			plt.text(0.13, 0.85, sample, fontsize=33, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
	# 		else:
	# 			plt.text(0.13, 0.85, sample + ": " + "{:,}".format(int(total_count)) + " transcribed indels", fontsize=33, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

	# 		panel1.set_yticklabels(ylabels, fontsize=30)
	# 		plt.gca().yaxis.grid(True)
	# 		plt.gca().grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
	# 		panel1.set_xlabel('')
	# 		panel1.set_ylabel('')
	# 		plt.legend(handles=[trans, untrans], prop={'size':20})

	# 		if percentage:
	# 			plt.ylabel("Percentage of Indels", fontsize=35, fontname="Times New Roman", weight = 'bold')
	# 		else:
	# 			plt.ylabel("Number of Indels", fontsize=35, fontname="Times New Roman", weight = 'bold')

	# 		panel1.tick_params(axis='both',which='both',\
	# 						   bottom=False, labelbottom=False,\
	# 						   left=False, labelleft=True,\
	# 						   right=False, labelright=False,\
	# 						   top=False, labeltop=False,\
	# 						   direction='in', length=25, colors='gray', width=2)

	# 		[i.set_color("black") for i in plt.gca().get_yticklabels()]





	# 		pp.savefig(plot1)
	# 		plt.close()
	# 	pp.close()
	# ===================================================================================================================================



	elif plot_type == 'INDEL_simple' or plot_type == 'simple_INDEL' or plot_type == 'ID_simple' or plot_type == 'simple_ID' or plot_type == '28':
		with open(matrix_path) as f:
			next(f)
			first_line = f.readline()
			first_line = first_line.strip().split()
			mutation_type = first_line[0]
			mutation_type_list = mutation_type.split(":")
			if len(mutation_type_list) != 4:
				sys.exit("The matrix does not match the correct SBS96 format. Please check you formatting and rerun this plotting function.")
		pp = PdfPages(output_path + 'ID_simple_plots_' + project + '.pdf')

		indel_types = ['1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5', '1:Del:C:6'
					   '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5', '1:Del:T:6'
					   '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
					   '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5',
					   'long_Del', 'long_Ins', 'MH', 'complex']

		mutations = OrderedDict()

		try:
			with open (matrix_path) as f:
				first_line = f.readline()
				samples = first_line.strip().split("\t")
				samples = samples[1:]
				for sample in samples:
					mutations[sample] = OrderedDict()
					mutations[sample]['1DelC'] = [0,0,0,0,0,0]
					mutations[sample]['1DelT'] = [0,0,0,0,0,0]
					mutations[sample]['1InsC'] = [0,0,0,0,0,0]
					mutations[sample]['1InsT'] = [0,0,0,0,0,0]
					mutations[sample]['long_Del'] = [0]
					mutations[sample]['long_Ins'] = [0]
					mutations[sample]['MH'] = [0]
					mutations[sample]['complex'] = [0]

				for lines in f:
					line = lines.strip().split()
					categories = line[0].split(":")
					if len(categories) < 2:
						mut_type = categories[0]
						repeat_size = 0
					else:
						mut_type = categories[0] + categories[1] + categories[2]
						repeat_size = int(categories[3])
					sample_index = 1

					for sample in samples:
						#if mut_type in mutations[sample].keys():
							if percentage:
								mutCount = float(line[sample_index])
								if mutCount < 1 and mutCount > 0:
									sig_probs = True
							else:
								mutCount = int(line[sample_index])
							mutations[sample][mut_type][repeat_size] = mutCount

						# else:
						# 	if percentage:
						# 		mutCount = float(line[sample_index])
						# 	else:
						# 		mutCount = int(line[sample_index])
						# 	if int(mut_type[0]) > 1:
						# 		repeat_size = 0
						# 		if categories[2] == 'M':
						# 			mut_type = 'MH'
						# 			mutations[sample][mut_type][repeat_size] += mutCount
						# 		else:
						# 			if categories[1] == 'Del':
						# 				mut_type = 'Del'
						# 				mutations[sample][mut_type][repeat_size] += mutCount
						# 			else:
						# 				mut_type = 'Ins'
						# 				mutations[sample][mut_type][repeat_size] += mutCount

							#continue
							sample_index += 1

			for sample in mutations:
				total_count = sum(sum(nuc) for nuc in mutations[sample].values())
				plt.rcParams['axes.linewidth'] = 2
				plot1 = plt.figure(figsize=(15,13))
				plt.rc('axes', edgecolor='black')
				panel1 = plt.axes([0.12, 0.12, 0.8, 0.77])
				xlabels = []

				x = 0.4
				ymax = 0
				colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256],
						  #[188/256,25/256,26/256],
						  [23/256,100/256,171/256],[98/256,64/256,155/256], [98/256,64/256,155/256]]

				i = 0
				for key in mutations[sample]:
					l = 1
					for seq in mutations[sample][key]:
						xlabels.append(l)
						if percentage:
							if total_count > 0:
								plt.bar(x, seq/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
								if seq/total_count*100 > ymax:
									ymax = seq/total_count*100
						else:
							plt.bar(x, seq,width=0.4,color=colors[i],align='center', zorder=1000)
							if seq > ymax:
									ymax = seq
						x += 1
						l += 1
					if i < 4:
						i += 1
				x = .126
				y_top = .9
				y_bottom = .075
				y = int(ymax*1.25)
				y2 = y+2

				for i in range(0, 4, 1):
					panel1.add_patch(plt.Rectangle((x,y_top), .154, .037, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
					panel1.add_patch(plt.Rectangle((x,y_bottom), .154, .037, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
					x += .1715

				x -= .001
				panel1.add_patch(plt.Rectangle((x,y_top), .098, .037, facecolor=colors[i+1], clip_on=False, transform=plt.gcf().transFigure))
				panel1.add_patch(plt.Rectangle((x,y_bottom), .098, .037, facecolor=colors[i+1], clip_on=False, transform=plt.gcf().transFigure))

				yText = y_top + .0055
				plt.text(.185, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
				plt.text(.36, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
				plt.text(.53, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
				plt.text(.705, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)

				yText_labels_top = yText + .045
				yText_labels_bottom = y_bottom - .03
				yText_labels_bottom_sec = yText_labels_bottom -.025

				plt.text(.2, yText_labels_top, '1bp Deletion', fontsize=35, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
				plt.text(.54, yText_labels_top, '1bp Insertion', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.155, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=30, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
				plt.text(.505, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=30, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.827, yText_labels_top, '>1bp', fontsize=30, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.83, yText_labels_bottom_sec, 'Type', fontsize=30, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

				x = .127
				yText_labels_bottom = y_bottom - 0.025

				for l in range (0, 4, 1):
					if l < 2:
						for i in range(1, 6, 1):
							plt.text(x, yText_labels_bottom, str(i), fontsize=25, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
							x += 0.028
						x -= 0.005
						plt.text(x, yText_labels_bottom, '6+', fontsize=25, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
						x += 0.037
					else:
						if l == 2:
							x += 0
						for i in range(0, 5, 1):
							plt.text(x, yText_labels_bottom, str(i), fontsize=25, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
							x += 0.028
						x -= 0.005
						plt.text(x, yText_labels_bottom, '5+', fontsize=25, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
						x += 0.037

				yText_labels_bottom += 0.01
				plt.text(x, yText_labels_bottom, 'Del', fontsize=17, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure, rotation='vertical')
				x += 0.026
				plt.text(x, yText_labels_bottom, 'Ins', fontsize=17, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure, rotation='vertical')
				x += 0.0295
				yText_labels_bottom += 0.003
				plt.text(x, yText_labels_bottom, 'MH', fontsize=17, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure, rotation='vertical')
				x += 0.0295
				yText_labels_bottom += 0.005
				plt.text(x, yText_labels_bottom, 'COMP', fontsize=10, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure, rotation='vertical')

				if y <= 4:
					y += 4

				while y%4 != 0:
					y += 1
				ytick_offest = int(y/4)

				if percentage:
					ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
					ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%",
							  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
				else:
					ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
					ylabels= [0, ytick_offest, ytick_offest*2,
							  ytick_offest*3, ytick_offest*4]


				if not percentage:
					ylabels = ['{:,}'.format(int(x)) for x in ylabels]
					if len(ylabels[-1]) > 3:
						ylabels_temp = []
						if len(ylabels[-1]) > 7:
							for label in ylabels:
								if len(label) > 7:
									ylabels_temp.append(label[0:-8] + "m")
								elif len(label) > 3:
									ylabels_temp.append(label[0:-4] + "k")
								else:
									ylabels_temp.append(label)

						else:
							for label in ylabels:
								if len(label) > 3:
									ylabels_temp.append(label[0:-4] + "k")
								else:
									ylabels_temp.append(label)
						ylabels = ylabels_temp




				panel1.set_xlim([0, 28])
				panel1.set_ylim([0, y])
				panel1.set_yticks(ylabs)

				if sig_probs:
					plt.text(0.13, 0.85, sample, fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
				else:
					plt.text(0.13, 0.85, sample + ": " + "{:,}".format(int(total_count)) + " indels", fontsize=40, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

				panel1.set_yticklabels(ylabels, fontsize=30)
				plt.gca().yaxis.grid(True)
				plt.gca().grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
				panel1.set_xlabel('')
				panel1.set_ylabel('')

				if percentage:
					plt.ylabel("Percentage of Indels", fontsize=35, fontname="Times New Roman", weight = 'bold')
				else:
					plt.ylabel("Number of Indels", fontsize=35, fontname="Times New Roman", weight = 'bold')

				panel1.tick_params(axis='both',which='both',\
								   bottom=False, labelbottom=False,\
								   left=False, labelleft=True,\
								   right=False, labelright=False,\
								   top=False, labeltop=False,\
								   direction='in', length=25, colors='gray', width=2)

				[i.set_color("black") for i in plt.gca().get_yticklabels()]


				pp.savefig(plot1)
				plt.close()
			pp.close()

		except:
			print("There may be an issue with the formatting of you matrix file.")
			os.remove(output_path + 'ID_simple_plots_' + project + '.pdf')

	elif plot_type == '96' or plot_type == 'ID96' or plot_type == '96ID' or plot_type == 'IDSB' or plot_type == '415':
		with open(matrix_path) as f:
			next(f)
			first_line = f.readline()
			first_line = first_line.strip().split("\t")
			mutation_type = first_line[0]
			mutation_type_list = mutation_type.split(":")
			if len(mutation_type_list) != 5:
				print(mutation_type_list)
				sys.exit("The matrix does not match the correct ID-96 format. Please check you formatting and rerun this plotting function.")

		pp = PdfPages(output_path + 'ID_TSB_plots_' + project + '.pdf')

		indel_types_tsb = []
		tsb_I = ['T','U','N','B','Q']
		indel_types = ['1:Del:C:0', '1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5',
					   '1:Del:T:0', '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5',
					   '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
					   '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5',
							# >1bp INDELS
					   '2:Del:R:0', '2:Del:R:1', '2:Del:R:2', '2:Del:R:3', '2:Del:R:4', '2:Del:R:5',
					   '3:Del:R:0', '3:Del:R:1', '3:Del:R:2', '3:Del:R:3', '3:Del:R:4', '3:Del:R:5',
					   '4:Del:R:0', '4:Del:R:1', '4:Del:R:2', '4:Del:R:3', '4:Del:R:4', '4:Del:R:5',
					   '5:Del:R:0', '5:Del:R:1', '5:Del:R:2', '5:Del:R:3', '5:Del:R:4', '5:Del:R:5',
					   '2:Ins:R:0', '2:Ins:R:1', '2:Ins:R:2', '2:Ins:R:3', '2:Ins:R:4', '2:Ins:R:5',
					   '3:Ins:R:0', '3:Ins:R:1', '3:Ins:R:2', '3:Ins:R:3', '3:Ins:R:4', '3:Ins:R:5',
					   '4:Ins:R:0', '4:Ins:R:1', '4:Ins:R:2', '4:Ins:R:3', '4:Ins:R:4', '4:Ins:R:5',
					   '5:Ins:R:0', '5:Ins:R:1', '5:Ins:R:2', '5:Ins:R:3', '5:Ins:R:4', '5:Ins:R:5',
							#MicroHomology INDELS
					   '2:Del:M:1', '3:Del:M:1', '3:Del:M:2', '4:Del:M:1', '4:Del:M:2', '4:Del:M:3',
					   '5:Del:M:1', '5:Del:M:2', '5:Del:M:3', '5:Del:M:4', '5:Del:M:5']

		for indels in indel_types:
			for tsbs in tsb_I:
				indel_types_tsb.append(tsbs + ":" + indels)

		sig_probs = False
		mutations = OrderedDict()
		try:
			with open (matrix_path) as f:
				first_line = f.readline()
				if pcawg:
					samples = first_line.strip().split(",")
					samples = samples[4:]
					samples = [x.replace('"','') for x in samples]
				else:
					samples = first_line.strip().split("\t")
					samples = samples[1:]
				for sample in samples:
					mutations[sample] = OrderedDict()
					mutations[sample]['1DelC'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['1DelT'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['1InsC'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['1InsT'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['2DelR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['3DelR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['4DelR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['5DelR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['2InsR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['3InsR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['3InsR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['4InsR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['5InsR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
					mutations[sample]['2DelM'] = [[0,0]]
					mutations[sample]['3DelM'] = [[0,0],[0,0]]
					mutations[sample]['4DelM'] = [[0,0],[0,0],[0,0]]
					mutations[sample]['5DelM'] = [[0,0],[0,0],[0,0],[0,0],[0,0]]

				for lines in f:
					if pcawg:
						line = lines.strip().split(",")
						line = [x.replace('"','') for x in line]
						if line[1] == 'repeats':
							mut_type = line[2][0] + line[0][0] + line[0][1].lower() + line[0][2].lower() + "R"
						else:
							mut_type = line[2][0] + line[0][0] + line[0][1].lower() + line[0][2].lower() + line[1][0]
						try:
							repeat_size = int(line[3])
						except:
							repeat_size = int(line[3][0])
						if line[1] == 'MH':
							repeat_size -= 1
						sample_index = 4
					else:
						line = lines.strip().split()
						if line[0] not in indel_types_tsb:
							continue
						categories = line[0].split(":")
						bias = categories[0]
						if bias == 'B' or bias == 'N' or bias == 'Q':
							continue
						mut_type = categories[1] + categories[2] + categories[3]

						repeat_size = int(categories[4])
						if categories[3] == 'M':
							repeat_size -= 1
						sample_index = 1

					for sample in samples:
						if mut_type in mutations[sample].keys():
							if percentage:
								mutCount = float(line[sample_index])
								if mutCount < 1 and mutCount > 0:
									sig_probs = True
							else:
								mutCount = int(line[sample_index])
							if bias == 'T':
								mutations[sample][mut_type][repeat_size][0] = mutCount
							else:
								mutations[sample][mut_type][repeat_size][1] = mutCount
						else:
							continue
						sample_index += 1

			sample_count = 0
			for sample in mutations.keys():
				total_count = sum(sum(sum(tsb) for tsb in nuc) for nuc in mutations[sample].values())
				plt.rcParams['axes.linewidth'] = 2
				plot1 = plt.figure(figsize=(43.93,12))
				plt.rc('axes', edgecolor='black')
				panel1 = plt.axes([0.045, 0.17, 0.92, 0.65])
				xlabels = []

				x = 0.4
				ymax = 0
				colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256],
						  [253/256,202/256,181/256], [252/256,138/256,106/256], [241/256,68/256,50/256], [188/256,25/256,26/256],
						  [208/256,225/256,242/256], [148/256,196/256,223/256], [74/256,152/256,201/256], [23/256,100/256,171/256],
						  [226/256,226/256,239/256], [182/256,182/256,216/256], [134/256,131/256,189/256], [98/256,64/256,155/256]]

				i = 0
				for key in mutations[sample]:
					l = 1
					for seq in mutations[sample][key]:
						xlabels.append(l)
						if percentage:
							if total_count > 0:
								trans = plt.bar(x, seq[0]/total_count*100,width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
								x += 0.2
								untrans = plt.bar(x, seq[1]/total_count*100,width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
								if seq[0]/total_count*100 > ymax:
										ymax = seq[0]/total_count*100
								if seq[1]/total_count*100 > ymax:
										ymax = seq[1]/total_count*100

						else:
							trans = plt.bar(x, seq[0],width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
							x += 0.2
							untrans = plt.bar(x, seq[1],width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
							if seq[0] > ymax:
									ymax = seq[0]
							if seq[1] > ymax:
									ymax = seq[1]

						x += 0.799
						l += 1
					i += 1

				x = .0475
				y_top = .827
				y_bottom = .114
				y = int(ymax*1.25)
				y2 = y+2
				for i in range(0, 12, 1):
					panel1.add_patch(plt.Rectangle((x,y_top), .0595, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
					panel1.add_patch(plt.Rectangle((x,y_bottom), .0595, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
					x += .0665

				panel1.add_patch(plt.Rectangle((x-.001,y_top), .006, .05, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
				panel1.add_patch(plt.Rectangle((x-.001,y_bottom), .006, .05, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
				x +=.011
				panel1.add_patch(plt.Rectangle((x,y_top), .0155, .05, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
				panel1.add_patch(plt.Rectangle((x,y_bottom), .0155, .05, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
				x += .022
				panel1.add_patch(plt.Rectangle((x,y_top), .027, .05, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
				panel1.add_patch(plt.Rectangle((x,y_bottom), .027, .05, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
				x += .0335
				panel1.add_patch(plt.Rectangle((x,y_top), .049, .05, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
				panel1.add_patch(plt.Rectangle((x,y_bottom), .049, .05, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))




				yText = y_top + .01
				plt.text(.072, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
				plt.text(.1385, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
				plt.text(.205, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
				plt.text(.2715, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
				plt.text(.338, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.4045, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.471, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.5375, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
				plt.text(.604, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.6705, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.737, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.8035, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
				plt.text(.844, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.861, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.888, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.93, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

				yText_labels_top = yText + .075
				yText_labels_bottom = y_bottom - .03
				yText_labels_bottom_sec = yText_labels_bottom - .045

				plt.text(.08, yText_labels_top, '1bp Deletion', fontsize=40, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
				plt.text(.21, yText_labels_top, '1bp Insertion', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.375, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.64, yText_labels_top, '>1bp Insertions at Repeats\n       (Insertion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.85, yText_labels_top, ' Mircohomology\n(Deletion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

				plt.text(.058, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=35, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
				plt.text(.19, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.39, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.65, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				plt.text(.85, yText_labels_bottom_sec, 'Mircohomology Length', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

				x = .0477
				for i in range (0, 8, 1):
					if i != 2 and i != 3:
						plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
					else:
						plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

					x += .0665

				for i in range (0, 4, 1):
					plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
					x += .0665

				plt.text(x, yText_labels_bottom, '1', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				x += .011
				plt.text(x, yText_labels_bottom, '1  2', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				x += .022
				plt.text(x, yText_labels_bottom, '1  2  3', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				x += .0335
				plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)


				if y <= 4:
					y += 4

				while y%4 != 0:
					y += 1
				ytick_offest = int(y/4)

				if percentage:
					ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
					ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%",
							  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
				else:
					ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
					ylabels= [0, ytick_offest, ytick_offest*2,
							  ytick_offest*3, ytick_offest*4]

				labs = np.arange(0.375,83.375,1)

				if not percentage:
					ylabels = ['{:,}'.format(int(x)) for x in ylabels]
					if len(ylabels[-1]) > 3:
						ylabels_temp = []
						if len(ylabels[-1]) > 7:
							for label in ylabels:
								if len(label) > 7:
									ylabels_temp.append(label[0:-8] + "m")
								elif len(label) > 3:
									ylabels_temp.append(label[0:-4] + "k")
								else:
									ylabels_temp.append(label)

						else:
							for label in ylabels:
								if len(label) > 3:
									ylabels_temp.append(label[0:-4] + "k")
								else:
									ylabels_temp.append(label)
						ylabels = ylabels_temp

				panel1.set_xlim([0, 83])
				panel1.set_ylim([0, y])
				panel1.set_xticks(labs)
				panel1.set_yticks(ylabs)

				if sig_probs:
					plt.text(0.0475, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
				else:
					plt.text(0.0475, 0.75, sample + ": " + "{:,}".format(int(total_count)) + " indels", fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

				custom_text_upper_plot = ''
				try:
					custom_text_upper[sample_count]
				except:
					custom_text_upper = False
				try:
					custom_text_middle[sample_count]
				except:
					custom_text_middle = False
				try:
					custom_text_bottom[sample_count]
				except:
					custom_text_bottom = False

				if custom_text_upper:
					plot_custom_text = True
					if len(custom_text_upper[sample_count]) > 40:
						print("To add a custom text, please limit the string to <40 characters including spaces.")
						plot_custom_text = False
				if custom_text_middle:
					if len(custom_text_middle[sample_count]) > 40:
						print("To add a custom text, please limit the string to <40 characters including spaces.")
						plot_custom_text = False

				if plot_custom_text:
					x_pos_custom = 0.95
					if custom_text_upper and custom_text_middle:
						custom_text_upper_plot = custom_text_upper[sample_count] + "\n" + custom_text_middle[sample_count]
						if custom_text_bottom:
							custom_text_upper_plot += "\n" + custom_text_bottom[sample_count]

					if custom_text_upper and not custom_text_middle:
						custom_text_upper_plot = custom_text_upper[sample_count]
						panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

					elif custom_text_upper and custom_text_middle:
						if not custom_text_bottom:
							panel1.text(x_pos_custom, 0.72, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')
						else:
							panel1.text(x_pos_custom, 0.68, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

					elif not custom_text_upper and custom_text_middle:
						custom_text_upper_plot = custom_text_middle[sample_count]
						panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')




				panel1.set_yticklabels(ylabels, fontsize=30)
				plt.gca().yaxis.grid(True)
				plt.gca().grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
				panel1.set_xlabel('')
				panel1.set_ylabel('')
				panel1.legend(handles=[trans, untrans], prop={'size':30})


				if percentage:
					plt.ylabel("Percentage of Indels", fontsize=35, fontname="Times New Roman", weight = 'bold')
				else:
					plt.ylabel("Number of Indels", fontsize=35, fontname="Times New Roman", weight = 'bold')

				panel1.tick_params(axis='both',which='both',\
								   bottom=False, labelbottom=False,\
								   left=False, labelleft=True,\
								   right=False, labelright=False,\
								   top=False, labeltop=False,\
								   direction='in', length=25, colors='gray', width=2)

				[i.set_color("black") for i in plt.gca().get_yticklabels()]

				pp.savefig(plot1)
				plt.close()
				sample_count += 1
			pp.close()

		except:
			print("There may be an issue with the formatting of you matrix file.")
			os.remove(output_path + 'ID_TSB_plots_' + project + '.pdf')



	else:
		print("The provided plot_type: ", plot_type, " is not supported by this plotting function")

def plotDBS(matrix_path, output_path, project, plot_type, percentage=False, custom_text_upper=None, custom_text_middle=None, custom_text_bottom=None):

	# if 'roman' in matplotlib.font_manager.weight_dict:
	# 	del matplotlib.font_manager.weight_dict['roman']
	# 	matplotlib.font_manager._rebuild()
	plot_custom_text = False
	pcawg = False
	sig_probs = False
	buff_list = dict()

	if plot_type == '78' or plot_type == '78DBS' or plot_type == 'DBS78':

		first_line=matrix_path.iloc[0,:]
		mutation_type = first_line[0]
		if first_line[0][2] != ">":
			pcawg = True
		if len(mutation_type) != 5 and first_line[0][2] == ">":
			sys.exit("The matrix does not match the correct DBS96 format. Please check you formatting and rerun this plotting function.")

		dinucs = ['TT>GG','TT>CG','TT>AG','TT>GC','TT>CC','TT>AC','TT>GA','TT>CA','TT>AA','AC>CA','AC>CG','AC>CT','AC>GA',
				  'AC>GG','AC>GT','AC>TA','AC>TG','AC>TT','CT>AA','CT>AC','CT>AG','CT>GA','CT>GC','CT>GG','CT>TG','CT>TC',
				  'CT>TA','AT>CA','AT>CC','AT>CG','AT>GA','AT>GC','AT>TA','TG>GT','TG>CT','TG>AT','TG>GC','TG>CC','TG>AC',
				  'TG>GA','TG>CA','TG>AA','CC>AA','CC>AG','CC>AT','CC>GA','CC>GG','CC>GT','CC>TA','CC>TG','CC>TT','CG>AT',
				  'CG>GC','CG>GT','CG>TC','CG>TA','CG>TT','TC>GT','TC>CT','TC>AT','TC>GG','TC>CG','TC>AG','TC>GA','TC>CA',
				  'TC>AA','GC>AA','GC>AG','GC>AT','GC>CA','GC>CG','GC>TA','TA>GT','TA>CT','TA>AT','TA>GG','TA>CG','TA>GC']

		revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
		mutations = OrderedDict()

		first_line=matrix_path.iloc[0,:]
		if pcawg:
			samples = matrix_path.columns[:]
			samples = samples[2:]
			samples = [x.replace('"','') for x in samples]
		else:
			samples = matrix_path.columns[:]
			samples = samples[1:]
		for sample in samples:
			mutations[sample] = OrderedDict()
			mutations[sample]['AC'] = OrderedDict()
			mutations[sample]['AT'] = OrderedDict()
			mutations[sample]['CC'] = OrderedDict()
			mutations[sample]['CG'] = OrderedDict()
			mutations[sample]['CT'] = OrderedDict()
			mutations[sample]['GC'] = OrderedDict()
			mutations[sample]['TA'] = OrderedDict()
			mutations[sample]['TC'] = OrderedDict()
			mutations[sample]['TG'] = OrderedDict()
			mutations[sample]['TT'] = OrderedDict()


		for lines_tmp in range(0,matrix_path.shape[0]):
			if pcawg:
				line = matrix_path.iloc[lines_tmp,:]
				line = [x.replace('"','') for x in line]
				mut = line[0] + ">" + line[1]
				nuc = line[1]
				mut_type = line[0]
				if mut not in dinucs:
					nuc = revcompl(nuc)
					mut_type = revcompl(mut_type)
				sample_index = 2
			else:
				line = matrix_path.iloc[lines_tmp,:]
				mut = line[0]
				nuc = line[0][3:]
				mut_type = line[0][0:2]
				if mut not in dinucs:
					nuc = revcompl(nuc)
					mut_type = revcompl(mut_type)
				sample_index = 1

			for sample in samples:
				if percentage:
					mutCount = float(line[sample_index])
					if mutCount < 1 and mutCount > 0:
						sig_probs = True
				else:
					mutCount = int(line[sample_index])
				mutations[sample][mut_type][nuc] = mutCount
				sample_index += 1

		sample_count = 0
		for sample in mutations.keys():
			total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
			plt.rcParams['axes.linewidth'] = 4
			plot1 = plt.figure(figsize=(43.93,9.92))
			plt.rc('axes', edgecolor='grey')
			panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
			xlabels = []

			x = 0.4
			ymax = 0
			colors = [[3/256,189/256,239/256], [3/256,102/256,204/256],[162/256,207/256,99/256],
					  [1/256,102/256,1/256], [255/256,153/256,153/256], [228/256,41/256,38/256],
					  [255/256,178/256,102/256], [255/256,128/256,1/256], [204/256,153/256,255/256],
					  [76/256,1/256,153/256]]
			i = 0
			for key in mutations[sample]:
				muts = mutations[sample][key].keys()
				muts = sorted(muts)
				for seq in muts:
					xlabels.append(seq)
					if percentage:
						if total_count > 0:
							plt.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
							if mutations[sample][key][seq]/total_count*100 > ymax:
									ymax = mutations[sample][key][seq]/total_count*100
					else:
						plt.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
						if mutations[sample][key][seq] > ymax:
								ymax = mutations[sample][key][seq]
					x += 1
				i += 1

			x = .043
			y3 = .87
			y = ymax*1.25
			y2 = y+2
			i = 0
			panel1.add_patch(plt.Rectangle((.043,y3), .101, .05, facecolor=colors[0], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((.151,y3), .067, .05, facecolor=colors[1], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((.225,y3), .102, .05, facecolor=colors[2], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((.334,y3), .067, .05, facecolor=colors[3], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((.408,y3), .102, .05, facecolor=colors[4], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((.517,y3), .067, .05, facecolor=colors[5], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((.591,y3), .067, .05, facecolor=colors[6], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((.665,y3), .102, .05, facecolor=colors[7], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((.774,y3), .102, .05, facecolor=colors[8], clip_on=False, transform=plt.gcf().transFigure))
			panel1.add_patch(plt.Rectangle((.883,y3), .102, .05, facecolor=colors[9], clip_on=False, transform=plt.gcf().transFigure))

			yText = y3 + .06
			plt.text(.07, yText, 'AC>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.163, yText, 'AT>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.255, yText, 'CC>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.345, yText, 'CG>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.435, yText, 'CT>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.527, yText, 'GC>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.6, yText, 'TA>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.69, yText, 'TC>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.8, yText, 'TG>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.915, yText, 'TT>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

			if y <= 4:
				y += 4

			if percentage:
				ytick_offest = round((y/4), 1)
				ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
				ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%",
						  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
			else:
				# if y < 10:
				# 	if y/4 - int(y/4) > 0.5:
				# 		ytick_offest = int(y/4) + 1
				# 	else:
				# 		ytick_offest = int(y/4)
				if y < 4:
					y = 4
				#else:
				ytick_offest = int(y/4)
				if ytick_offest == 0:
					ytick_offest = 1
				ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
				ylabels= [0, ytick_offest, ytick_offest*2,
						  ytick_offest*3, ytick_offest*4]
			if y < 4:
				y = 4

			if sig_probs:
				plt.text(0.045, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
			else:
				plt.text(0.045, 0.75, sample + ": " + "{:,}".format(int(total_count)) + " double subs", fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

			custom_text_upper_plot = ''
			try:
				custom_text_upper[sample_count]
			except:
				custom_text_upper = False
			try:
				custom_text_middle[sample_count]
			except:
				custom_text_middle = False
			try:
				custom_text_bottom[sample_count]
			except:
				custom_text_bottom = False

			if custom_text_upper:
				plot_custom_text = True
				if len(custom_text_upper[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False
			if custom_text_bottom:
				if len(custom_text_bottom[sample_count]) > 40:
					print("To add a custom text, please limit the string to <40 characters including spaces.")
					plot_custom_text = False

			if plot_custom_text:
				x_pos_custom = 0.98
				if custom_text_upper and custom_text_middle:
					custom_text_upper_plot = custom_text_upper[sample_count] + "\n" + custom_text_middle[sample_count]
					if custom_text_bottom:
						custom_text_upper_plot += "\n" + custom_text_bottom[sample_count]

				if custom_text_upper and not custom_text_middle:
					custom_text_upper_plot = custom_text_upper[sample_count]
					panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

				elif custom_text_upper and custom_text_middle:
					if not custom_text_bottom:
						panel1.text(x_pos_custom, 0.75, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')
					else:
						panel1.text(x_pos_custom, 0.7, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

				elif not custom_text_upper and custom_text_middle:
					custom_text_upper_plot = custom_text_middle[sample_count]
					panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')


			if not percentage:
				ylabels = ['{:,}'.format(int(x)) for x in ylabels]
				if len(ylabels[-1]) > 3:
					ylabels_temp = []
					if len(ylabels[-1]) > 7:
						for label in ylabels:
							if len(label) > 7:
								ylabels_temp.append(label[0:-8] + "m")
							elif len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)

					else:
						for label in ylabels:
							if len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)
					ylabels = ylabels_temp


			labs = np.arange(0.44,78.44,1)
			panel1.set_xlim([0, 78])
			panel1.set_ylim([0, y])
			panel1.set_xticks(labs)
			panel1.set_yticks(ylabs)
			panel1.set_xticklabels(xlabels, rotation='vertical', fontsize=30, color='grey', fontname='Courier New', verticalalignment='top', fontweight='bold')

			panel1.set_yticklabels(ylabels, fontsize=25)
			plt.gca().yaxis.grid(True)
			plt.gca().grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
			panel1.set_xlabel('')
			panel1.set_ylabel('')

			if percentage:
				plt.ylabel("Percentage of Double Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')
			else:
				plt.ylabel("Number of Double Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')

			panel1.tick_params(axis='both',which='both',\
							   bottom=False, labelbottom=True,\
							   left=True, labelleft=True,\
							   right=True, labelright=False,\
							   top=False, labeltop=False,\
							   direction='in', length=25, colors='lightgray', width=2)

			[i.set_color("black") for i in plt.gca().get_yticklabels()]
			[i.set_color("grey") for i in plt.gca().get_xticklabels()]

			buffer = io.BytesIO()
			plt.savefig(buffer, format="png")
			buff_list[sample]=buffer			
			plt.close()
			sample_count += 1
		return buff_list
		#pp.close()
	# except:
	# 	print("There may be an issue with the formatting of you matrix file.")
	# 	os.remove(output_path + 'DBS_78_plots_' + project + '.pdf')

	elif plot_type == '312' or plot_type == '78SB' or plot_type == 'SB78' or plot_type == '186':
		with open(matrix_path) as f:
			next(f)
			first_line = f.readline()
			first_line = first_line.strip().split()
			mutation_type = first_line[0]
			if len(mutation_type) != 7 and mutation_type[1] != ':':
				sys.exit("The matrix does not match the correct SBS96 format. Please check you formatting and rerun this plotting function.")

		#pp = PdfPages(output_path + 'DBS_186_plots_' + project + '.pdf')

		dinucs = ['TT>GG','TT>CG','TT>AG','TT>GC','TT>CC','TT>AC','TT>GA','TT>CA','TT>AA',
				  'CT>AA','CT>AC','CT>AG','CT>GA','CT>GC','CT>GG','CT>TG','CT>TC','CT>TA',
				  'CC>AA','CC>AG','CC>AT','CC>GA','CC>GG','CC>GT','CC>TA','CC>TG','CC>TT',
				  'TC>GT','TC>CT','TC>AT','TC>GG','TC>CG','TC>AG','TC>GA','TC>CA','TC>AA']

		revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','>':'>'}[B] for B in x][::-1])
		mutations = OrderedDict()

		with open (matrix_path) as f:
			first_line = f.readline()
			samples = first_line.strip().split("\t")
			samples = samples[1:]
			for sample in samples:
				mutations[sample] = OrderedDict()
				mutations[sample]['CC'] = OrderedDict()
				mutations[sample]['CT'] = OrderedDict()
				mutations[sample]['TC'] = OrderedDict()
				mutations[sample]['TT'] = OrderedDict()

			for lines in f:
				line = lines.strip().split()
				mut = line[0][2:]
				nuc = line[0][5:]
				mut_type = line[0][2:4]
				bias = line[0][0]
				if bias == 'N' or bias == 'B' or bias == 'Q':
					continue
				else:
					if mut not in dinucs:
						if revcompl(mut) not in dinucs:
							continue
						nuc = revcompl(nuc)
						mut_type = revcompl(mut_type)
					sample_index = 1

					for sample in samples:
						if percentage:
							mutCount = float(line[sample_index])
							if mutCount < 1 and mutCount > 0:
								sig_probs = True
						else:
							mutCount = int(line[sample_index])
						if nuc not in mutations[sample][mut_type]:
							mutations[sample][mut_type][nuc] = [0,0]
						if bias == 'T':
							mutations[sample][mut_type][nuc][0] = mutCount
						else:
							mutations[sample][mut_type][nuc][1] = mutCount
						sample_index += 1

			for sample in mutations.keys():
				total_count = sum(sum(sum(tsb) for tsb in nuc.values()) for nuc in mutations[sample].values())
				plt.rcParams['axes.linewidth'] = 2
				plot1 = plt.figure(figsize=(21,9.92))
				plt.rc('axes', edgecolor='lightgray')
				panel1 = plt.axes([0.07, 0.09, 0.92, 0.77])
				xlabels = []

				x = 0.3
				ymax = 0
				i = 0
				colors = [[3/256,189/256,239/256], [3/256,102/256,204/256],[162/256,207/256,99/256],
						  [1/256,102/256,1/256], [255/256,153/256,153/256], [228/256,41/256,38/256],
						  [255/256,178/256,102/256], [255/256,128/256,1/256], [204/256,153/256,255/256],
						  [76/256,1/256,153/256]]
				for key in mutations[sample]:
					muts = mutations[sample][key].keys()
					muts = sorted(muts)
					for seq in muts:
						xlabels.append(seq)
						if percentage:
							try:
								trans = plt.bar(x, mutations[sample][key][seq][0]/total_count*100,width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
								x += 0.2
								untrans = plt.bar(x, mutations[sample][key][seq][1]/total_count*100,width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
								x += .8
								if mutations[sample][key][seq][0]/total_count*100 > ymax:
										ymax = mutations[sample][key][seq][0]/total_count*100
								if mutations[sample][key][seq][1]/total_count*100 > ymax:
										ymax = mutations[sample][key][seq][1]/total_count*100
							except:
								trans = plt.bar(x, 0,width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
								untrans = plt.bar(x, 0, width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')

						else:
							trans = plt.bar(x, mutations[sample][key][seq][0],width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
							x += 0.2
							untrans = plt.bar(x, mutations[sample][key][seq][1],width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
							x += .8
							if mutations[sample][key][seq][0] > ymax:
									ymax = mutations[sample][key][seq][0]
							if mutations[sample][key][seq][1] > ymax:
									ymax = mutations[sample][key][seq][1]
					i += 1


				y3 = .87
				y = int(ymax*1.25)

				panel1.add_patch(plt.Rectangle((.075,y3), .218, .05, facecolor=colors[0], clip_on=False, transform=plt.gcf().transFigure))
				panel1.add_patch(plt.Rectangle((.302,y3), .218, .05, facecolor=colors[2], clip_on=False, transform=plt.gcf().transFigure))
				panel1.add_patch(plt.Rectangle((.532,y3), .218, .05, facecolor=colors[4], clip_on=False, transform=plt.gcf().transFigure))
				panel1.add_patch(plt.Rectangle((.765,y3), .218, .05, facecolor=colors[7], clip_on=False, transform=plt.gcf().transFigure))

				yText = y3 + .06
				plt.text(.13, yText, 'CC>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
				plt.text(.37, yText, 'CT>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
				plt.text(.59, yText, 'TC>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
				plt.text(.83, yText, 'TT>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

				if y <= 4:
					y += 4

				while y%4 != 0:
					y += 1
				ytick_offest = int(y/4)

				x_shaded = 0
				panel1.add_patch(plt.Rectangle((x_shaded,0), 8.9, y, facecolor=colors[0], zorder=0, alpha = 0.25, edgecolor='grey'))
				x_shaded += 8.9
				panel1.add_patch(plt.Rectangle((x_shaded,0), 9, y, facecolor=colors[2], zorder=0, alpha = 0.25, edgecolor='grey'))
				x_shaded += 9
				panel1.add_patch(plt.Rectangle((x_shaded,0), 9, y, facecolor=colors[4], zorder=0, alpha = 0.25, edgecolor='grey'))
				x_shaded += 9
				panel1.add_patch(plt.Rectangle((x_shaded,0), 9.1, y, facecolor=colors[7], zorder=0, alpha = 0.25, edgecolor='grey'))


				if percentage:
					ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
					ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%",
							  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
				else:

					if ytick_offest == 0:
						ytick_offest = 1
					ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
					ylabels= [0, ytick_offest, ytick_offest*2,
							  ytick_offest*3, ytick_offest*4]


				if sig_probs:
					plt.text(0.08, 0.8, sample, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
				else:
					plt.text(0.08, 0.8, sample + ": " + "{:,}".format(int(total_count)) + " transcribed double subs", fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

				if not percentage:
					ylabels = ['{:,}'.format(int(x)) for x in ylabels]
					if len(ylabels[-1]) > 3:
						ylabels_temp = []
						if len(ylabels[-1]) > 7:
							for label in ylabels:
								if len(label) > 7:
									ylabels_temp.append(label[0:-8] + "m")
								elif len(label) > 3:
									ylabels_temp.append(label[0:-4] + "k")
								else:
									ylabels_temp.append(label)

						else:
							for label in ylabels:
								if len(label) > 3:
									ylabels_temp.append(label[0:-4] + "k")
								else:
									ylabels_temp.append(label)
						ylabels = ylabels_temp


				labs = np.arange(0.55,36.44,1)
				panel1.set_xlim([0, 36])
				panel1.set_ylim([0, y])
				panel1.set_xticks(labs)
				panel1.set_yticks(ylabs)
				panel1.set_xticklabels(xlabels, rotation='vertical', fontsize=30, color='grey', fontname='Courier New', verticalalignment='top', fontweight='bold')

				panel1.set_yticklabels(ylabels, fontsize=25)
				plt.gca().yaxis.grid(True)
				plt.gca().grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
				panel1.set_xlabel('')
				panel1.set_ylabel('')
				panel1.legend(handles=[trans, untrans], prop={'size':30})

				if percentage:
					plt.ylabel("Percentage of Double Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')
				else:
					plt.ylabel("Number of Double Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')

				panel1.tick_params(axis='both',which='both',\
								   bottom=False, labelbottom=True,\
								   left=True, labelleft=True,\
								   right=True, labelright=False,\
								   top=False, labeltop=False,\
								   direction='in', length=25, colors='lightgray', width=2)

				[i.set_color("black") for i in plt.gca().get_yticklabels()]
				[i.set_color("grey") for i in plt.gca().get_xticklabels()]

				panel1.set_xlim([0, 36])
				buffer = io.BytesIO()
				plt.savefig(buffer, format="png")
				buff_list[sample]=buffer
				plt.close()
		return buff_list	
			#pp.close()
		# except:
		# 	print("There may be an issue with the formatting of you matrix file.")
		# 	os.remove(output_path + 'DBS_186_plots_' + project + '.pdf')

	else:
		print("The provided plot_type: ", plot_type, " is not supported by this plotting function")

# def main():
# # # 	#plotSBS("/Users/ebergstr/Desktop/AID/output/SBS/AID.SBS384.all", "/Users/ebergstr/Desktop/", "AID", '384', False, custom_text_upper=['Similarity to PCAWG: 0.98', 'Similarity to PCAWG: 0.9', 'Similarity to PCAWG: 0.7'], custom_text_bottom=['Stability: 0.98', 'Similarity to PCAWG: 0.9', 'Similarity to PCAWG: 0.7', 'Similarity to PCAWG: 0.9'])
# 	plotID("/Users/ebergstr/Desktop/BRCA/output/INDEL/BRCA.INDEL83.all", "/Users/ebergstr/Desktop/", "BRCA", '83', True, custom_text_upper=['Similarity to PCAWG: 0.98', 'Similarity to PCAWG: 0.9', 'Similarity to PCAWG: 0.7'], custom_text_middle=['Stability: 0.98', 'Similarity to PCAWG: 0.9', 'Similarity to PCAWG: 0.7', 'Similarity to PCAWG: 0.9'], custom_text_bottom=['Stability: 0.98', 'Similarity to PCAWG: 0.9', 'Similarity to PCAWG: 0.7', 'Similarity to PCAWG: 0.9'])
# # # 	plotDBS("/Users/ebergstr/Desktop/Mel/output/DINUC/Mel.DBS186.all", "/Users/ebergstr/Desktop/", "Mel", '186', False, custom_text_upper=['Similarity to PCAWG: 0.98', 'Similarity to PCAWG: 0.9', 'Similarity to PCAWG: 0.7'], custom_text_middle=['Stability: 0.98', 'Similarity to PCAWG: 0.9', 'Similarity to PCAWG: 0.7', 'Similarity to PCAWG: 0.9'], custom_text_bottom=['Stability: 0.98', 'Similarity to PCAWG: 0.9', 'Similarity to PCAWG: 0.7', 'Similarity to PCAWG: 0.9'])

# # # # 	#plotSBS("/Users/ebergstr/Desktop/BRCA/output/SBS/BRCA.SBS1536.all", "/Users/ebergstr/Desktop/", "BRCA", '1536', False, custom_text_upper=['Similarity to PCAWG: 0.98'], custom_text_middle= ['Similarity to PCAWG: 0.98'], custom_text_bottom=['Stability: 0.98'])
# # # 	#plotSBS("/Users/ebergstr/Desktop/BRCA/output/SBS/BRCA.SBS1536.all", "/Users/ebergstr/Desktop/", "BRCA", '1536', False, custom_text_upper=['Similarity to PCAWG: 0.98'], custom_text_middle= ['Similarity to PCAWG: 0.98'])

# # # 	# plotDBS("/Users/ebergstr/Downloads/Biliary-AdenoCA.dinucs.csv", "/Users/ebergstr/Desktop/", "test", '78', False)

# if __name__ == '__main__':
# 	main()
