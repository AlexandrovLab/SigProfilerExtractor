#!/usr/bin/env python3
"""
Created: February 21, 2020
@author: Mark Barnes
"""
import reportlab
import os
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, A4, landscape
from reportlab.lib import utils
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
import SigProfilerExtractor as cosmic
from PyPDF2 import PdfFileWriter, PdfFileReader, PdfFileMerger
# imports for saving plots to memory
import io
from PIL import Image
# imports for dashed line
from reportlab.lib.colors import black
paths = cosmic.__path__[0]
# Page Formatting
inch = 72
# USING LETTER LANDSCAPE DIMENSIONS
WIDTH_LETTER = 11 * inch
HEIGHT_LETTER = 8.5 * inch
MID_WIDTH_LETTER = 396
MID_HEIGHT_LETTER = HEIGHT_LETTER/2

WIDTH_GRAPH = 375
HEIGHT_GRAPH = 85

# Layout Formatting
HEIGHT_GAP = 93
WIDTH_GAP = 6
X_COORD = 0
Y_COORD = 1

# Coordinates for graphs on right side of plot
GRAPH_X_COORD = (WIDTH_LETTER) - WIDTH_GRAPH
GRAPH_Y_COORD = (HEIGHT_LETTER - HEIGHT_GAP)
TEXT_X_COORD = GRAPH_X_COORD + WIDTH_GRAPH - 50
TEXT_Y_COORD = (HEIGHT_LETTER - HEIGHT_GAP) + 63.75

reportlab.rl_config.TTFSearchPath.append(paths+'/src/Fonts/')
pdfmetrics.registerFont(TTFont('Arial-Bold', 'Arial Bold.ttf'))

# Pairs are (x-coordinate, y-coordinate)
LAYOUT_5_GRAPH = [
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 1) + 10),
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 2) + 5),
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 3)),
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 4) - 5),
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 5) - 10)]

LAYOUT_5_TEXT = [
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 1) + 10),
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 2) + 5),
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 3)),
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 4) - 5),
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 5) - 10)]


LAYOUT_4_GRAPH = [
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 1) - HEIGHT_GRAPH/2 + 10),
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 2) - HEIGHT_GRAPH/2 + 5),
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 3) - HEIGHT_GRAPH/2 - 5),
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 4) - HEIGHT_GRAPH/2 - 10)]

LAYOUT_4_TEXT = [
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 1)- HEIGHT_GRAPH/2 + 10),
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 2)- HEIGHT_GRAPH/2 + 5),
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 3)- HEIGHT_GRAPH/2 - 5),
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 4)- HEIGHT_GRAPH/2 - 10)]


LAYOUT_3_GRAPH = [
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 2) + 5),
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 3)),
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 4) - 5)]

LAYOUT_3_TEXT = [
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 2) + 5),
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 3)),
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 4) - 5)]


LAYOUT_2_GRAPH = [
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 2) - HEIGHT_GRAPH/2 + 5),
(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 3) - HEIGHT_GRAPH/2 - 5)]

LAYOUT_2_TEXT = [
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 2)- HEIGHT_GRAPH/2 + 5),
(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 3)- HEIGHT_GRAPH/2 - 5)]


LAYOUT_1_GRAPH = [(GRAPH_X_COORD - WIDTH_GAP, (GRAPH_Y_COORD - HEIGHT_GRAPH * 3))]

LAYOUT_1_TEXT = [(TEXT_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 3))]

# Pairs (position, height), organized from 1 plot to 5 plus plots
BRACKET_SIZES = [
(HEIGHT_LETTER * (3/8), HEIGHT_LETTER * (2/8)),
(HEIGHT_LETTER * (2.5/8), HEIGHT_LETTER * (3/8)),
(HEIGHT_LETTER * (2/8), HEIGHT_LETTER * (4/8)),
(HEIGHT_LETTER * (1.5/8), HEIGHT_LETTER * (5/8)),
(HEIGHT_LETTER * (1/8), HEIGHT_LETTER * (6/8))]

PLOT_NAME = 0
CONTRIBUTION = 1

# Helper functions for plotting the layout of a graph with 1-5 basis signatures
# Parameters:
#   bases 		- (List of Strings) The list of basis names
#   output_path - (String) The path to where the .png files are stored.
#   project 	- (String) The name of the project that is post-fixed to each file name.
#	c_draw 		- (Canvas) The canvas to draw the graph decomposition on.
def plot_1(bases, project, c_draw, denovo_plots_dict, basis_plots_dict):
	for i in range(0,1):
		image=basis_plots_dict[bases[i][0]]
		c_draw.drawImage(image, LAYOUT_1_GRAPH[i][X_COORD], LAYOUT_1_GRAPH[i][Y_COORD], width=WIDTH_GRAPH, height=HEIGHT_GRAPH)
		c_draw.drawString(LAYOUT_1_TEXT[i][X_COORD], LAYOUT_1_TEXT[i][Y_COORD], str(bases[i][1]) + "%")

def plot_2(bases, project, c_draw, denovo_plots_dict, basis_plots_dict):
	for i in range(0,2):
		image=basis_plots_dict[bases[i][0]]
		c_draw.drawImage(image, LAYOUT_2_GRAPH[i][X_COORD], LAYOUT_2_GRAPH[i][Y_COORD], width=WIDTH_GRAPH, height=HEIGHT_GRAPH)
		c_draw.drawString(LAYOUT_2_TEXT[i][X_COORD], LAYOUT_2_TEXT[i][Y_COORD], str(bases[i][1]) + "%")

def plot_3(bases, project, c_draw, denovo_plots_dict, basis_plots_dict):
	for i in range(0,3):
		image=basis_plots_dict[bases[i][0]]
		c_draw.drawImage(image, LAYOUT_3_GRAPH[i][X_COORD], LAYOUT_3_GRAPH[i][Y_COORD], width=WIDTH_GRAPH, height=HEIGHT_GRAPH)
		c_draw.drawString(LAYOUT_3_TEXT[i][X_COORD], LAYOUT_3_TEXT[i][Y_COORD], str(bases[i][1]) + "%")

def plot_4(bases, project, c_draw, denovo_plots_dict, basis_plots_dict):
	for i in range(0,4):
		image=basis_plots_dict[bases[i][0]]
		c_draw.drawImage(image, LAYOUT_4_GRAPH[i][X_COORD], LAYOUT_4_GRAPH[i][Y_COORD], width=WIDTH_GRAPH, height=HEIGHT_GRAPH)
		c_draw.drawString(LAYOUT_4_TEXT[i][X_COORD], LAYOUT_4_TEXT[i][Y_COORD], str(bases[i][1]) + "%")

def plot_5(bases, project, c_draw, denovo_plots_dict, basis_plots_dict):
	for i in range(0,5):
		image=basis_plots_dict[bases[i][0]]
		c_draw.drawImage(image, LAYOUT_5_GRAPH[i][X_COORD], LAYOUT_5_GRAPH[i][Y_COORD], width=WIDTH_GRAPH, height=HEIGHT_GRAPH)
		c_draw.drawString(LAYOUT_5_TEXT[i][X_COORD], LAYOUT_5_TEXT[i][Y_COORD], str(bases[i][1]) + "%")

def plot_6_plus(bases, project, c_draw, denovo_plots_dict, basis_plots_dict):
	for i in range(0,5):
		image=basis_plots_dict[bases[i][0]]
		c_draw.drawImage(image, LAYOUT_5_GRAPH[i][X_COORD], LAYOUT_5_GRAPH[i][Y_COORD], width=WIDTH_GRAPH, height=HEIGHT_GRAPH)
		c_draw.drawString(LAYOUT_5_TEXT[i][X_COORD], LAYOUT_5_TEXT[i][Y_COORD], str(bases[i][1]) + "%")

	extra_sigs = "* "
	for i in range(5, len(bases)-1):
		extra_sigs += str(bases[i][0]) + " (" + str(bases[i][1]) + "%), "

	extra_sigs += bases[len(bases)-1][0] + " (" + str(bases[len(bases)-1][1]) + "%)"
	c_draw.drawString(GRAPH_X_COORD, (TEXT_Y_COORD - HEIGHT_GRAPH * 6) - 10, extra_sigs)

# Helper function to add calculations to layout
# Parameters:
#	c_draw 		- (Canvas) The canvas to draw the graph decomposition on.
#	statistics 	- (Pandas Dataframe) Dataframe w/ calculations
def draw_statistics(c_draw, statistics, sig_version, custom_text):
	
	cos_sim = statistics["Cosine Similarity"][0]
	cor_coeff = statistics["Correlation Coefficient"][0]
	l1_norm_percent = statistics["L1 Norm %"][0]
	l2_norm_percent = statistics["L2 Norm %"][0]
	kl_divergence = statistics["KL Divergence"][0]


	c_draw.drawString(WIDTH_GAP+15, LAYOUT_2_TEXT[1][Y_COORD]-90, \
		"Cosine Similarity: " + str(cos_sim))
	c_draw.drawString(WIDTH_GAP+15, LAYOUT_2_TEXT[1][Y_COORD]-100, \
		"Correlation: " + str(cor_coeff))
	c_draw.drawString(WIDTH_GAP+105, LAYOUT_2_TEXT[1][Y_COORD]-90, \
		"L1 Error %: " + str(l1_norm_percent) + "%")
	c_draw.drawString(WIDTH_GAP+105, LAYOUT_2_TEXT[1][Y_COORD]-100, \
		"L2 Error %: " + str(l2_norm_percent) + "%")
	c_draw.drawString(WIDTH_GAP+195, LAYOUT_2_TEXT[1][Y_COORD]-90, \
		"KL Divergence: " + str(kl_divergence))

	if sig_version is not None:
		c_draw.drawString(WIDTH_GAP+195, LAYOUT_2_TEXT[1][Y_COORD]-100, \
			"Signature Version: " + str(sig_version))
	if custom_text is not None:
		c_draw.drawString(WIDTH_GAP+15, LAYOUT_2_TEXT[1][Y_COORD]-120, \
			str(custom_text))

# Helper function to resize bracket depending on number of bases plotted
# Parameters:
#	num_bases 	- (Integer) The number of bases to be plotted.
#	c_draw 		- (Canvas) The canvas to draw the graph decomposition on.
def draw_bracket(num_bases, c_draw):
	num_plts = num_bases - 1
	if(num_bases >= 5):
		num_plts = 4
	c_draw.drawImage(paths+"/src/Accolade_fermante.png", MID_WIDTH_LETTER - 15, \
		BRACKET_SIZES[num_plts][0], width = 20, height = BRACKET_SIZES[num_plts][1], mask='auto')

def crop_margins(pdf_to_edit, num_bases):
	pdf_to_edit.seek(0)
	pdf_file = PdfFileReader(pdf_to_edit, "rb")
	page = pdf_file.getPage(0)
	writer = PdfFileWriter()
	output_plot_buff = io.BytesIO()
	
	if (num_bases == 1):
		page.mediaBox.lowerRight = (792,155)
		page.mediaBox.lowerLeft = (0,155)
		page.mediaBox.upperRight = (792,402)
		page.mediaBox.upperLeft = (0,402)
		writer.addPage(page)
		writer.write(output_plot_buff)
	elif (num_bases == 2):
		page.mediaBox.lowerRight = (792,155)
		page.mediaBox.lowerLeft = (0,155)
		page.mediaBox.upperRight = (792,422)
		page.mediaBox.upperLeft = (0,422)
		writer.addPage(page)
		writer.write(output_plot_buff)
	elif (num_bases == 3):
		page.mediaBox.lowerRight = (792,150)
		page.mediaBox.lowerLeft = (0,150)
		page.mediaBox.upperRight = (792,462)
		page.mediaBox.upperLeft = (0,462)
		writer.addPage(page)
		writer.write(output_plot_buff)
	elif (num_bases == 4):
		page.mediaBox.lowerRight = (792,112)
		page.mediaBox.lowerLeft = (0,112)
		page.mediaBox.upperRight = (792,498)
		page.mediaBox.upperLeft = (0,498)
		writer.addPage(page)
		writer.write(output_plot_buff)
	elif (num_bases == 5):
		page.mediaBox.lowerRight = (792,75)
		page.mediaBox.lowerLeft = (0,75)
		page.mediaBox.upperRight = (792,537)
		page.mediaBox.upperLeft = (0,537)
		writer.addPage(page)
		writer.write(output_plot_buff)
	elif (num_bases > 5):
		page.mediaBox.lowerRight = (792,50)
		page.mediaBox.lowerLeft = (0,50)
		page.mediaBox.upperRight = (792,537)
		page.mediaBox.upperLeft = (0,537)
		writer.addPage(page)
		writer.write(output_plot_buff)
	return output_plot_buff

# Parameters:
#   de_novo_name 				(String) 			The name of the denovo signature.
#   basis_names					(List of Strings)	The names of the basis signatures
#	output_path 				(String)			Path to where to save the output.
#	project						(String)			The project name that is appended to file names.
#	c							(Canvas)			The canvas that is being drawn on.
#	reconstruction				(Boolean) 			True to create reconstruction
#	statistics					(Pandas Dataframe) 	If reconstructing, then include statistics.
#	sig_version					(String) 			The version of the Cosmic Signatures used
#	denovo_plots_dict			(Dictionary) 		Signatures are keys, ByteIO plots are values
#	basis_plots_dict			(Dictionary) 		Signatures are keys, ByteIO plots are values
#	reconstruction_plot_dict	(Dictionary) 		Signatures are keys, ByteIO plots are values
#
# Output:
#   A graph of the de_novo signature's decomposition.
def gen_plot(de_novo_name, bases, output_path, project, c, reconstruction, \
	statistics, sig_version, custom_text, denovo_plots_dict, basis_plots_dict, \
	reconstruction_plot_dict):

	# THIS IS THE RIGHT SIDE, IT CHANGES
	num_bases = len(bases)
	# print("num bases is: " + str(num_bases))
	if (num_bases == 1):
		plot_1(bases, project, c, denovo_plots_dict, basis_plots_dict)
	elif (num_bases == 2):
		plot_2(bases, project, c, denovo_plots_dict, basis_plots_dict)
	elif (num_bases == 3):
		plot_3(bases, project, c, denovo_plots_dict, basis_plots_dict)
	elif (num_bases == 4):
		plot_4(bases, project, c, denovo_plots_dict, basis_plots_dict)
	elif (num_bases == 5):
		plot_5(bases, project, c, denovo_plots_dict, basis_plots_dict)
	elif (num_bases > 5):
		plot_6_plus(bases, project, c, denovo_plots_dict, basis_plots_dict)

	recon_image=reconstruction_plot_dict[de_novo_name]
	denovo_image=denovo_plots_dict[de_novo_name]
	# THIS IS THE LEFT SIDE
	if reconstruction:
		c.drawImage(denovo_image, WIDTH_GAP, LAYOUT_2_GRAPH[0][Y_COORD], width=WIDTH_GRAPH, height=HEIGHT_GRAPH)
		c.drawImage(recon_image, WIDTH_GAP, LAYOUT_2_GRAPH[1][Y_COORD], width=WIDTH_GRAPH, height=HEIGHT_GRAPH)
		c.drawString(WIDTH_GRAPH - WIDTH_GAP - 25, LAYOUT_2_TEXT[0][Y_COORD], "Original")
		c.drawString(WIDTH_GRAPH - WIDTH_GAP - 45, LAYOUT_2_TEXT[1][Y_COORD], "Reconstructed")
		draw_statistics(c, statistics, sig_version, custom_text)

		# Draw dashed line
		c.setLineWidth(2)
		c.setDash(25,5)
		c.setStrokeColor(black)
		c.setFillColorRGB(255,255,255)
		p = c.beginPath()
		p.moveTo(WIDTH_GAP,MID_HEIGHT_LETTER)
		p.lineTo(WIDTH_GRAPH + 2,MID_HEIGHT_LETTER)
		c.drawPath(p, stroke=1, fill=1)
	else:
		c.drawImage(denovo_image, WIDTH_GAP, MID_HEIGHT_LETTER - HEIGHT_GRAPH/2, width=WIDTH_GRAPH, height=HEIGHT_GRAPH)
	draw_bracket(len(bases), c)

	c.showPage()

# iterate over the csv file and generate plots accordingly
def gen_decomposition(denovo_name, basis_names, weights, output_path, project, \
	denovo_plots_dict, basis_plots_dict, reconstruction_plot_dict, reconstruction, \
	statistics,  sig_version=None, custom_text=None):
	
	buff = io.BytesIO()
	c = canvas.Canvas(buff, pagesize=letter)
	c.setPageSize(landscape(letter))
	c.setFont("Arial-Bold", 7.19)

	basis_plots = []
	for i in range(0,len(basis_names)):
		basis_plots.append([basis_names[i], weights[i]])


	# create for loop to iterate through list, then change second value in list of lists
	# Otherwise sorts strings and then 5.14% > 48.54%
	for j in range(0, len(basis_names)):
		basis_plots[j][1] = float(basis_plots[j][1].strip("%"))
	sorted_list = sorted(basis_plots, key=lambda tup: tup[1], reverse=True)

	gen_plot(denovo_name, sorted_list, output_path, project, c, reconstruction, \
		statistics, sig_version, custom_text, denovo_plots_dict, basis_plots_dict, \
		reconstruction_plot_dict)

	c.save()
	
	# Take the plot and crop the margins
	byte_plot = crop_margins(buff, len(basis_names))
	buff.close()
	return byte_plot
