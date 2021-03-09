#!/usr/bin/env python3

"""
Created: March 10th, 2020
@author: Mark Barnes

PlotDecomposition works with matrix formats SigProfiler SBS-96, SBS-1536, DBS-78,
and ID-83. This program is intended to take two matrices.

(1) Sample matrix - A SigProfiler formatted SBS-96, SBS-1536, DBS-78, or ID-83
matrix.
(2) Basis matrix - A SigProfiler formatted SBS-96, SBS-1536, DBS-78, or ID-83
matrix that is the decomposition of (1).

When running the function 'run_PlotDecomposition' a plot of the decomposition will
be generated and saved to the output folder. Refer to the function below to learn
more about the parameters required to generate the decomposition plot.
"""

import os
import pandas as pd
import numpy as np
import scipy.stats
import sigProfilerPlotting as pltCNV
from SigProfilerExtractor import SigProfilerPlottingMatrix as sigPlt
from SigProfilerExtractor import PlotDecomposition_SBS96 as spd_96
from SigProfilerExtractor import PlotDecomposition_SBS288 as spd_288
from SigProfilerExtractor import PlotDecomposition_SBS1536 as spd_1536
from SigProfilerExtractor import PlotDecomposition_DBS78 as spd_78
from SigProfilerExtractor import PlotDecomposition_ID83 as spd_83
from SigProfilerExtractor import PlotDecomposition_CNV48 as cnv_48
from SigProfilerExtractor import subroutines as sub
# imports for working with plots in memory
import io
from PIL import Image
from reportlab.lib.utils import ImageReader
# Global Variables
SBS_CONTEXTS = ["6", "24", "96", "288", "384", "1536", "6144"]
DBS_CONTEXTS = ["78", "186", "1248", "2976"]
ID_CONTEXTS = ["28", "83", "415"]
CNV_CONTEXTS = ["48"]
mtype_options = ["6", "24", "96", "384", "1536", "6144", "28", "83", "415", "78", "186", "1248", "2976"]

# Helper function for converting BytesIO to image so it can be plotted by reportlab
def bytes_to_img(byte_png):
	byte_png.seek(0)
	tmp_im=Image.open(byte_png)
	image = ImageReader(tmp_im)
	return image
	
# Helper function to convert byte array to image array
def open_byte_to_img_dict(byte_dict):
	img_dict = dict()
	for name in byte_dict.keys():
		tmp_img = bytes_to_img(byte_dict[name])
		img_dict[name] = tmp_img
	return img_dict
		
	

def calculate_similarities(denovo, denovo_name, est_denovo):
	from numpy import inf

	# If matrix is 1536 context, then collapse it to 96 format
	if denovo.shape[0]==1536:
		index = denovo.iloc[:,0]
		denovo_tmp = pd.DataFrame(denovo, index=index)
		denovo_tmp = denovo.groupby(denovo_tmp.index.str[1:8]).sum()
		denovo = pd.DataFrame(denovo_tmp)
		denovo = denovo.reset_index()
	elif denovo.shape[0]==288:
		index = denovo.iloc[:,0]
		denovo_tmp = pd.DataFrame(denovo, index=index)
		denovo_tmp = denovo.groupby(denovo_tmp.index.str[2:9]).sum()
		denovo = pd.DataFrame(denovo_tmp)
		denovo = denovo.reset_index()


	sample_names = [denovo_name]

	if  sample_names is False:
		sample_names = ["None"]*denovo.shape[1]

	cosine_similarity_list = []
	cosine_distance_list = []
	correlation_list = []
	correlation_distance_list = []
	kl_divergence_list = []
	l1_norm_list = []
	l2_norm_list = []
	relative_l1_list = []
	relative_l2_list = []

	p_i = denovo[denovo_name]
	q_i = est_denovo
	
	cosine_similarity_list.append(round(sub.cos_sim(p_i,q_i ),3))
	cosine_distance_list.append(round(scipy.spatial.distance.cosine(p_i, q_i),3))
	correlation_list.append(round(scipy.stats.pearsonr(p_i,q_i)[0],3))
	correlation_distance_list.append(round(1-scipy.stats.pearsonr(p_i,q_i)[0],3))
	kl_divergence_list.append(round(scipy.stats.entropy(p_i,q_i),4))
	l1_norm_list.append(round(np.linalg.norm(p_i-q_i , ord=1),2))
	relative_l1_list.append(round((l1_norm_list[-1]/np.linalg.norm(p_i, ord=1))*100,3))
	l2_norm_list.append(round(np.linalg.norm(p_i-q_i , ord=2),2))
	relative_l2_list.append(round((l2_norm_list[-1]/np.linalg.norm(p_i, ord=2))*100,3))
	kl_divergence_list = np.array(kl_divergence_list)
	kl_divergence_list[kl_divergence_list == inf] =1000
	
	similarities_dataframe = pd.DataFrame({"Sample Names": sample_names, \
										   "Cosine Similarity": cosine_similarity_list, \
										   "Cosine Distance": cosine_distance_list, \
										   "Correlation Distance": correlation_distance_list, \
										   "Correlation Coefficient": correlation_list, \
										   "L1 Norm": l1_norm_list, \
										   "L1 Norm %":relative_l1_list, \
										   "L2 Norm": l2_norm_list, \
										   "L2 Norm %": relative_l2_list, \
										   "KL Divergence": kl_divergence_list})
	similarities_dataframe = similarities_dataframe.set_index("Sample Names")

	return similarities_dataframe

# Determine if the matrix matches the format indicated by mtype
def matrix_is_formatted(mtx, mtype):

	# Files to be checked against
	FormatFiles_path = "src/FormatFiles/"
	example_SBS6 = FormatFiles_path + "Sample_Files.SBS6.all"
	example_SBS24 = FormatFiles_path + "Sample_Files.SBS24.all"
	example_SBS96 = FormatFiles_path + "Sample_Files.SBS96.all"
	example_SBS384 = FormatFiles_path + "Sample_Files.SBS384.all"
	example_SBS1536 = FormatFiles_path + "Sample_Files.SBS1536.all"
	example_SBS6144 = FormastFiles_path + "Sample_Files.SBS6144.all"
	example_ID28 = FormatFiles_path + "Sample_Files.ID28.all"
	example_ID83 = FormatFiles_path + "Sample_Files.ID83.all"
	example_ID415 = FormatFiles_path + "Sample_Files.ID415.all"
	example_DBS78 = FormatFiles_path + "Sample_Files.DBS78.all"
	example_DBS186 = FormatFiles_path + "Sample_Files.DBS186.all"
	example_DBS1248 = FormatFiles_path + "Sample_Files.DBS1248.all"
	example_DBS2976 = FormatFiles_path + "Sample_Files.DBS2976.all"

	example_files = [example_SBS6, example_SBS24, example_SBS96, example_SBS384, \
	example_SBS1536, example_SBS6144, example_ID28, example_ID83, example_ID415, \
	example_DBS78, example_DBS186, example_DBS1248, example_DBS2976]

	if mtype not in mtype_options:
		raise Exception('Input context format does not match any of the' +
						' existing supported contexts. The full list is: '
						+ str(mtype_options))
		return False

	# check that the input matrix has the correct format
	f1_names = mtx.iloc[:,0]

	file_index = mtype_options.index(mtype)
	ref_file = pd.read_csv(example_files[file_index], sep="\t")
	ref_names = ref_file.iloc[:,0]

	if (f1_names.equals(ref_names)):
		return True
	else:
		return False


def genSBS_pngs(denovo_mtx, basis_mtx, output_path, project, mtype):
	denovo_plots = dict()
	basis_plots = dict()
	if mtype == "1536" or mtype == "288":
		denovo_plots = sigPlt.plotSBS(denovo_mtx, output_path, project, mtype, True)
		basis_plots = sigPlt.plotSBS(basis_mtx, output_path, project, "96", True)
	elif mtype == "96":
		denovo_plots = sigPlt.plotSBS(denovo_mtx, output_path, project, mtype, True)
		basis_plots = sigPlt.plotSBS(basis_mtx, output_path, project, mtype, True)
	return denovo_plots,basis_plots

def genDBS_pngs(denovo_mtx, basis_mtx, output_path, project, mtype):
	denovo_plots = dict()
	basis_plots = dict()
	denovo_plots = sigPlt.plotDBS(denovo_mtx, output_path, project, mtype, True)
	basis_plots = sigPlt.plotDBS(basis_mtx, output_path, project, mtype, True)
	return denovo_plots,basis_plots

def genID_pngs(denovo_mtx, basis_mtx, output_path, project, mtype):
	denovo_plots = dict()
	basis_plots = dict()
	denovo_plots = sigPlt.plotID(denovo_mtx, output_path, project, mtype, True)
	basis_plots = sigPlt.plotID(basis_mtx, output_path, project, mtype, True)
	return denovo_plots,basis_plots
	
def genCNV_pngs(denovo_mtx, basis_mtx, output_path, project, mtype):
	denovo_plots = dict()
	basis_plots = dict()
	denovo_plots = pltCNV.plotCNV(denovo_mtx, output_path, project, plot_type="pdf", percentage=True, aggregate=False, read_from_file=False, write_to_file=False)
	basis_plots = pltCNV.plotCNV(basis_mtx, output_path, project, plot_type="pdf", percentage=True, aggregate=False, read_from_file=False, write_to_file=False)
	return denovo_plots,basis_plots

# signames, weights
def gen_sub_plots(denovo_mtx, basis_mtx, output_path, project, mtype):

	if mtype in SBS_CONTEXTS:
		if not os.path.exists(output_path):
			os.makedirs(output_path)
		denovo_plots,basis_plots=genSBS_pngs(denovo_mtx, basis_mtx, output_path, project, mtype)
		return denovo_plots,basis_plots
	elif mtype in DBS_CONTEXTS:
		if not os.path.exists(output_path):
			os.makedirs(output_path)
		denovo_plots,basis_plots=genDBS_pngs(denovo_mtx, basis_mtx, output_path, project, mtype)
		return denovo_plots,basis_plots
	elif mtype in ID_CONTEXTS:
		if not os.path.exists(output_path):
			os.makedirs(output_path)
		denovo_plots,basis_plots=genID_pngs(denovo_mtx, basis_mtx, output_path, project, mtype)
		return denovo_plots,basis_plots
	elif mtype in CNV_CONTEXTS:
		if not os.path.exists(output_path):
			os.makedirs(output_path)
		denovo_plots, basis_plots=genCNV_pngs(denovo_mtx, basis_mtx, output_path, project, mtype)
		return denovo_plots,basis_plots
	else:
		print("ERROR: mtype is " + mtype + " and is not yet supported.")

# generate the plot for the reconstruction
def gen_reconstructed_png(denovo_name, basis_mtx, basis_names, weights, output_path, project, mtype):
	reconstruction_plot=dict()
	col_names=[denovo_name]
	mut_col = basis_mtx.iloc[:,0]

	recon_plot = basis_mtx[basis_names[0]]*float(weights[0].strip("%"))/100

	for i in range(1,len(weights)):
		recon_plot = recon_plot + basis_mtx[basis_names[i]]*(float(weights[i].strip("%"))/100)

	recon_plot = pd.Series(recon_plot, name=denovo_name)
	reconstruction_mtx = pd.concat([mut_col, recon_plot], axis=1)
	if mtype in SBS_CONTEXTS:
		if mtype == "1536" or mtype == "288":
			reconstruction_plot=sigPlt.plotSBS(reconstruction_mtx, output_path, "reconstruction_" + project, "96", True)
		else:
			reconstruction_plot=sigPlt.plotSBS(reconstruction_mtx, output_path, "reconstruction_" + project, mtype, True)
	elif mtype in DBS_CONTEXTS:
		reconstruction_plot=sigPlt.plotDBS(reconstruction_mtx, output_path, "reconstruction_" + project, mtype, True)
	elif mtype in ID_CONTEXTS:
		reconstruction_plot=sigPlt.plotID(reconstruction_mtx, output_path, "reconstruction_" + project, mtype, True)
	elif mtype in CNV_CONTEXTS:
		 reconstruction_plot = pltCNV.plotCNV(reconstruction_mtx, output_path, "reconstruction_"+project, plot_type="pdf", \
		 	percentage=True, aggregate=False, read_from_file=False, write_to_file=False)
	else:
		print("ERROR: mtype is " + mtype + " and is not yet supported.")

	return reconstruction_mtx,reconstruction_plot


def gen_decomposition(denovo_name, basis_names, weights, output_path, project, \
	mtype, denovo_plots_dict, basis_plots_dict, reconstruction_plot_dict, \
	reconstruction=False, statistics=None, sig_version=None, custom_text=None):

	"""
	Generate the correct plot based on mtype.

	Parameters:
	----------
	denovo_name: 				(String) 			Name of denovo signature
	basis_names: 				(List of Strings) 	Names of basis signatures
	weights:					(List of Strings) 	Percentile contribution for each basis signature
	output_path: 				(String) 			Path to existing output directory
	project: 					(String) 			Project name appended to file names
	mtype: 						(String) 			The context 'mtype_options' has valid values
	denovo_plots_dict			(Dictionary)		Signatures are keys, ByteIO plots are values
	basis_plots_dict			(Dictionary)		Signatures are keys, ByteIO plots are values
	reconstruction_plot_dict	(Dictionary)		Signatures are keys, ByteIO plots are values
	reconstruction: 			(Boolean) 			True to generate plot w/ reconstruction
	statistics: 				(Pandas Dataframe) 	Output from calculate_similarities()
	"""

	if mtype == "6":
		print("Need to add support for SBS6 Decomposition")
	elif mtype == "24":
		print("Need to add support for SBS24 Decomposition")
	elif mtype == "96":
		byte_plot=spd_96.gen_decomposition(denovo_name, basis_names, weights, output_path, \
			project, denovo_plots_dict, basis_plots_dict, reconstruction_plot_dict, \
			reconstruction, statistics, sig_version, custom_text)
		return byte_plot
	elif mtype == "288":
		byte_plot=spd_288.gen_decomposition(denovo_name, basis_names, weights, output_path, \
			project, denovo_plots_dict, basis_plots_dict, reconstruction_plot_dict, \
			reconstruction, statistics, sig_version, custom_text)
		return byte_plot
	elif mtype == "384":
		print("Need to add support for SBS24 Decomposition")
	elif mtype == "1536":
		byte_plot=spd_1536.gen_decomposition(denovo_name, basis_names, weights, output_path, \
			project, denovo_plots_dict, basis_plots_dict, reconstruction_plot_dict, \
			reconstruction, statistics, sig_version, custom_text)
		return byte_plot
	elif mtype == "6144":
		print("Need to add support for SBS6144 Decomposition")
	elif mtype == "28":
		print("Need to add support for ID28 Decomposition")
	elif mtype == "83":
		byte_plot=spd_83.gen_decomposition(denovo_name, basis_names, weights, output_path, \
			project, denovo_plots_dict, basis_plots_dict, reconstruction_plot_dict, \
			reconstruction, statistics, sig_version, custom_text)
		return byte_plot
	elif mtype == "415":
		print("Need to add support for ID415 Decomposition")
	elif mtype == "78":
		byte_plot=spd_78.gen_decomposition(denovo_name, basis_names, weights, output_path, \
			project, denovo_plots_dict, basis_plots_dict, reconstruction_plot_dict, \
			reconstruction, statistics, sig_version, custom_text)
		return byte_plot
	elif mtype == "186":
		print("Need to add support for DBS186 Decomposition")
	elif mtype == "1248":
		print("Need to add support for DBS1248 Decomposition")
	elif mtype == "2976":
		print("Need to add support for DBS2976 Decomposition")
	elif mtype == "48":
		byte_plot=cnv_48.gen_decomposition(denovo_name, basis_names, weights, output_path, \
			project, denovo_plots_dict, basis_plots_dict, reconstruction_plot_dict, \
			reconstruction, statistics, sig_version, custom_text)
		return byte_plot



def run_PlotDecomposition(denovo_mtx, denovo_name, basis_mtx, basis_names, \
		weights, nonzero_exposures, output_path, project, mtype, sig_version=None, custom_text=None):
	"""
	Generates a decomposition plot of the denovo_mtx using the basis_mtx.

	Parameters:
	----------

	denovo_mtx: Pandas Dataframe. This format represents the catalog of mutations seperated by tab.

	denovo_name: String. The name of the one sample in denovo_mtx to decompose.

	basis_mtx: Pandas Dataframe. This format represents the catalog of mutations seperated by tab.

	basis_names: List of Strings. The names of the samples in denovo_mtx that
	the denovo_name sample from denovo_mtx is decomposed into.
	ie. basis_names=["SBS1", "SBS5", "SBS15", "SBS20"]

	weights: List of Strings. The percentile weight corresponding to each basis
	in basis_names. Refer to example function call below for more detail.
	ie. weights=["11.58%", "42.38%", "16.46%", "29.58%"]

	output_path: String. Path to where to store the output.

	project: String. This string is appended to the file name output.

	mtype: String. The context of the data. Valid values include: "96", "1536","78", and "83".

	sig_version: String. The version of signatures being used.

	custom_text: String. A custom message displayed on decomposition plot.

	Returns:
	-------
	None.
	"""
	
	denovo_plots_dict,basis_plots_dict=gen_sub_plots(denovo_mtx, basis_mtx, output_path, project, mtype)
	reconstructed_mtx,reconstruction_plot_dict = gen_reconstructed_png(denovo_name, basis_mtx, basis_names, \
		weights, output_path, project, mtype)
	
	present_sigs=np.array(basis_mtx[basis_names])
	reconstructed_mtx = np.dot(present_sigs,nonzero_exposures)
	
	denovo_plots_dict = open_byte_to_img_dict(denovo_plots_dict)
	basis_plots_dict = open_byte_to_img_dict(basis_plots_dict)
	reconstruction_plot_dict = open_byte_to_img_dict(reconstruction_plot_dict)
	
	statistics=calculate_similarities(denovo_mtx, denovo_name, reconstructed_mtx)
	byte_plot = gen_decomposition(denovo_name, basis_names, weights, output_path, project, \
		mtype, denovo_plots_dict, basis_plots_dict, reconstruction_plot_dict, \
		reconstruction=True, statistics=statistics, sig_version=sig_version, \
		custom_text=custom_text)
	

	return byte_plot
