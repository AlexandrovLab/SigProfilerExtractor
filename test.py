#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:45:38 2019

@author: mishugeb
"""
from SigProfilerExtractor import sigpro as sig
import SigProfilerExtractor as spe_mod
import os

def run_text():
    data = sig.importdata("text")
    sig.sigProfilerExtractor("text", "test_text_output", data,
            minimum_signatures=3, maximum_signatures=3, nmf_replicates=5,
            min_nmf_iterations=100, max_nmf_iterations=1000, nmf_test_conv=100)

def run_vcf():
    vcf_data = os.path.join(spe_mod.__path__[0], "data/VCFInput/")
    sig.sigProfilerExtractor("vcf", "test_vcf_output", vcf_data,
            minimum_signatures=3, maximum_signatures=3, nmf_replicates=5,
            min_nmf_iterations=100, max_nmf_iterations=1000, nmf_test_conv=100)

def run_cnv():
    data = sig.importdata("seg:BATTENBERG")
    sig.sigProfilerExtractor("seg:BATTENBERG", "test_segCNV_output", data,
            minimum_signatures=3, maximum_signatures=3, nmf_replicates=5,
            min_nmf_iterations=100, max_nmf_iterations=1000, nmf_test_conv=100)

def run_matobj():
    data = sig.importdata("matobj")
    sig.sigProfilerExtractor("matobj", "test_matobj_output", data,
            minimum_signatures=3, maximum_signatures=3, nmf_replicates=5,
            min_nmf_iterations=100, max_nmf_iterations=1000, nmf_test_conv=100)

def run_csv():
    data = sig.importdata("csv")
    sig.sigProfilerExtractor("csv", "test_csv_output", data,
            minimum_signatures=3, maximum_signatures=3, nmf_replicates=5,
            min_nmf_iterations=100, max_nmf_iterations=1000, nmf_test_conv=100)


if __name__ == '__main__':
    run_text()
    # run_vcf()
    run_cnv()
    # run_matobj()
    # run_csv()