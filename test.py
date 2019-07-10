#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:45:38 2019

@author: mishugeb
"""

from sigproextractor import sigpro as sig
data = sig.importdata("vcf")
sig.sigProfilerExtractor("vcf", "example_output", data, startProcess=1, endProcess=2, totalIterations=3)