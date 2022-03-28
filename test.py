#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:45:38 2019

@author: mishugeb
"""
from SigProfilerExtractor import sigpro as sig
def main():
    data = sig.importdata("text")
    sig.sigProfilerExtractor("text", "example_output_modified", data, minimum_signatures=1, maximum_signatures=3, nmf_replicates=3)

if __name__ == '__main__':
    main()
