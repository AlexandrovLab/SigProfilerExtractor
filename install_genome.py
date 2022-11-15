#!/usr/bin/env python3
import sys
from SigProfilerMatrixGenerator import install as genInstall

def install_ref(ref_path):
  genInstall.install('GRCh37', offline_files_path=ref_path)

if __name__=="__main__":
  ref_path=sys.argv[1]
  install_ref(ref_path)