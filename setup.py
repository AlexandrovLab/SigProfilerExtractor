from setuptools import setup
import shutil
import os
import sys
import subprocess

#remove the dist folder first if exists
if os.path.exists("dist"):
    shutil.rmtree("dist")

VERSION = '1.1.0'


with open('README.md') as f:
	long_description = f.read()

def write_version_py(filename='SigProfilerExtractor/version.py'):
    # Copied from numpy setup.py
    cnt = """
# THIS FILE IS GENERATED FROM SIGPROFILEREXTRACTOR SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
Update = 'Heirarchy option deleted, clustering deleted and signatures orders by the mutation burden'
    
    """
    fh = open(filename, 'w')
    fh.write(cnt % {'version': VERSION,})
    fh.close()
requirements=[
          'matplotlib>=3.3.0',
          'scipy>=1.3.3',
          'torch==1.5.1',
          'numpy>=1.17.4', 
          'pandas>=0.25.3', 
          'nimfa>=1.1.0', 
          'SigProfilerMatrixGenerator>=1.1.23', 
          'sigProfilerPlotting>=1.1.8', 
          'pillow',
          'statsmodels>=0.9.0',
          'scikit-learn>=0.22',
	  'psutil>=5.6.1',
          'reportlab>=3.5.42',
          'PyPDF2>=1.26.0',
      'xlrd==1.2.0'
           ]

operating_system = sys.platform   
print(operating_system)
if operating_system  in ['win32','cygwin','windows']:
    requirements.remove('matplotlib>=3.3.0')
    requirements.remove('torch==1.5.1')
    print('Trying to install pytorch!')
    code = 1
    try:
        code = subprocess.call(['pip', 'install', 'torch===1.5.1+cpu',  '-f', 'https://download.pytorch.org/whl/torch_stable.html'])
        if code != 0:
            raise Exception('Torch  instalation failed !')
    except:
        try:
            code = subprocess.call(['pip3', 'install', 'torch===1.5.1+cpu',  '-f', 'https://download.pytorch.org/whl/torch_stable.html'])
            if code != 0:
                raise Exception('Torch instalation failed !')
        except:
            print('Failed to install pytroch, please install pytroch manually be following the simple instructions over at: https://pytorch.org/get-started/locally/')
    if code == 0:
        print('Successfully installed pytorch version! (If you need the GPU version, please install it manually, checkout the mindsdb docs and the pytroch docs if you need help)')
    
    
write_version_py()
setup(name='SigProfilerExtractor',
      version=VERSION,
      description='Extracts mutational signatures from mutational catalogues',
      long_description=long_description,
      long_description_content_type='text/markdown',  # This is important!	
      url="https://github.com/AlexandrovLab/SigProfilerExtractor.git",
      author='S Mishu Ashiqul Islam',
      author_email='m0islam@ucsd.edu',
      license='UCSD',
      packages=['SigProfilerExtractor'],
      install_requires=requirements,
      include_package_data=True,      
      zip_safe=False)
