from setuptools import setup
import shutil
import os

#remove the dist folder first if exists
if os.path.exists("dist"):
    shutil.rmtree("dist")

VERSION = '0.0.5.50'


with open('README.md') as f:
	long_description = f.read()

def write_version_py(filename='sigproextractor/version.py'):
    # Copied from numpy setup.py
    cnt = """
# THIS FILE IS GENERATED FROM SIGPROEXTRACTOR SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
    
    """
    fh = open(filename, 'w')
    fh.write(cnt % {'version': VERSION,})
    fh.close()
    
write_version_py()
setup(name='sigproextractor',
      version=VERSION,
      description='Extracts mutational signatures from mutational catalogues',
      long_description=long_description,
      long_description_content_type='text/markdown',  # This is important!	
      url="https://github.com/AlexandrovLab/SigProfilerExtractor.git",
      author='S Mishu Ashiqul Islam',
      author_email='m0islam@ucsd.edu',
      license='UCSD',
      packages=['sigproextractor'],
      install_requires=[
          'matplotlib>=2.2.2',
          'scipy>=1.1.0', 
          'numpy>=1.14.4', 
          'pandas>=0.23.4', 
          'nimfa>=1.1.0', 
          'SigProfilerMatrixGenerator>=1.0.1', 
          'sigProfilerPlotting>=1.0.1', 
          'pillow',
          'statsmodels>=0.9.0',
          'scikit-learn>=0.20.2'
           ],
      include_package_data=True,      
      zip_safe=False)
