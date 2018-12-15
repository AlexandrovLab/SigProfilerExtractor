from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='sigproextractor',
      version='0.5',
      description='Extracts mutational signatures from mutational catalogues',
      url="https://github.com/AlexandrovLab/SigProfilerExtractor.git",
      author='S Mishu Ashiqul Islam',
      author_email='m0islam@ucsd.edu',
      license='UCSD',
      packages=['sigproextractor'],
      include_package_data=True,      
      zip_safe=False)
