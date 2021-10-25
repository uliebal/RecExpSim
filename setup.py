import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='biolabsim',
    version='0.1.0',
    author='Ulf Liebal',
    author_email='ulf.liebal@rwth-aachen.de',
    description='biolabsim is a collection of workflows using the silvio software.',
    keywords='biotechnology, microbiology, virtual cell, systems biology',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://git.rwth-aachen.de/ulf.liebal/biolabsim.git',
    packages=[], # contains no packages, only notebooks.
    package_dir={'': 'src'},
    package_data = {'': ['*.csv','*.pkl','*.xml']},
    python_requires='~=3.9',
    install_requires=[
        'biopython ~= 1.79',
        'cobra ~= 0.22',
        'joblib ~= 1.0',
        'matplotlib ~= 3.3',
        'numpy ~= 1.18',
        'pandas ~= 1.2',
        'pickleshare ~= 0.7',
        'scipy ~= 1.6',
    ],
)
