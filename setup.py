from pathlib import Path

from setuptools import setup


setup(
    name='phylotreeclus',
    version='1.0.2',
    author='Elisa Billard',
    author_email='elisabillard1905@gmail.com',
    description='Clustering algorithm for phylogenetic cancer evolution trees',
    long_description=(Path(__file__).parent / "README.md").read_text(),
    long_description_content_type='text/markdown',
    url='https://bitbucket.org/schwarzlab/phylotreeclustering',
    classifiers=[
                "Programming Language :: Python :: 3",
                "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent",
    ],
    packages=['PhyloTreeClustering'],
    scripts=['phylotreeclus'],
    install_requires=[
        'numpy',
        'pandas',
        'biopython',
        'matplotlib',
        'scikit-learn'
    ]
)