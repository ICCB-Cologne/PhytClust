from setuptools import setup

setup(
    name='phylotreeclus',
    version='1.0.0',
    author='Elisa Billard',
    author_email='elisabillard1905@gmail.com',
    description='Clustering algorithm for phylogenetic cancer evolution trees',
    url='https://bitbucket.org/schwarzlab/phylotreeclustering',
    classifiers=[
                "Programming Language :: Python :: 3",
                "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent",
    ],
    packages=['PhyloTreeClustering'],
    scripts=['phylotreeclus'],
    install_requires=[
        'numpy>=1.25.0',
        'pandas>=1.5.3',
        'biopython>=1.78',
        'matplotlib>=3.7.1',
        'scikit-learn>=1.2.2'
    ]
)