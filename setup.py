from pathlib import Path

from setuptools import setup


setup(
    name="phytclust",
    version="1.0.1",
    author="Katyayni Ganesan",
    author_email="katyayni.ganesan@iccb-cologne.org",
    description="Monophyletic Clustering Algorithm for Phylogenetic Trees",
    long_description=(Path(__file__).parent / "README.md").read_text(),
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/schwarzlab/phylotreeclustering",  # To be changed
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=["PhytClust"],
    install_requires=[
        "numpy>=1.20.1",
        "pandas>=1.2.2",
        "biopython>=1.78",
        "matplotlib>=3.3.4",
        "scikit-learn>=1.2.2",
        "scipy>=1.4.0",
    ],
)
