from setuptools import setup, find_packages

long_desc = """A  pipeline to rapidly reconstruct the SSU rRNAs and explore
phylogenetic composition of an Illumina (meta)genomic or transcriptomic
dataset
"""

setup(
    name="phyloflash",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A pipeline to rapidly reconstruct the SSU rRNAs",
    long_description=long_desc,
    long_description_content_type="text/markdown",
    url="https://github.com/nick-youngblut/phyloflash",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires='>=3.7',
    #install_requires=[   ],
    package_data={
        "phyloflash": [
            "data/barrnap-HGV/all/*.hmm", 
            "data/barrnap-HGV/lsu/*.hmm", 
            "data/barrnap-HGV/ssu/*.hmm",
            "data/barrnap-HGV/nhmmer-darwin",
            "data/barrnap-HGV/nhmmer-linux"
        ]
    },
    entry_points={
        "console_scripts": [
            "phyloflash=phyloflash.cli:main"
        ],
    },
)
