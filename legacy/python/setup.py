from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="nilhmm",
    version="0.1.0",
    author="Fausto Rodriguez-Zapata",
    author_email="faustovrz@gmail.com",
    description="HMM-based introgression calling for maize NIL populations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sawers-rellan-labs/nilhmm",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "scikit-learn>=1.0.0",
        "hmmlearn>=0.2.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.900",
        ],
        "genomics": [
            "pysam>=0.19.0",
            "cyvcf2>=0.30.0",
        ],
        "docs": [
            "sphinx>=4.0",
            "sphinx_rtd_theme>=1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "call-bzea-introgressions=scripts.call_bzea_introgressions:main",
            "nilhmm-optimize=scripts.parameter_tuning:main",
            "nilhmm-preprocess=scripts.preprocess_vcf:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
