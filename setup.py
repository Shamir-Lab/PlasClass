import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="plasclass-dpellow",
    version="0.1",
    author="David Pellow",
    author_email="dpellow@mail.tau.ac.il",
    description="Classification of plasmid sequences",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dpellow/plasclass",
    packages=['plasclass'],
    scripts=['classify_fasta.py'],
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    requires=['python (<3.0)'],
    install_requires=[
        'scipy<1.3',
        'scikit-learn<=0.20',
        'joblib',
        'numpy<1.17'],
    package_data={'plasclass': ['data/*']},    
    include_package_data=True,
)
