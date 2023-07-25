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
    url="https://github.com/Shamir-Lab/PlasClass",
    packages=['plasclass'],
    scripts=['classify_fasta.py', 'plasclass/train.py'],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'scipy==1.10.0',
        'scikit-learn==0.21.3',
        'joblib==1.2.0',
        'numpy==1.22.0'],
    package_data={'plasclass': ['data/*']},
    include_package_data=True,
)
