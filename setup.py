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
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'scipy==1.3.2',
        'scikit-learn==0.21.3',
        'joblib==0.14',
        'numpy==1.17'],
    package_data={'plasclass': ['data/*']},
    include_package_data=True,
)
