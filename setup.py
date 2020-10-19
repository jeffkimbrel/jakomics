import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="jakomics",
    version="0.10.11",
    author="Jeff Kimbrel",
    author_email="jakpot@gmail.com",
    description="Various omics tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['jakomics'],
    # classifiers=[
    #     "Programming Language :: Python :: 3",
    #     "License :: OSI Approved :: MIT License",
    #     "Operating System :: OS Independent",
    # ],
    python_requires='>=3.6',
)
