import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="maginator", 
    version="0.1.8",
    author="Jakob Russel",
    author_email="russel2620@gmail.com",
    description="MAGinator: Abundance, strain, and functional profiling of MAGs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Russel88/MAGinator",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Development Status :: 4 - Beta"],
    python_requires='>=3.5',
    install_requires=[
        "setuptools"],
    entry_points = {
        'console_scripts': ['maginator = maginator.main:cli']
    },
    include_package_data=True,
    zip_safe=False,
    package_data = {'maginator': ['workflow/*']}
)
