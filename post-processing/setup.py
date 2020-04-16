from setuptools import setup, find_packages

project_name = "legopy"
required_packages = ['numpy', 'matplotlib', 'f90nml']
version = '0.1'

setup(
    name = project_name,
    version = version,
    description = 'Post-processing framework for the Legolas code',
    author = 'Niels Claes',
    author_email = 'niels.claes@kuleuven.be',
    keywords = 'interface data-analysis',
    python_requires = '>=3.6',
    install_requires = required_packages,
    packages = find_packages()
)
