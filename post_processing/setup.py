from setuptools import setup, find_packages
from pathlib import Path

package_name = "pylbo"
required_packages = ['numpy', 'matplotlib', 'f90nml', 'tqdm', 'psutil']

version_filepath = (Path(__file__).parent / 'pylbo/_version.py').resolve()
VERSION = None
with open(version_filepath) as f:
    for line in f.readlines():
        if line.strip().startswith('__version__'):
            VERSION = line.split("=")[-1].strip().strip('"')

setup(
    name = package_name,
    version = VERSION,
    description = 'Post-processing framework for the Legolas code',
    author = ['Niels Claes', 'Jordi De Jonghe'],
    author_email = ['niels.claes@kuleuven.be', 'jordi.dejonghe@kuleuven.be'],
    keywords = 'interface data-analysis',
    python_requires = '>=3.6',
    install_requires = required_packages,
    packages = find_packages()
)
