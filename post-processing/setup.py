from setuptools import setup, find_packages

package_name = "pylbo"
required_packages = ['numpy', 'matplotlib', 'f90nml']
version = '0.1'

setup(
    name = package_name,
    version = version,
    description = 'Post-processing framework for the Legolas code',
    author = ['Niels Claes', 'Jordi De Jonghe'],
    author_email = ['niels.claes@kuleuven.be', 'jordi.dejonghe@kuleuven.be'],
    keywords = 'interface data-analysis',
    python_requires = '>=3.6',
    install_requires = required_packages,
    packages = find_packages()
)
