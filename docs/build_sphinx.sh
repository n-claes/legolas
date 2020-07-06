rm -rf _autosummary
rm -rf content/src-docs/sphinx
rm -rf sphinx_source/_rst_files
sphinx-apidoc -o sphinx_source/_rst_files ../post_processing/pylbo
sphinx-autogen -o _autosummary sphinx_source/_rst_files/*.rst
sphinx-build -b html sphinx_source content/src-docs/sphinx
