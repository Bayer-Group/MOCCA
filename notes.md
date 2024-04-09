# Setting up automatic PyPI publishing

1. Create PyPI account
2.  Fill in the from at https://pypi.org/manage/account/publishing/

The code should be automatically published on every push to `main`.

# Setting up docs

To build docs locally, `pip install sphinx sphinx-rtd-theme` and then run `make html` in the `docs` directory.

### Automatic deployment to github pages
The docs are deployed automatically to GitHub pages on push to main. This is done by the `deploy_pages.yaml` workflow.

# Packaging example data
Example chromatogramd data don't fit into PyPI package size limit. The workflow `package_data.yaml` compresses them and pushes on `example-data` branch.

The data can be then downloaded into the package using `python -m mocca2 --download-data`.

If master repo changes, the link needs to be changed in `__main__.py` and user name and email should be changed in the workflow.

# Automatic tests
All automatic tests are in the `tests/` directory. They are automatically run using the `ci.yaml` workflow.

# Publishing the package to PyPI and GitHub releases
This is done by the `publish to pypu.yaml` workflow.

If master repo changes, it needs to be updated in the PyPI registry.