# Setting up automatic PyPI publishing

1. Create PyPI account
2.  Fill in the from at https://pypi.org/manage/account/publishing/

The code should be automatically published on every push to `main`.

# Setting up docs

To build docs locally, `pip install sphinx sphinx-rtd-theme` and then run `make html` in the `docs` directory.

Change the user name and email in the github workflow. Enable github pages in repo>settings>pages and set up "deploy from branch".