.. include:: /defs.hrst
.. _how_to_cite:

How to cite
===========

|Finesse| can be cited in publications, technical notes, or any other work by
referencing the projec level Zenodo DOI https://doi.org/10.5281/zenodo.821363.
This link will always take you to the latest version of the DOI. Version
specific ones can also be found on Zenodo.


BibTeX
------

For those using BibTeX to manage their references, the following BibTeX
entry can be used to cite the project:

.. jupyter-execute::
    :hide-code:

    import textwrap
    import requests
    headers = {'accept': 'application/x-bibtex'}
    response = requests.get('https://zenodo.org/api/records/12662017',
                            headers=headers)
    response.encoding = 'utf-8'
    zenodo_record = response.text

    print(zenodo_record)
