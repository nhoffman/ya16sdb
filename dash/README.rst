===============================================
NCBI Species Level 16s Reference Sequence Plots
===============================================

Web Site: https://deenurp.org/plots
Github: https://github.com/nhoffman/ya16sdb

Authors: Noah Hoffman, Chris Rosenthal

.. image:: screenshot.png

summary of features
===================

Users can browse Bacteria and Archea species using the Genus and
Species dropdown menus.  Users can also search for sequences by
seqname, species taxonomy id, or species taxonomy name using the
search box or by url arguments.  Users can modify scatter plot markers
by color and shape and use the Selection and Visibility drop down menus
to increase the size or visibility of the markers respectively.  Both
the Selection menu and plot selection tool will update the data table
below the plot layout and Axes dropdown menus.

standalone dev environment
==========================

Note that the dev application uses data in
``filter_details.feather.gz`` by default.

Set up::

  python3 -m venv py3-env
  source py3-env/bin/activate
  pip install -U pip
  pip install -r requirements.txt

Run the app::

  python app.py

The application will be served at ``http://127.0.0.1:8050``

Docker
======

Building a Docker image::

  docker build --platform=linux/amd64 . -t ya16sdb

Run the application from Docker::

  docker run --platform=linux/amd64 --rm -p 8000:8000 ya16sdb

