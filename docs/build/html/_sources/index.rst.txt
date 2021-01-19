.. MAZE Project documentation master file, created by
   sphinx-quickstart on Fri Oct  9 10:52:14 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MAZE Project's documentation!
===========================================
This project aims to extend ASE to more naturally represent the properties of zeolites and facilitate the calculations required to determine their properties. The main functionality of this code comes about by creating classes which represent zeolites. These zeolite classes inherit from the ASE Atoms object, and can be treated as one. They also include additional functionality such as tagging the unique sites when building from a properly labeled cif file. For a complete description of the code see the docs. If you would prefer to learn from example please see the demo section of the docs. For more information on the importance of zeolite and the motivation for writing this code see the background section. 


.. toctree::
   :maxdepth: 3

   background
   demos
   design_pattern
   key_classes
   source
   modules
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
