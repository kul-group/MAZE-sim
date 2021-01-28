.. MAZE Project documentation master file, created by
   sphinx-quickstart on Fri Oct  9 10:52:14 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MAZE Project's documentation!
===========================================
This project aims to extend ASE to more naturally represent the properties of zeolites and facilitate the calculations used to determine their properties. The main functionality of this code creates classes which represent zeolites. These zeolite classes inherit from the ASE Atoms object, and can be treated the same way. They also include additional functionalities such as tagging the unique crytallographic sites when building from a properly labeled cif file. For a complete description of the code see the design pattern section. If you would prefer to learn how to use the code please see the demo section. For more information on the importance of zeolites and the motivation for writing this code see the background section.


.. toctree::
   :maxdepth: 3

   installation
   background
   tutorials
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
