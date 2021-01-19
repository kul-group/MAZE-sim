**M**ultiscale **A**tomic **Z**eolite Simulation **E**nvironment (**MAZE**)
========================================================
Introduction
===============
This project aims to extend the Atomic Simulation Environment (ASE) to more naturally represent the properties of zeolites and facilitate the calculations required to determine their properties. 

The main functionality of this code comes about by creating classes which represent zeolites. These zeolite classes inherit from the ASE Atoms object, and can be treated as ASE Atoms objects. The following code demonstrates how a zeotype object can be treated as an ASE Atoms object. 

Cif files downloaded from the iza structure database (http://www.iza-structure.org/databases/ ) contain additional 
information about the unique sites in a zeolite. This information is lost when loading a cif file using
ase.Atoms.io.read method. Thus, a new ```build_from_cif_with_labels``` static method was created, which loads a cif 
file while retaining information about the atom labels. These labels can be accessed by using the Atom object tags or 
the dictionaries ``site_to_atom_indices`` and ``atom_indices_to_site``. 

Installation 
=================

1. Clone this git repository 
    ``` 
    git clone https://github.com/kul-group/MAZE-sim.git
    ```
2. Navigate into the cloned directory 
    ```
    cd MAZE-sim
    ```
3. Install the package with pip
    ```
    pip install . 
    ```
4. Verify the package is installed by running a python shell and importing the ``maze`` package. 
    ```
    python
   >>> import maze 
   >>> maze.zeotypes.Zeotype()
   Zeotype(symbols='', pbc=False)
    ```


Features 
=======
This improved, unique atom site labeling is just one of the many features added by the Zeotype class. Additional features include

-	Support to easily extract clusters from zeotype objects and then reincorporate them in once they have been optimized 
-	Cap atoms in clusters so that calculations can be performed on them 
-	Identify different types of atoms in a zeotype 
-	Download cif files from the zeotype database by calling a python function 
-	Add adsorbate and remove adsorbate atoms into the zeotype 

Check out the documentation here, for additional details on how to use this package. 


Contribute
----------

- Source Code: https://github.com/kul-group/MAZE-sim
- Issue Tracker: https://github.com/kul-group/MAZE-sim/issues

Support
-------

If you are having issues, please let us know.
We have a mailing list located at: dexter.d.antonio@gmail.com

License
-------

Copyright 2021 Dexter Antonio

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.