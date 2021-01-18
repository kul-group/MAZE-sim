Multiscale Atomic Zeolite Simulation Environment (MAZE)
========================================================

This project aims to extend the Atomic Simulation Environment (ASE) to more naturally represent the properties of zeolites and facilitate the calculations required to determine their properties. 

The main functionality of this code comes about by creating classes which represent zeolites. These zeolite classes inherit from the ASE Atoms object, and can be treated as ASE Atoms objects. The following code demonstrates how a zeotype object can be treated as an ASE Atoms object. 
```
import sys
sys.path.insert(0, "<absolute path to zeotype folder>")
from zeotype import Zeotype   # import Zeotype class  
from ase import Atoms 
from ase.io import view 
zeo = Zeotype.build_from_cif_with_labels(“CHA.cif”)
if issubclass(type(zeo), Atoms):
    print(“zeo is a subclass of ase.Atom’s object”)   # this will print

# you can use zeo as you would any other ASE atoms object 

print(zeo.get_positions()) # prints out a position of all of the Atoms 
view(zeo)  # you can view your zeotype object as you would an Atoms object 
```

Cif files downloaded from the iza structure database (http://www.iza-structure.org/databases/ ) contain additional 
information about the unique sites in a zeolite. This information is lost when loading a cif file using
ase.Atoms.io.read method. Thus, a new ```build_from_cif_with_labels``` static method was created, which loads a cif 
file while retaining information about the atom labels. These labels can be accessed by using the Atom object tags or 
the dictionaries```site_to_atom_indices`` and ```atom_indices_to_site``` . 

Features 
=======
This improved, unique atom site labeling is just one of the many features added by the Zeotype class. Additional features include

-	Support to easily extract clusters from zeotype objects and then reincorporate them in once they have been optimized 
-	Cap atoms in clusters so that calculations can be performed on them 
-	Identify different types of atoms in a zeotype 
-	Download cif files from the zeotype database by calling a python function 
-	Add adsorbate and remove adsorbate atoms into the zeotype 


Installation
------------
Clone this github repo
Add the cloned repo to your python path with 
```python
import sys
sys.path.insert(0, "<absolute path to parent zeotype folder>")
sys.path.insert(0, "<absolute path to parent zeotype>")
```

For example,
```python
import sys
sys.path.insert(0, "/Users/myusername/Code")
sys.path.insert(0, "/Users/myusername/Code/zeotype")
```


Contribute
----------

- Issue Tracker: github.com/kul-group/zeotype/issues
- Source Code: https://github.com/kul-group/zeotype

Support
-------

If you are having issues, please let us know.
We have a mailing list located at: ddantonio@ucdavis.edu

License
-------

The project licensed -------
