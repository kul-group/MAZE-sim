=========================
ZeoSE Index Mapper  
=========================

*********
Motivation 
*********

The zeotype code is built off of ASE Atoms’ object. It is important to understand where this code comes from and the motivation for it. 

The main problem that was faced when building the Zeotype code, was that the ase.Atoms object does not store Atom objects in itself. Instead, it stores numpy arrays that contain the information needed to build Atom objects when they are needed. This “on-the-fly” Atom object creation saves memory, but also leads to some unintuitive behavior. For example, every time the brackets are used to subscript an Atoms object, a different Atom object is created.  Furthermore, even when using the less stringent double equals (==) to check equality between two Atom objects, False is returned. 

An example of this unintuitive behavior is shown below. 

Creating Atoms object 

	>>> from ase import Atoms
	>>> d = 1.1
	>>> co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)])
	>>> co
	Atoms(symbols='CO', pbc=False)

Getting one copy of the first atom 

	>>> first_atom = co[0]
	>>> print(first_atom)
	Atom('C', [0.0, 0.0, 0.0], index=0)
	>>> id(first_atom)
	140571043655264

Getting a second copy of the first atom 

	>>> first_atom_2 = co[0]
	>>> print(first_atom_2)
	Atom('C', [0.0, 0.0, 0.0], index=0)
	>>> id(first_atom_2)
	140571043655424

checking if they have the same object id 

	>>> first_atom is first_atom_2
	False


checking if they are equal as defined in the Atom class 

	>>> first_atom == first_atom_2
	False


The biggest issue that this creates is that it is challenging to relate Atom objects and Atoms objects created from a “parent” Atoms object back to the “parent”. 

*********
Wordy Solution Description 
*********

Zeotype simulations workflows frequently involve extracting atoms and adding atoms. This is challenging with the atomic simulation environment because unique identities of the atoms are not stored. This can be partially solved by using the tagging feature, to assign unique tags to every atom, but utilizing the tagging feature for this purpose restricts its use to simply tracking atoms. 

The Zeotype project aims to solve this problem by creating an IndexMapper, which is a table that stores the mapping between the indices of a Zeolite object and all of the ImperfectZeolite objects derived from the parent Zeotype. The IndexMapper.main_index can be thought of as a table that looks like this 


+------+---------+-------------------+----------+
| main | parent  |  ImperfectZeotype1| Cluster2 |
+======+=========+===================+==========+
| 0    | 0       | 0                 |     0    |
+------+---------+-------------------+----------+
| 1    |   1     |      1            |     None |
+------+---------+-------------------+----------+
| 2    |   2     |        2          |     2    |
+------+---------+-------------------+----------+
| ...  |    ...  |   ...             |    ...   |
+------+---------+-------------------+----------+
| 100  |   None  |         99        |    None  |
+------+---------+-------------------+----------+
 

The implementation of the IndexMapper.main_index is a dictionary of dictionaries, where the key for the first dictionary is the main index, and the keys for the second dictionary is the names of the ImperfectZeolites. The values of the second dictionary are the indices of the parent. For example, the above table would be represented as the following nested dictionary


	{0: {‘parent’:0, ‘ImperfectZeotype1’:0, ‘Cluster2’:None},
	1: {‘parent’:1, ‘ImperfectZeotype1’:1, ‘Cluster2’:None},
	2: {‘parent’:2, ‘ImperfectZeotype1’:2, ‘Cluster2’:None},
	….
	100: {‘parent’: None, ‘ImperfectZeotype1’:99, ‘Cluster2’:None}}


To keep this mapping straight, a functional programing-like interface is added for creating and removing atoms from a Zeolite object. The main idea is that when Atoms are added or removed from the Zeolite object a copy of the object being operated on is returned, rather than modifying the original object. This has several advantages, which I will not go into here, but the main consequence is that the ``add_atoms`` and ``delete_atoms`` method of the ImperfectZeolite return a new ImperfectZeolite object with the applied modifications. These methods also take care of adding another column to the main_index corresponding to the newly created ImperfectZeolite. 

After no references to an object exist, Python’s garbage collected deletes it. During this deleting process, the Zeolite is deregistered from the main_index table. 

The additional functionality of the Zeolite code is based off of the bedrock of the ``add_atoms`` method and the ``delete_atoms`` method. 


*********
The ``delete_atoms`` Method 
*********

The ``delete_atoms`` method is simplier than the 




