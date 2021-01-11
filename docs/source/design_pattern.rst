=========================
ZeoSE Index Mapper  
=========================

*********************************************************
Motivation
*********************************************************
This zeotype code is built off of atomic simulation environment (ASE)'s Atoms object (not to be confused with ASE's Atom object). It is important to understand the motivation behind the code's design decisions.

One feature of the ase.Atoms object is it stores numpy arrays containing information needed to build Atom objects only when the Atom objects are needed. This “on-the-fly” Atom object creation saves memory, but unfortunately, it makes relating different Atoms objects to each other challenging.

One example of this unintuitive behavior is when square brackets are used to access a specific Atom object from an Atoms object (i.e. ``first_atom = co[0]``). Every time this operation is performed a new Atom object is created. A demonstration of this behavior is shown below:

First we create a simple Atoms object.

	>>> from ase import Atoms
	>>> d = 1.1
	>>> co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)])
	>>> co
	Atoms(symbols='CO', pbc=False)

Then we use the ``[index]`` (i.e. ``__getitem__``) operation to get the first atom from the co Atoms object. We then print the ID of the ``first_atom`` object.

	>>> first_atom = co[0]
	>>> print(first_atom)
	Atom('C', [0.0, 0.0, 0.0], index=0)
	>>> id(first_atom)
	140571043655264

After performing this operation again, and checking the ID of ``first_atom_2`` we notice that ``first_atom`` and ``first_atom_2`` have different IDs.

	>>> first_atom_2 = co[0]
	>>> print(first_atom_2)
	Atom('C', [0.0, 0.0, 0.0], index=0)
	>>> id(first_atom_2)
	140571043655424

This behavior is unintuitive because if ``co`` was a numpy array filled with objects, then both ``first_atom`` and ``first_atom_2`` would point towards the same object and ``id(first_atom) == id(first_atom_2)``. The reason for this strange behavior is that every time the ``__getitem__`` operation is performed a new Atom object is created, thus the ID's of the two Atom objects are different.

We can check equality between objects in two different ways. One is with the ``is`` operation, which checks to see if the variables reference the same object. The other is with the ``==`` operation, which uses the ``__eq__`` method defined in the left class to check equality. If we use the ``is`` operation on our two Atom objects, we get ``False``.

	>>> first_atom is first_atom_2
	False

We already knew that ``id(first_atom) != id(first_atom_2)`` , thus this result shouldn't be surprising. The result of the ``==`` operation is surprising.

	>>> first_atom == first_atom_2
	False


Even though ``first_atom`` and ``first_atom_2`` came from the same Atoms object and share the same values for all of their properties, the ``==`` operation returns ``False``! This is more unexpected behavior.

The ``__getitem__`` operation can be used to select a subset of the parent Atoms object, if a list of indices is pasted in, the same equality behavior is encountered:

    >>> from ase import Atoms
    >>> d = 1.1
    >>> co2 = Atoms('CO', positions=[(-d/2, d/2, 0), (0, 0, 0), (d/2, d/2, 0)])
    Atoms(symbols='CO2', pbc=False)
    >>> co = co2[0:2]
    >>> co
    Atoms(symbols='CO', pbc=False)
    >>> co[0] == co2[0]
    False


The above example showed some unexpected behavior, but it has not yet been made clear why this is an issue. To demonstrate why this is an issue, let us take a look at a typical Zeolite simulation workflow:

#. Load a Zeolite cif file into an Atoms object
#. Identify T sites in that Atoms object
#. For each unique T site:
    a. Create a new Atoms object consisting of T site and the Atoms surrounding the T site
    #. Remove a Si from the new Atoms object
    #. Optimize the structure of the new Atoms object
    #. Integrate the optimized structure back into the original Zeolite Atoms object
    #. Save the altered Zeolite object as a .traj file with a unique name

Steps 1 and 2 are challenging, and the ZeoSE package presents a helpful solution in the form of the ``Zeolite.build_from_cif_with_labels`` method, which reads a ``cif`` file and keeps track of the unique T sites in a dictionary. Achieving this with the base ASE package is not easy.

In part 3's sub-steps the problem with the "on-the-fly" object creation emerges. Part a,b c are doable with the base ASE package.  Part d is not, because there is no way to map the optimized Atoms structure back into its parent Zeolite structure. The ZeoSE project solves this sub-Atoms mapping issue through the use of an custom Index Mapper.

*********************************************************
The Index Mapper Solution Solution
*********************************************************

Zeotype simulation workflows frequently involve extracting atoms and adding atoms. This is challenging with  ASE because unique identities of the atoms are not stored. This code solves the identity storage problem by creating an IndexMapper, which is a table that stores the mapping between indices of a Zeolite object and all ImperfectZeolite objects derived from the parent Zeotype. The IndexMapper.main_index can be thought of as a table that looks like this:


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
 

The implementation of the IndexMapper.main_index is a dictionary of dictionaries, where the keys for the parent dictionary are the main indices, and the keys for the sub-dictionaries are the names of the ImperfectZeolites. The values of the sub-dictionaries are the indices of the parent. For example, the above table would be represented as the following nested dictionary:


	{0: {‘parent’:0, ‘ImperfectZeotype1’:0, ‘Cluster2’:None},
	1: {‘parent’:1, ‘ImperfectZeotype1’:1, ‘Cluster2’:None},
	2: {‘parent’:2, ‘ImperfectZeotype1’:2, ‘Cluster2’:None},
	….
	100: {‘parent’: None, ‘ImperfectZeotype1’:99, ‘Cluster2’:None}}


To keep this mapping straight, a functional programing-like interface is added for creating and removing atoms from a Zeolite object. When atoms are added or removed from the Zeolite object, a copy of the object being operated on is returned, rather than modifying the original object. Thus, the ``add_atoms`` and ``delete_atoms`` methods of the ImperfectZeolite return new ImperfectZeolite objects with the user-specified modifications. These methods also add another column to the main_index corresponding to the newly created ImperfectZeolite.

When Python's garbage collector deletes an ImperfectZeolite object, the object is deregistered from the main_index table.

The additional functionality of the Zeolite code is based off of the bedrock of the ``add_atoms`` method and the ``delete_atoms`` method. 


*********************************************************
The ``delete_atoms`` Method
*********************************************************

The ``delete_atoms`` method is simpler than the ``add_atoms`` method so it will be described first.

The delete_atoms method returns a copy of the original ImpefectZeolite with the specified atoms deleted.







