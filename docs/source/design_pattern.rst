=========================
MAZE Index Mapper
=========================
*********************************************************
Terminology
*********************************************************
The following documentation frequently references ASE's Atoms object, ASE's Atom object and the individual atoms within ASE's Atoms objects. Since this can get confusing, clarity the following definitions are provided for clarity:

* ase: short for the atomic simulation environment, which is what the MAZE package is based off of.
* ase.Atoms class: Refers to the class definition of ``Atoms`` in the ase package.
* ase.Atom class: Refers to the class definition of ``Atom`` in the ase package.
* ase.Atoms object: Refers to an instance of the class ``ase.Atoms``.
* ase.Atom object: Refers to an instance of the class ``ase.Atom``.
* atoms: (lowercase) refers to the atoms represented by an ``Atoms`` object or a ``Zeotype`` object.

*********************************************************
Motivation
*********************************************************
The class ``Zeotype``, inherits from the atomic simulation environment (ase)'s Atoms class. It also adds additional features which make tracking atoms easier.  To understand the MAZE code, it is essential to understand the motivation behind the code's design decisions.

One feature of the ase.Atoms class is that an instance of the class stores numpy arrays containing information needed to build Atom objects when the Atom objects are needed. This “on-the-fly” Atom object creation saves memory, but unfortunately, it makes relating different Atoms objects to each other challenging.

An example of this problematic behavior occurs when square brackets are used to access a specific Atom object from an Atoms object (i.e. ``first_atom = co[0]``). Every time this operation is performed, a new Atom object is created. A demonstration is shown below:

First we create a simple Atoms object.

	>>> from ase import Atoms
	>>> d = 1.1
	>>> co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)])
	>>> co
	Atoms(symbols='CO', pbc=False)

Next, we use the ``[index]`` (i.e. ``__getitem__``) operation to get the first atom from the co Atoms object. We then print the ID of the ``first_atom`` object.

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

This behavior is unintuitive because if ``co`` was a numpy array filled with objects, then both ``first_atom`` and ``first_atom_2`` would point towards the same object and ``id(first_atom) == id(first_atom_2)``. The reason for this strange behavior is that every time the ``__getitem__`` operation is performed a new Atom object is created, and thus the ID's of the two Atom objects are different.

We can check equality between objects in two different ways. One is with the ``is`` operation, which checks to see if the variables reference the same object. The other is with the ``==`` operation, which uses the ``__eq__`` method defined in the left class to check equality. If we use the ``is`` operation on our two Atom objects, we get ``False``.

	>>> first_atom is first_atom_2
	False

We already knew that ``id(first_atom) != id(first_atom_2)`` , so this result isn't surprising. However, the result of the ``==`` operation is surprising.

	>>> first_atom == first_atom_2
	False


Even though ``first_atom`` and ``first_atom_2`` came from the same Atoms object and share the same values for all of their properties, the ``==`` operation returns ``False``! This is more unexpected behavior.

The ``__getitem__`` operation can be used to select a subset of the parent Atoms object. If a list of indices is pasted in, the same equality behavior is encountered:

    >>> from ase import Atoms
    >>> d = 1.1
    >>> co2 = Atoms('CO', positions=[(-d/2, d/2, 0), (0, 0, 0), (d/2, d/2, 0)])
    Atoms(symbols='CO2', pbc=False)
    >>> co = co2[0:2]
    >>> co
    Atoms(symbols='CO', pbc=False)
    >>> co[0] == co2[0]
    False


The above example demonstrated some unusual behavior, now we will examine why this is an issue. Let's take a look at a typical Zeolite simulation workflow:

#. Load a Zeolite cif file into an Atoms object
#. Identify the T sites in that Atoms object
#. For each unique T site:


    a. Create a new Atoms object consisting of a T site and the Atoms adjacent to the T site (e.g. a 5T cluster)
    #. Remove an Si atom from the new Atoms object
    #. Optimize the structure of the new Atoms object
    #. Integrate the optimized structure back into the original Zeolite Atoms object
    #. Save the altered Zeolite object as a .traj file with a unique name

Steps 1 and 2 are challenging, and the MAZE package presents solution in the form of the ``Zeolite.make`` method, which reads a ``cif`` file and keeps track of the unique T sites in a dictionary. Achieving this with the base ASE package is not trivial.

In the sub-steps of part 3, the problem with "on-the-fly" object creation emerges. Part a,b c are can be completed with the base ASE package. However, part d is not, because there is no way to map the optimized Atoms structure back into its parent Zeolite structure. The MAZE project solves this sub-Atoms mapping problem with the use of a custom Index Mapper.

*********************************************************
The Index Mapper Solution
*********************************************************

Zeotype simulation workflows frequently involve extracting atoms and adding atoms. This is challenging with ase because unique identities of the atoms are not stored. MAZE solves the identity storage problem by creating an ``IndexMapper`` object, which is a table that stores the mapping between indices of a parent ``PerfectZeolite`` object and all ``Zeolite`` objects derived from the parent ``PerfectZeolite``. The IndexMapper.main_index can be thought of as a table that looks like this:


+------+---------+-------------------+----------+
| main | parent  |  Zeolite_2        |Cluster_3 |
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
 

The implementation of the ``IndexMapper.main_index`` is a dictionary of dictionaries, where the keys for the parent dictionary are the main indices, and the keys for the sub-dictionaries are the names of the ``Zeolite`` objects. The values of the sub-dictionaries are the indices of the parent. For example, the above table would be represented as the following nested dictionary:

.. code-block:: json

	{0: {‘parent’:0, 'Zeolite_2':0, ‘Cluster_3’:None},
	1: {‘parent’:1, 'Zeolite_2':1, ‘Cluster_3’:None},
	2: {‘parent’:2, 'Zeolite_2':2, ‘Cluster_3’:None},
	….
	100: {‘parent’: None, ‘Zeolite_2’:99, ‘Cluster_3’:None}}


To keep this mapping straight, a functional programing-like interface is added for creating and removing atoms from a ``Zeolite`` object. When atoms are added or removed from the ``Zeolite`` object, a copy of the object being operated on is returned, rather than modifying the original object. Thus, the ``add_atoms`` and ``delete_atoms`` methods of the ``Zeolite`` class return new ``Zeolite`` objects with the user-specified modifications. These methods also add another column to the main_index corresponding to the newly created ``Zeolite``.

Note that when Python's garbage collector deletes an ``Zeolite`` object, the object is deregistered from the ``main_index`` table.

The additional functionality of the MAZE code builds on the ``add_atoms`` method and the ``delete_atoms`` method. The ``delete_atoms`` method is described in detail in the following section.


*********************************************************
The ``delete_atoms`` Method
*********************************************************
Z
The delete_atoms method returns a copy of the original ``Zeolite`` with the specified atoms deleted. An example of the delete method is shown below:

.. code-block:: bash

   >>> import maze
   >>> from maze import Zeolite
   >>> my_z = Zeolite.make('BEA')
   >>> atom_indices_to_delete = [i for i in range(0, 50)]  # make a list from 0 to 49
   >>> my_new_z = my_z.delete_atoms(atom_indices_to_delete)  # make a new iz with the first 50 atoms deleted
   >>> print('my_z has', len(my_z), 'atoms in it')
    my_z has 192 atoms in it
   >>> print('my_new_z has', len(my_new_z), 'atoms in it')
    my_new_z has 142 atoms in it


The ``my_new_iz`` object has 50 less atoms than the original ``my_iz`` object. At first, this does not appear any different than using the ``del`` operator on an ase.Atoms object.


However, there is now a new ``index_mapper`` object, which shows the relationship between all of the zeolite objects in the program.

.. code-block:: python

    # insert script above here
    import pandas as pd
    index_mapping_dataframe = pd.DataFrame(my_new_iz.index_mapper.main_index).T
    zeolites = [my_zeolite, my_iz, my_new_iz]
    zeolites_names = ['my_zeolite', 'my_iz', 'my_new_iz']
    for name, var in zip(zeolites_names, zeolites):
        print(name + '.name', var.name)
    print("DATAFRAME")
    print(index_mapping_dataframe)

output

.. code-block:: bash

    my_zeolite.name parent
    my_iz.name ImperfectZeotype_1
    my_new_iz.name ImperfectZeotype_2
    DATAFRAME
   ========  ====================  ====================
   parent    ImperfectZeotype_1    ImperfectZeotype_2
   ========  ====================  ====================
   0         0                     NaN
   1         1                     NaN
   2         2                     2
   3         3                     3
   4         4                     4
   ...       ...                   ...
   187       187                   137
   188       188                   138
   189       189                   139
   190       190                   140
   191       191                   141
   ========  ====================  ====================


With this mapping, we can alter the Imperfect Zeolite with fewer atoms and then integrate it back into the larger Zeolite.

To offer further insight into how the ``delete_atoms`` method works, let's examine the source code:

.. code-block:: python

       def delete_atoms(self, indices_to_delete) -> 'ImperfectZeotype':
           """Delete atoms from imperfect zeotype by returning a copy with atoms deleted

           :param indices_to_delete: Indices of atoms in current zeotype to delete
           :return: a copy of self with atoms deleted
           """
           new_self_a = ase.Atoms(self)
           del new_self_a[indices_to_delete]
           new_self = self.__class__(new_self_a)
           self.set_attrs_source(new_self, self)
           old_to_new_map = self._get_old_to_new_map(self, new_self)
           self.index_mapper.register(self.name, new_self.name, old_to_new_map)
           return new_self

We will now go through this line-by-line. The first line uses ``ase.Atoms`` initialization method to build an Atoms object that contains all of the atoms of the imperfect zeolite being operated on, but none of the additional information encoded in the imperfect zeolite object. The point of this step is to create a simple copy of ``self``, without all of the complexities added by the ``ImperfectZeolite`` object. The ``ase.Atoms`` initialization method is analogous to ``deepcopy`` in that there is no shared information between ``self`` and ``new_self_a``.

The next line deletes the atoms using the ``del`` operation on the new_self_a object. The side effects of this operation are contained to the ``new_self_a object``. After the ``del`` operation, a ``new_self`` is built using the ``self.__class__`` method. This is used so that a subclass will return another copy of itself rather than an ``ImperfectZeolite`` object.

After ``new_self`` is created its attributes are set to that of its source. It is important that ``new_self`` share the same ``index_mapper`` and ``parent_zeotype`` as its source. The one attribute difference will be its name, which is uniquely set during initialization.

Now comes the registration. First an ``old_to_new_map`` is created which maps the indices in ``self`` to those in ``new_self``. This mapping is done based on the position of the atoms, which have not changed during the delete operation. Second, this ``old_to_new_map`` is used in conjunction with the ``self.index_mapper.register`` method to add another column to the table corresponding to the ``new_self`` object. After registration, this ``new_self`` object is finally returned.

This ``delete_atoms`` method is used in the initialization of ``Cluster`` and ``OpenDefect`` objects.


*********************************************************
Conclusion
*********************************************************

Hopefully this guide on the Index Mapper Design pattern was useful. For more details read through the corresponding source code.