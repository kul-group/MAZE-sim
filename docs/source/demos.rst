===============
Demos
===============

To showcase the capabilities of the ZeoSE code various demos will be presented, along with comments describing how to use the functions shown in the demos.

******************************************************
Installing
******************************************************

As of writing this, the ZeoSE code cannot be installed via pip or conda. Thus, to use the ZeoSE code you must clone the ZeoSE repo and then manually add the location of the cloned repo to Python's path so that it can be imported.

First, pick a folder to clone the repo into. In this tutorial I will use the ``Code`` folder in my home directory. Then clone the ZeoSE the repo with the git clone command.

.. code-block:: bash

    $ cd ~/Code
    $ git clone https://github.com/kul-group/zeotype.git

If cloned corectly, you should receive similar output

.. code-block:: bash

    Cloning into 'zeotype'...
    remote: Enumerating objects: 196, done.
    remote: Counting objects: 100% (196/196), done.
    remote: Compressing objects: 100% (129/129), done.
    remote: Total 1061 (delta 80), reused 171 (delta 67), pack-reused 865
    Receiving objects: 100% (1061/1061), 5.88 MiB | 7.17 MiB/s, done.
    Resolving deltas: 100% (631/631), done.

Now perform the following bash commands to get the directory of the folder containing the ZeoSE repo, and the path of the cloned repo.

.. code-block:: bash

    $ pwd
    /Users/dda/Code
    $ cd zeotype
    $ pwd
    /Users/dda/Code/zeotype

Now create a new python script in a different directory and append both of those paths to your python.

You should now be able to import the zeolite modules.


.. code-block:: bash

    import sys
    repo_parent_dir = "/Users/dda/Code"
    repo_dir = "/Users/dda/Code/zeotype"
    sys.path.insert(0, repo_parent_dir)   # folder containing zeotype folder
    sys.path.insert(0, repo_dir)          # zeotype folder

    import zeotype
    print(zeotype)

My output is  ``<module 'zeotype' from '/Users/dda/Code/zeotype/__init__.py'>`` . Your output should be something similar.


******************************************************
Cif Fetching from Database of Zeolite Structures
******************************************************

The `Database of Zeolite Structures <http://www.iza-structure.org/databases/>`_ is a useful resource for Zeolite simulation experiments. It contains cif files for all of the synthesized zeolites, organized by their three letter zeolite code. Downloading them from the website is easy when working on a local machine, but challenging when working on a remote machine. To facilitate smoother workflows, a simple python function which downloads cif files from the database was created. An example of using this to download a few different cif files is shown below. Try it out!

First import the ZeoSE package.

    >>> zeose_repo_location = "/Users/dda/Code"
    >>> import sys
    >>> import os
    >>> sys.path.insert(0, zeose_repo_location)                          # folder containing zeotype folder
    >>> sys.path.insert(0, os.path.join(zeose_repo_location, "zeotype"))  # zeotype folder
    >>> import zeotype
    >>> import glob

Next we define a helper function which prints out all of the directories in the current working directory. This will help us to visualize what the function is doing.


    >>> def print_dirs():
    ...     print('dirs in cwd',  glob.glob('**/'))
    ...

Now we print the directories in our current directory using our newly created helper function.

    >>> print_dirs()
    dirs in cwd []

Now let us try downloading a cif file, with the zeotype.download_cif function. We will pick the zeolite `GOO`. By default the directory that the cif file is downloded to is `data`. If this 'data' directory doesn't exist, it is created.
    >>> zeotype.download_cif("GOO") # downloads "GOO.cif" to data/GOO.cif
    >>> print_dirs()
    dirs in cwd ['data/']
    >>> print('files in data dir', glob.glob("data/*"))
    files in data dir ['data/GOO.cif']

The function worked, and now the cif file is just where we want it. If we want to download it to a custom location, we can do that as well

    >>> zeotype.download_cif("OFF", data_dir="my_other_data")
    >>> print_dirs()
    dirs in cwd ['my_other_data/', 'data/']
    >>> print('files in my_other_data dir', glob.glob("my_other_data/*"))
    files in my_other_data dir ['my_other_data/OFF.cif']


******************************************************
Building a Zeolite from a cif file
******************************************************

A cif file downloaded from the IZA-SC Database of Zeolite Strucutres looks like this:

.. code-block:: text

    data_CHA
    #**************************************************************************
    #
    # CIF taken from the IZA-SC Database of Zeolite Structures
    # Ch. Baerlocher and L.B. McCusker
    # Database of Zeolite Structures: http://www.iza-structure.org/databases/
    #
    # The atom coordinates and the cell parameters were optimized with DLS76
    # assuming a pure SiO2 composition.
    #
    #**************************************************************************

    _cell_length_a                  13.6750(0)
    _cell_length_b                  13.6750(0)
    _cell_length_c                  14.7670(0)
    _cell_angle_alpha               90.0000(0)
    _cell_angle_beta                90.0000(0)
    _cell_angle_gamma              120.0000(0)

    _symmetry_space_group_name_H-M     'R -3 m'
    _symmetry_Int_Tables_number         166
    _symmetry_cell_setting             trigonal

    loop_
    _symmetry_equiv_pos_as_xyz
    '+x,+y,+z'
    '2/3+x,1/3+y,1/3+z'
    '1/3+x,2/3+y,2/3+z'
    '-y,+x-y,+z'
    '2/3-y,1/3+x-y,1/3+z'
    '1/3-y,2/3+x-y,2/3+z'
    '-x+y,-x,+z'
    '2/3-x+y,1/3-x,1/3+z'
    '1/3-x+y,2/3-x,2/3+z'
    '-y,-x,+z'
    '2/3-y,1/3-x,1/3+z'
    '1/3-y,2/3-x,2/3+z'
    '-x+y,+y,+z'
    '2/3-x+y,1/3+y,1/3+z'
    '1/3-x+y,2/3+y,2/3+z'
    '+x,+x-y,+z'
    '2/3+x,1/3+x-y,1/3+z'
    '1/3+x,2/3+x-y,2/3+z'
    '-x,-y,-z'
    '2/3-x,1/3-y,1/3-z'
    '1/3-x,2/3-y,2/3-z'
    '+y,-x+y,-z'
    '2/3+y,1/3-x+y,1/3-z'
    '1/3+y,2/3-x+y,2/3-z'
    '+x-y,+x,-z'
    '2/3+x-y,1/3+x,1/3-z'
    '1/3+x-y,2/3+x,2/3-z'
    '+y,+x,-z'
    '2/3+y,1/3+x,1/3-z'
    '1/3+y,2/3+x,2/3-z'
    '+x-y,-y,-z'
    '2/3+x-y,1/3-y,1/3-z'
    '1/3+x-y,2/3-y,2/3-z'
    '-x,-x+y,-z'
    '2/3-x,1/3-x+y,1/3-z'
    '1/3-x,2/3-x+y,2/3-z'

    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
        O1    O     0.9020    0.0980    0.1227
        O2    O     0.9767    0.3101    0.1667
        O3    O     0.1203    0.2405    0.1315
        O4    O     0.0000    0.2577    0.0000
        T1    Si    0.9997    0.2264    0.1051




An important piece of information in this file is the _atom_site_label (01, 02, ... T1, T2.. ect.) that is located in the first column of the cif file near the atom position information. This information about the atoms identities is lost when ``ase.io.read`` function is used to build an Atoms object form a cif file. Knowing the identity of the T sites is critical for zeolite simulation experiments. This issue inspired the creation of a custom cif reading function for the zeotype object, ``build_from_cif_with_labels`` which creates a zeolite object and labels the unique atoms, by tagging them, and storing the mapping between the ``atom_site_label`` and the atom indices in the dictionaries ``self.site_to_atom_indices`` and ``self.atom_indices_to_site``.

To demonstrate this feature, let us try building a Zeotype object from a cif file.

First import the zeotype pacakage
    >>> zeose_repo_location = "/Users/dda/Code"
    >>> import sys
    >>> import os
    >>> sys.path.insert(0, zeose_repo_location)                          # folder containing zeotype folder
    >>> sys.path.insert(0, os.path.join(zeose_repo_location, "zeotype"))  # zeotype folder
    >>> import zeotype

Download a cif file
    >>> zeotype.download_cif('CHA', data_dir='data') # Download CHA.cif
Then use the static method ``build_from_cif_with_labels``
    >>> my_zeolite = zeotype.Zeotype.build_from_cif_with_labels('data/CHA.cif')  # build from code

The Zeotype has been built. The atom idenity information is now stored in two dictionaries. Let's take a look at them:

    >>> print('site_to_atom_indices map', my_zeolite.site_to_atom_indices, sep='\n\n')
    site_to_atom_indices map

.. code-block:: json

    {'O1': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    'O2': [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
    'O3': [36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53],
    'O4': [54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71],
    'T1': [72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107]}

.. code-block:: python

    >>> print('atom indices to site map', my_zeolite.atom_indices_to_site, sep='\n\n')
    atom indices to site map

.. code-block:: json

    {0: 'O1', 1: 'O1', 2: 'O1', 3: 'O1', 4: 'O1', 5: 'O1', 6: 'O1', 7: 'O1', 8: 'O1', 9: 'O1', 10: 'O1', 11: 'O1', 12: 'O1', 13: 'O1', 14: 'O1', 15: 'O1', 16: 'O1', 17: 'O1', 18: 'O2', 19: 'O2', 20: 'O2', 21: 'O2', 22: 'O2', 23: 'O2', 24: 'O2', 25: 'O2', 26: 'O2', 27: 'O2', 28: 'O2', 29: 'O2', 30: 'O2', 31: 'O2', 32: 'O2', 33: 'O2', 34: 'O2', 35: 'O2', 36: 'O3', 37: 'O3', 38: 'O3', 39: 'O3', 40: 'O3', 41: 'O3', 42: 'O3', 43: 'O3', 44: 'O3', 45: 'O3', 46: 'O3', 47: 'O3', 48: 'O3', 49: 'O3', 50: 'O3', 51: 'O3', 52: 'O3', 53: 'O3', 54: 'O4', 55: 'O4', 56: 'O4', 57: 'O4', 58: 'O4', 59: 'O4', 60: 'O4', 61: 'O4', 62: 'O4', 63: 'O4', 64: 'O4', 65: 'O4', 66: 'O4', 67: 'O4', 68: 'O4', 69: 'O4', 70: 'O4', 71: 'O4', 72: 'T1', 73: 'T1', 74: 'T1', 75: 'T1', 76: 'T1', 77: 'T1', 78: 'T1', 79: 'T1', 80: 'T1', 81: 'T1', 82: 'T1', 83: 'T1', 84: 'T1', 85: 'T1', 86: 'T1', 87: 'T1', 88: 'T1', 89: 'T1', 90: 'T1', 91: 'T1', 92: 'T1', 93: 'T1', 94: 'T1', 95: 'T1', 96: 'T1', 97: 'T1', 98: 'T1', 99: 'T1', 100: 'T1', 101: 'T1', 102: 'T1', 103: 'T1', 104: 'T1', 105: 'T1', 106: 'T1', 107: 'T1'}

Depending on the situation one dictionary might be more useful than the other.

