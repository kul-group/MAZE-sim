******************************************************
Installation
******************************************************

First, make sure you have Atomic Simulation Environment installed.

Next, clone the MAZE-sim repo with the git clone command.

.. code-block:: bash

    $ git clone https://github.com/kul-group/MAZE-sim.git

Now install the package using pip, by typing the following commands:

.. code-block:: bash

    cd MAZE-sim
    pip install .

Finally, verify that the install was successful by running a python shell in your terminal and importing the MAZE package.

.. code-block:: bash

    python
    >>> import maze
    >>> maze.zeotypes.PerfectZeolite()
    Zeotype(symbols='', pbc=False)
    python
    >>> import maze
    >>> maze.zeotypes.PerfectZeolite()
    Zeotype(symbols='', pbc=False)
    python
    >>> import maze
    >>> maze.zeotypes.Zeotype()
    Zeotype(symbols='', pbc=False)
