Background
====================================================================


Zeolites are microporous aluminosilicate materials with important applications in catalysis and chemical separations. The atomic structure of a zeolite determines its chemical properties, which means that finding zeolite structures ideal for important applications is essential. This is complicated by the large number of known zeolite structures. Currently, over 245 unique zeolites have been synthesized, which is  a small fraction of the millions theoretically stable zeolites that have been identified. Exploring this vast space of different zeolite structures can be accelerated using chemical simulations, which can predict the properties of different zeolites in order to identify promising targets for laboratory experiments.

The Atomic Simulation Environment (ASE) is a set of tools that facilitates chemical simulations. ASE offers an object-oriented interface to represent and manipulate groups of atoms while providing an interface to other chemical simulation software. Many successful zeolite screening projects have been carried out with this ASE to date. However, existing chemical simulation software is not well suited to zeolite-specific tasks.

The **M**\ ultiscale **A**\ tomistic **Z**\ eotype Simulation **E**\ nvironment (MAZE) aims to streamline zeolite simulation workflows, by building upon ASE. MAZE utilizes an object-oriented programming approach built around an "ImperfectZeolite" class that

1. overcomes the cumbersome indexing issues associated with adding and removing atoms (e.g., creating defects, inserting TMs),

2. simplifies multiscale calculations through custom-built wrappers for different codes and

3. maintains the provenance of each calculation.

The third functionality is especially important for modern computational studies, as every calculation step (e.g. DFT, AIMD, force field etc.) can be easily traced back to the IZA CIF file.

For examples on the MAZE package please see the Demos section.