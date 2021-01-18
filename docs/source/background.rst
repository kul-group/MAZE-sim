Background
====================================================================


Zeolites are microporous materials integral to the modern world. Their roles include (but are not limited to) gas separators and catalysts. The structure of a zeolite determines its viability for industrial applications. Searching for optimal zeolite structures for specific applications is complicated by the vastness of the zeolite chemical space. Currently, over 245 unique zeolites have been synthesized, which is  minor compared to the one million theoretical zeolites that have been identified. Exploring this expansive space can be achieved by chemical simulation experiments, which can predict the properties of both synthesized and theoretical zeolites and identify promising targets for laboratory experiments.

The Atomic Simulation Environment (ASE) is a set of tools that facilitates molecular simulation experiments. ASE offers an object-oriented interface to represent individual and groups of atoms and calculators. Many successful zeolite screening projects have been carried out with this software to date. Nevertheless,As the currently available general-purpose python codes are not ideally-suited for zeolite-specific tasks.

The Multiscale Atomistic Zeotype Simulation Environment (MAZE) aims to streamline zeolite experiment workflows. MAZE utilizes an object-oriented programming approach built around an "ImperfectZeolite" class that (1) overcomes the cumbersome indexing issues associated with adding and removing atoms (e.g., creating defects, inserting TMs), (2) simplifies multiscale calculations through custom-built wrappers for different codes and (3) maintains the provenance of each calculation. The latter aspect is especially important for modern computational studies, as every calculation step (e.g. DFT, AIMD, force field etc.) can be easily traced back to the IZA CIF file.

For examples on the MAZE package please see the Demos section.