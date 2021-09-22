# Yukawa Complexation Langevin Dynamics

The Yukawa potential can be used to simulate implicit salt effects, by using a screened potential akin to the Debye-Huckel approximation allows for this.

To do this however some new features (at the time of writing c. September 2021) need to be added to `starpolymers`. At the moment it is only possible to produce diblock copolymers and LAMMPS configurations using the `full` atom style, where charged particles are allocated an explicit charge. However instead each molecule must be given a subset of particle types to choose from representing each charge that is present within the molecule (or rather molecule type).

This space allows for experimentation and comparison between traditional kspace electrostatic force evaluations and equivalent systems that take advantage of the screened electrostatics that are often found in biological systems.