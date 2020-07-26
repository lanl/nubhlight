nubhlight open-sourcing description
===

# Description

nubhlight is a code for simulating disks of materials around black
holes formed after the in-spiral and merger of two-black holes. It is
based on the open-source code bhlight. nubhlight models all relevant
physics required to solve the post-merger disk problem. This includes:
general relativity in a static spacetime, neutrino and weak physics,
relativistic magnetohydrodynamics, and the ability to use open-sourced
tabulated equations of state provided by stellarcollapse.org.

## Details

nubhlight is based on the open-source code bhlight. (bhlight is
available here: https://github.com/AFD-Illinois/ebhlight) bhlight has
the following capabilities and methods:

- Finite volumes are used to model general relativistic
  magnetohydrodynamics as described in Gammie, McKinney, and Toth,
  2006 (see: https://doi.org/10.1086/374594).

- A Monte Carlo method to model radiation, as described in Dolence et
  al., 2009 (see: https://doi.org/10.1088/0067-0049/184/2/387).

- The radiation and fluid sectors are are coupled together via
  operator splitting as described in Ryan, Dolence, and Gammie, 2015
  (see: https://doi.org/10.1088/0004-637X/807/1/31).

The following changes were made at LANL to turn bhlight into
nubhlight:

- The equation of state for the relativistic magnetohydrodynamics was
  extended to conserve electron number and to include open-source
  tabulated data provided by stellarcollapse.org

- Radiation packets were changed to model neutrinos of three flavors:
  electron neutrino, electron antineutrino, all other flavors grouped
  together as "heavy" neutrinos.
  
- Tracer particles were added

- All details of required changes can be found in Miller, Ryan, and
  Dolence, 2019 (see: https://doi.org/10.3847/1538-4365/ab09fc).

## Applications

nubhlight is useful for astrophysics research on black hole
accretion. It is a research code and has no other uses and has no
commercial applications.

# Development Status and History

## When

nubhlight was developed by author Jonah Miller during his postdoc at
LANL. Development work began in February 2018 and finished in March,
2019.

## Where

The software additions needed (to create the nubhlight product from
the existing bhlight product) were created at LANL.

## Mission or purpose for which the software was developed

Software additions needed to create nubhlight from bhlight were
undertaken to support fundamental research (as defined in NSDD-189)
for modelling disks of material formed around black holes after
neutron star mergers. The nubhlight code aids in open astrophysics
work performed by the code author.  This work supports an LDRD grant
on this fundamental research topic, as well as research for an Office
of Science SciDAC grant. The original SciDAC program plan can be found
here: https://www.scidac.gov/documents/SciDAC.pdf.

## Further detail

nubhlight was created from the open source software project bhlight by
extending it in several ways that were called out above in the code
details section.  The link for bhlight is also listed above and
details of its licensing can be found there.

## Funding

- DOE office of science SciDAC project

- LDRD through the Center for Nonlinear Studies

# Software Details

## Operating system

Linux only. Tested on Debian 4.0 and above.

## Compilers

GNU compiler only.

## Languages

nubhlight is written in C99 and Fortran90. However, the build system
and data analysis utilities uses python 3.

## Dependencies

nubhlight has the following software dependencies, which it links
against. (Source code for dependencies **not** included. Dependencies
must be pre-installed on target machine.)

- GNU scientific library (https://www.gnu.org/software/gsl/)

- Any implementation of the MPI protocol, such as OpenMPI (https://www.open-mpi.org/)

- The hierarchical data format, HDF5 (https://www.hdfgroup.org/downloads/hdf5/)

- python 3 (https://www.python.org/)

Additionally, data analysis tools provided along with nubhlight use
the scientific python software stack and rely on the following
libraries:

- h5py (https://www.h5py.org/)

- numpy (https://numpy.org/)

- scipy (https://www.scipy.org/)

- matplotlib (https://matplotlib.org/)
