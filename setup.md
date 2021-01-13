---
title: Setup
---

LICHEM requires an underlying installation of the programs you intend to
do calculations with.

QM Wrappers
- Gaussian (`Gaussian`, `g09`, `g16`)
- PSI4 (`PSI4`)
- NWChem (`NWChem`)

MM Wrappers
- Tinker (`TINKER` keyword will work for any 7 or 8 series; Tinker HP requires
    `poledit` from another version)
- LAMMPS (`LAMMPS`)

In the case of this tutorial, you will need access to
[Gaussian (either `09` or `16`)](https://gaussian.com/) and
[Tinker7.1](https://dasher.wustl.edu/tinker/).

## Installing TINKER

The Tinker executables and source code are
[available on their site](https://dasher.wustl.edu/tinker/).

Install Tinker 7.1.3.

> ## WARNING
>
> Tinker parameter files are different between 7.1 and various versions of
> Tinker 8!
> You **MUST** modify the parameter files associated with the version you are
> going to use for QM/MM!
{: .warning}

> ## Tip
>
> If you have multiple versions of Tinker installed, export a global variable
> to the bin of each edition name.
>
> Add this to your `~/.bashrc`
> ~~~
> export Tinker7=/home/$USER/bin/tinker7/bin
> ~~~
> {: .source}
> and then call the program with `$Tinker7/xyzedit`.
{: .callout}

## Gaussian

If you are using Gaussian on an HPC cluster, simply load the module for the
version you would like to use.
You may need to be added to the user's group first.

For TACC resources, this means you will need to fill out a
[usage agreement](https://portal.tacc.utexas.edu/software/gaussian).

## Installing LICHEM

The instructions for installing LICHEM are included here for convenience.

~~~
$ mkdir LICHEM
$ git clone https://github.com/CisnerosResearch/LICHEM.git ./LICHEM/
$ ./configure
$ make install
~~~
{: .language-bash}

## Python 3

This tutorial also uses Python 3 code to make setting up LICHEM simulations
easier.
We recommend an [Anaconda](https://www.anaconda.com/products/individual)
installation of Python 3.

- Video of [Windows Install](https://www.youtube.com/watch?v=xxQ0mzZ8UvA)
- Video of [MacOS Install](https://www.youtube.com/watch?v=TcSAln46u9U)
- [Using `conda` to install packages](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-pkgs.html#installing-packages)

The scripts we will use require the [`parmed`](https://github.com/ParmEd/ParmEd),
[`numpy`](https://numpy.org/install/), and
[`pandas`](https://pandas.pydata.org/getting_started.html)
packages to be installed.
The Anaconda installation should include `numpy` and `pandas`.

~~~
$ conda install -c omnia parmed
$ conda install numpy
$ conda install pandas
~~~
{: .source}

{% include links.md %}
