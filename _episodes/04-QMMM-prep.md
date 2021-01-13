---
title: "QM/MM Structure Preparation"
teaching: 60
exercises: 0
questions:
- ""
objectives:
- "???"
keypoints:
- "Parameters are, and always be, the downfall of the computational chemist."
- "File numbering is a really common trouble spot, so pay close attention to
what you're working with and where."
---

**Highlights**
* TOC
{:toc}

> ## Overview of the QM/MM Process
>
> 1. Set-up and run MD
> 2. Cluster the trajectory (see [phase 3](03-MD-cluster/index.html))
> 3. Write out an unstripped frame centered on the origin
>    - Build a parameter file!
> 4. Convert to TINKER XYZ
> 5. Build the `regions.inp`, `connect.inp`, `tinker.key`, and `BASIS` file.
> 6. Run an SP
> 7. Do a DFP for the reactant ([phase 5](05-DFP-prods/index.html))
> 8. Build the product from the optimized reactant ([phase 5](05-DFP-prods/index.html))
> 9. Compare the product and optimized reactant ([phase 5](05-DFP-prods/index.html))
> 10. QSM ([phase 6](06-QSM/index.html))
{: .checklist}

## The Trials ~~and Tribulations~~ of PDB to XYZ

TINKER uses a specific type of XYZ file, aptly referred to as a TINKER XYZ.

```
EXAMPLE LINE
```

You need to start with a Tinker XYZ file to build the relevant LICHEM XYZ and
connect files.
These two files are read into LICHEM (with some other information), and
LICHEM will then handle any inter-program file conversions.

If you're using AMOEBA as a force field after running AMBER MD, you will need
to have a set of AMOEBA parameters for any non-standard residues.
Since this tutorial involves the non-standard residue CO2, and we do not have
AMOEBA parameters for it, we will be using the point charge method for QM/MM.

You can use Tinker's `pdbxyz` program, but the program tends to crash
anytime something goes wrong, leaving you without any part of the new
Tinker XYZ, even if there was only one bad atom.
For this reason, several of us have coded ways around this program.
Here, we will use the
[`generate_TINKER_parameters.py`](https://github.com/emleddin/pdbxyz-xyzpdb)
script in conjuction with the
[`pdbxyz4amber-pmd-params.py`](https://github.com/emleddin/pdbxyz-xyzpdb/blob/main/pdbxyz4amber-pmd-params.py).
This script will use the original `prmtop` file to pull the same paramters that
were used for the AMBER MD and write them into a custom Tinker parameter file.
Because it takes the existing parameters, the non-standard residues are covered.

I have a few other means for the conversion and back-conversion of PDBs and
Tinker XYZs in my
[`pdbxyz-xyzpdb` repository](https://github.com/emleddin/pdbxyz-xyzpdb).
Mark has also written a Python package,
[`PDBTinker`](https://github.com/markahix/PDBTinker) to convert PDB files
to a Tinker XYZ using an AMOEBA force field.

## Generating the TINKER Parameter File

As stated, since we're using AMBER for the QM/MM, we will be using the
[`generate_TINKER_parameters.py`](https://github.com/emleddin/pdbxyz-xyzpdb)
script to build our parameter file.
This parameter file is independent of the TINKER version used.

Modify the options above the `Definitions` line, in this case `source_params`
and `param_file_name`.

~~~
## !!! This uses Pandas 1.0 !!! Without it, remove 'ignore_index' in the
## drop_duplicates lines.

import parmed as pmd
import pandas as pd
import numpy as np
import copy
from collections import OrderedDict

## Code to source a single parm file (not a leaprc)
# source_params = "parm99.dat"
# param_dat = pmd.load_file(source_params)

## Prmtop Method
source_params = "5Y2S_wat_fix.prmtop"
param_dat = pmd.load_file(source_params)

## It looks like anything sourced in the leaprc needs to have an absolute
## path to it, so you might need to modify the leaprc to incorporate the
## absolute path. It is STRONGLY RECOMMENDED that you copy the leaprc file to
## do this, and then reference that copy!!

## Leaprc Method
# source_params = "param_files/leaprc.ff14SB.OL15.tip3p"
# param_dat = pmd.amber.AmberParameterSet().from_leaprc(source_params)

## Leave the X dihedrals as X to fix by hand
## Setting this false will likely result in a MAXPRM issue with TINKER
## Because... well... you'll likely have over a million dihedral angles.
leave_as_X = True

## Give your FF a name (it will be preceded by AMBER-)
ff_name = "ff14SB"

## Give your new parameter file a name
param_file_name = "5Y2S_TINKER_params.prm"
~~~
{: .language-python}

Then run the script.
~~~
$ python3 generate_TINKER_parameters.py
~~~
{: .language-bash}

## Converting the PDB

Now that we have a functional parameter file, we can convert our PDB using
[`pdbxyz4amber-pmd-params.py`](https://github.com/emleddin/pdbxyz-xyzpdb/blob/main/pdbxyz4amber-pmd-params.py).
This script takes the cluster PDB file (`infile`) and the paramater file
(`param_file`) as input.

Like before, we will modify the first few lines with the file names.
~~~
import parmed as pmd
import numpy as np
import pandas as pd
import re
import sys

infile="5Y2S_subclust_c0_frame_6.pdb"
outfile="5Y2S_subclust_c0_frame_6.xyz"

param_file="5Y2S_TINKER_params.prm"
atom_lines="atom-lines.txt"

test_csv="test.csv"
~~~
{: .language-python}

Then run the script.
~~~
$ python3 pdbxyz4amber-pmd-params.py
~~~
{: .language-bash}

This script will usually take about 3 minutes to run because it consists of a
bunch of if statements that each atom has to go through.

> ## Potential Script Issue
>
> As Python updates, so do certain errors.
> It appears that `UnboundLocalError` may not be how Pandas refers to
> match errors now, instead using `NameError`.
> If you get:
> > ~~~
> > Traceback (most recent call last):
> >  File "pdbxyz4amber-pmd-params.py", line 256, in <module>
> >    system = convert_names(system, lines, AMOEBA)
> >  File "pdbxyz4amber-pmd-params.py", line 177, in convert_names
> >    print(residue.name, atom.name, test_RN, test_name)
> >NameError: name 'test_RN' is not defined
> > ~~~
> {: .language-python}
> Change `except UnboundLocalError:` to `except NameError:` on line 177.
{: .warning}

When the script runs properly, it should print:
~~~
Processing as AMBER parameters.
If I didn't find residues, they'll be listed here:
    Residue Name | Atom Name | Search ResName | Search Atom Name

ZN ZN ZN ZN
~~~
{: .output}

This alerts you that `ZN` is not being selected properly based on either your
PDB file or the parameter file.

A quick check of our PDB file shows that the `ZN` atom name is offset by a
column, but this doesn't actually affect the script.
```
ATOM   4054  OXT LYS   257     -20.102 -10.326  14.057  1.00  0.00           O
ATOM   4055 ZN   ZN6   258       3.179  -0.916  -1.326  1.00  0.00          ZN
ATOM   4056  H1  WT1   259       5.491   0.309  -0.943  1.00  0.00           H
```

Our parameter file also doesn't look particularly alarming.
```
atom        435  38    ZN    "ZN6 ZN"                         30     65.41    0 !! GUESSED CONNECTION
```

For some reason, it just doesn't like the `ZN6`, likely due to an internal
attempt at removing the final number to search for matches.
If we wanted to, we could build an exception in for when we run this script
in the future.
This is the easiest option if you're going to be testing more than one
snapshot, since you can just keep reusing the script.
Uncomment the first batch of lines under `## Address problem residues!` and
modify them accordingly.
```
        ## Address problem residues!
        if atom_test.empty == True:
            if residue.name in ('ZN'):     ## We know it's ZN because of the print-out
                test_RN = 'ZN6'            ## It needs to be looking for ZN6 in the PDB
                test_name = 'ZN'           ## ZN6 has an atom name of ZN
                res_test = lines[lines.ResName == test_RN]
                atom_test = res_test[res_test.AtomName == test_name]
```
{: .language-python}

If we don't want to modify the Python script, we can directly add the type to
the Tinker XYZ.
In this approach, we change:
```
  4055  ZN       3.179000   -0.916000   -1.326000     0 ATOM TYPE NOT FOUND
```
to
```
  4055  ZN       3.179000   -0.916000   -1.326000     435
```
in the Tinker XYZ that the script created.

## Checking the XYZ with `analyze`

After generating the Tinker XYZ, it is absolutely **crucial** to verify that
there are no missing parameters for TINKER.
Parameter issues cause major headaches, and it's so much easier to run `analyze`
than try and debug a QM/MM calculation.
You must use the `analyze` corresponding to whatever version of TINKER you
intend on using, as it changes across versions.

Before running `analyze`, however, we must create a `tinker.key` file.
The `tinker.key` file tells TINKER what parameter file it needs to use.
This key file only needs to contain one line (for now).
```
parameters 5Y2S_TINKER_params.prm
```

After you've created the `tinker.key`, you can run `analyze`.
```
$ analyze 5Y2S_subclust_c0_frame_6.xyz
```
{: .language-bash}

In our case, we're missing parameters! Oh, joy! Because they're water
parameters, it prints a gazillion and you can't see what the actual issue is.
So, you can modify your analyze command, or save it to an output file.
```
$ analyze 5Y2S_subclust_c0_frame_6.xyz | head -n 50
```
{: .language-bash}

Now we can determine what it doesn't like!
```
Atoms with an Unusual Number of Attached Atoms :

Type           Atom Name      Atom Type       Expected    Found

Valence        1417-N5           345              3         2
Valence        1454-N6           362              3         2
Valence        1797-N7           375              3         2
Valence        4059-c1           439              0         2
Valence        4060-o            440              0         1
Valence        4061-o            441              0         1
Valence        4062-c1           439              0         2
Valence        4063-o            440              0         1
Valence        4064-o            441              0         1

Undefined Angle Bending Parameters :

Type                  Atom Names                   Atom Classes

Angle        4067-HW   4066-OW   4068-HW           39   44   39
Angle        4070-HW   4069-OW   4071-HW           39   44   39
Angle        4073-HW   4072-OW   4074-HW           39   44   39
```
{: .output}

It's upset about the connectivity of several atoms (those in residues
`HD4`, `HD5`, `HE2`, and `CO2`, which makes sense because they're not typical
protein residues).
The 3 atoms listed as part of `HD4`, `HD5`, and `HE2` are all the nitrogen
coordinated to the zinc atom.
Since the zinc coordination itself is not considered a bond, we can decrease
the expected number of connections in the parameter file.
```
atom        345  32    N5    "HD4 NE2"                         7     14.01    2 !! GUESSED CONNECTION
atom        362  33    N6    "HD5 NE2"                         7     14.01    2 !! GUESSED CONNECTION
atom        375  34    N7    "HE2 ND1"                         7     14.01    2 !! GUESSED CONNECTION
```
[Note: you can remove the `!! GUESSED CONNECTION` if you want, since it's just
a comment anyway.]

The inferred `CO2` coordination is not correct, so its expected number of
connections also need to be adjusted in the parameter file.
```
atom        439  41    c1    "CO2 C"                           6     12.01    2 !! GUESSED CONNECTION
atom        440  42    o     "CO2 O1"                          8     16.00    1 !! GUESSED CONNECTION
atom        441  42    o     "CO2 O2"                          8     16.00    1 !! GUESSED CONNECTION
```

The undefined water parameters also make sense, as stated in the README for
the `generate_TINKER_parameters.py` script.
Building the parameter file drops the angle information from the solvent
(`angle HW OW HW 100.00 104.52`), so we have to append it to the angle section.
```
angle        39   44   39     100.0     104.52
```

After modifying the parameter file with these changes, we can re-run `analyze`.

```
Enter the Desired Analysis Types [G,P,E,A,L,D,M,V,C] :  E

Total Potential Energy :            -123334.5645 Kcal/mole

Intermolecular Energy :             -118769.0709 Kcal/mole

Energy Component Breakdown :           Kcal/mole      Interactions

Bond Stretching                         795.9634          30331
Angle Bending                          2286.9420          20568
Improper Torsion                        589.9997           4968
Torsional Angle                           0.0000          11012
Van der Waals                         16585.5434      146715508
Charge-Charge                       -143593.0130      940839611
```
{: .output}

Great! No problems this time, and a nice, negative Charge-Charge energy.
Huzzah!!!!!!

## Building `regions.inp` and Converting to LICHEM

Now that we have a Tinker XYZ, we can make use the
[`create_reg.py` script](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools)
to create our `regions.inp` file.

This script is a skeleton and will be active-site specific, so we need to
figure out some things first!
[Insert reference to selecting QM region]

Our QM region for 5Y2S will include the `ZN6`, `HD4`, `HD5`, `HE2`, `WT1`,
`CO2`, any waters within a specific cutoff of the `ZN6` and `CO2`, and
`HIS 64` in our active site.
In several QM/MM papers, `HIS 64` (the *biologically relevant* name for what is
actually **residue 61** in our PDB file) is important for the reaction.

> ## Selected QM/MM Citations for Human Carbonic Anhydrase
>
> - Zheng, Y. J. and Merz Jr., K. M. *J. Am. Chem. Soc.* **1992**, *114*, 26, 10498–10507.
DOI: [10.1021/ja00052a054](https://doi.org/10.1021/ja00052a054)
> - Chen, H.; Li, S.; Jiang, Y. *J. Phys. Chem. A* **2003**, *107*, 23, 4652–4660.
DOI: [10.1021/jp026788k](https://doi.org/10.1021/jp026788k)
> - Riccardi, D. and Cui, Q. *J. Phys. Chem. A* **2007**, *111*, 26, 5703–5711.
DOI: [doi.org/10.1021/jp070699w](https://doi.org/10.1021/jp070699w)
{: .callout}

Now that we know what to include in our active site, we can modify the
`create_reg.py` script.

First, we modify the header information, selecting the zinc residue as the
center of our shell around the QM region.
```
orig_pdb="5Y2S_subclust_c0_frame_6.pdb"
tink_xyz="5Y2S_subclust_c0_frame_6.xyz"

## Atom number for center of active atom shell
shell_center=4055
## Did you use the index from VMD for the shell_center? If yes, set True
VMD_index_shell=False
```
{: .language-python}

After changing the header, we have to re-write the
`select_QM(system, shell_center)` function.
You can use either VMD numbering (which starts are zero) or Tinker XYZ/PDB
numbering (which start at 1), but **you must** use the proper syntax for each.
It is *incredibly* easy to get the numbers wrong.

The table below contains the assignments for the QM atoms to use in the script.
It uses the atom numbers within the Tinker XYZ/PDB.

|                       | QM Residue(s)  | Number    | Pseudobond (CA, CB) | Boundary                  | Charge |
|-----------------------|----------------|-----------|---------------------|---------------------------|--------|
|                       | XYZ/PDB number |           | XYZ/PDB number      | XYZ/PDB number            |        |
| ZN6 (4055)            | 4055           | 1         |                     |                           | 2+     |
| HD4 (1405 - 1421)     | 1409 - 1419    | 11        | CA 1407             | C, HA, N, O, HN           | 0      |
| HD5 (1442 - 1458)     | 1446 - 1456    | 11        | CA 1444             | C, HA, N, O, HN           | 0      |
| HE2 (1789 - 1805)     | 1793 - 1803    | 11        | CA 1791             | C, HA, N, O, HN           | 0      |
| WT1 (4056 - 4058)     | 4056 - 4058    | 3         | -                   | -                         | 0      |
| CO2 (4062 - 4064)     | 4062 - 4064    | 3         | -                   | -                         | 0      |
| HID 61 (940 - 956)    | 944 - 954      | 11        | CA 942              | C, HA, N, O, HN           | 0      |
| GLU 103 (1574 - 1588) | 1581 - 1586    | 6         | CB 1578             | C, HA, N, O, HN, HB2, HB3 | 1-     |
| WAT                   | ???            | WAT       |                     |                           | 0      |
|                       | Total:         | 114 + WAT | 5                   | 27                        | 1+     |

```
$ python3 create_reg.py
```
{: .language-bash}

This file will create files named `regions.inp_backup` and `BASIS_list.txt`.
The reason for naming the `regions.inp` file with `_backup` is because running
the LICHEM conversion will replace any existing files named `connect.inp`,
`regions.inp`, or `xyzfile.xyz` in the current directory with the newly
generated one.
In the case of `regions.inp`, this new copy would only contain
```
QM_atoms: 0
Pseudobond_atoms: 0
Boundary_atoms: 0
Frozen_atoms: 0
```
getting rid of all the hard work you went through to create it!

To actually convert the Tinker XYZ, we use:
```
$ lichem -convert -t 5Y2S_subclust_c0_frame_6.xyz -k tinker.key > conversion-1.log
```
{: .language-bash}

Because we're using Gaussian, we need a `BASIS` file.
We can use another LICHEM command to create a skeleton `BASIS` file.
```
$ lichem -convert -b regions.inp
```
{: .language-bash}
This `BASIS` file can be modified using the information in `BASIS_list.txt`,
which maps the atoms between VMD numbering, PDB numbering, and the Gaussian
`BASIS` numbering.
All of this helps us ensure that the atoms we want to treat with a higher
basis, or that need have the pseudopotential atoms applied, are correct!

!!!!! DO THE BASIS SECTION

## Running a Single Point Energy Calculation


```
tag=5Y2S_subclust_c0_frame_6_out

## -n number of processors
## -x input xyz
## -c connectivity file
## -r regions file
## -o output xyz file
## -l output log file
lichem -n 20 -x xyzfile.xyz -c connect.inp -r regions.inp -o ${tag}.xyz -l ${tag}.log
```
{: .language-bash}

If you have a positive energy, kill the job.
- If it's in the QM region, something went wrong with the QM portion--debug that.
- If it's in the MM region, something went wrong with the MM portion--debug that.
- If it's in both, something catastrophic happened and everything failed--start
debugging in the QM.

{% include links.md %}
