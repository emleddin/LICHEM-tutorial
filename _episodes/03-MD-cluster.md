---
title: "Molecular Dynamics and Clustering"
teaching: 60
exercises: 60
questions:
- "How do I perform a molecular dynamics simulation with AMBER?"
- "What are some of the different approaches to clustering?"
- "How do I cluster a simulation based on specific criteria?"
- "How do I do subclustering?"
objectives:
- "Take an prmptop/inpcrd set through production."
- "Describe some benchmarks for determining whether a system is equilibrated."
- "Explain the difference between the k-means and DBscan clustering algorithms."
keypoints:
- "An MD simulation can be broken down into 4 phases: minimization, heating,
equilibration, and production."
- "DBscan can be very finicky, but it is a better approach to clustering."
---

**Highlights**
* TOC
{:toc}

## Molecular Dynamics Background

This where words would go, if I had them.

## Clustering

Clustering is a way to group together different snapshots of a trajectory.
It can be done through random assignment, or using specific characteristics,
such as a bond angle or distance.

`cpptraj` can be used to cluster trajectory data.
The k-means algorithm can be used to sort all trajectory points into a number
of defined clusters.
The DBscan algorithm can filter through outliers, but it requires a bit more
effort and testing to form clusters.
[David Sheehan](https://dashee87.github.io/data%20science/general/Clustering-with-Scikit-with-GIFs/)
has a great explanation of the different clustering algorithms.

When clustering, use

```
trajin /absolute/path/to/the/file/WT_protein_system_wat_imaged_1-50.nc
trajin /absolute/path/to/the/file/WT_protein_system_wat_imaged_1-100.nc
trajin /absolute/path/to/the/file/WT_protein_system_wat_imaged_1-100.nc

autoimage

## PA-O3' Distance
distance d1 :476@PA :457@O3' out WT_protein_sys_PaO_dist.dat

## O3'-PA-O5' Angle
angle a1 :457@O3' :476@PA :476@O5' out WT_protein_sys_OPO_angle.dat

## kdist test
cluster C0 dbscan kdist 9 data d1,a1 sieve 10

## DBSCAN based on PA-O distance and OPO angle
cluster c1 dbscan minpoints 10 epsilon 2.2 data d1,a1 \
 pairdist WT_protein_sys_PaO_OPO_db_clust_pairs.dat \
 loadpairdist \
 info WT_protein_sys_PaO_OPO_db_clust_detail_info.dat \
 out WT_protein_sys_PaO_OPO_db_clustnum_v_time.dat \
 summary WT_protein_sys_PaO_OPO_db_clust_summary.dat \
 avgout WT_protein_sys_PaO_OPO_db_clust avgfmt pdb \
 cpopvtime WT_protein_sys_PaO_OPO_db_popvtime.dat
#sieve 10 random
```

## Subclustering

Because clustering (via DBScan, at least) is used to group frames into like
categories, there should be one or two clusters that arise with the proper
characteristics for the proposed reaction.
Thus, to more randomly select frames that match these criteria, you can
recluster those within the matching cluster (what we refer to as subclustering).

To do this more easily, we use a `grep` command to separate out the list of
frames into their individual clusters.

```bash
$ grep " 0" WT_protein_system_H11_rms_clustnum_v_time.dat > clust_num_0.txt
$ grep " 1" WT_protein_system_H11_rms_clustnum_v_time.dat > clust_num_1.txt
$ grep " 2" WT_protein_system_H11_rms_clustnum_v_time.dat > clust_num_2.txt
...
```

From the `clust_num_X.txt` file with the closest match, you can use Python to
print a list of frames to then use in your `subclust.in` file with `cpptraj`.
You don't have to use Python, it's just one way we've chosen to highlight.

```python
import pandas as pd
clust = pd.read_csv("clust_num_1.txt", header=None, delim_whitespace=True)
test = clust.loc[:,0].values.tolist()
test2 = ''.join(str(i)+"," for i in test)
f=open("out_clust1.txt","w+")
f.write(test2)
f.close()
```

You can then write a `subclust.in` file that writes out the relevant frames
to a new trajectory file, and cluster those frames separately.

```
## Read in all the files used for the first set of clustering so that
## The proper frames are pulled
trajin /absolute/path/to/the/file/WT_protein_system_wat_imaged_1-50.nc
trajin /absolute/path/to/the/file/WT_protein_system_wat_imaged_1-100.nc
trajin /absolute/path/to/the/file/WT_protein_system_wat_imaged_1-100.nc

autoimage
## Center the non-solvent residues -- VERY IMPORTANT FOR TINKER
center :1-477 origin mass

## Write out the frames of a single cluster, as identified through the
## `grep` command
trajout WT_protein_system_wat_clustnum_X.nc netcdf onlyframes \
 8,15,26,29,36,40,41,53,60,71,72,75,76,77,80,93,98,113,124,132,134,\
 145,150,162,166,169,182,185,192,195,200,201,206,211,219,220,223,225,228,\
 230,237,249,257

autoimage

## Use this next block to create that file and clear any currently loaded
## files
#####################
go

clear trajin

go
####################

## Read in the file with specific frames from the cluster from the file you
## just wrote
trajin WT_protein_system_wat_clustnum_X.nc

## PA-O3' Distance
distance d1 :476@PA :457@O3'

## O3'-PA-O5' Angle
angle a1 :457@O3' :476@PA :476@O5'

## 10 clusters of k-means based on PaO distance and OPO angle
cluster popo kmeans clusters 10 data d1,a1 \
info WT_protein_system_cX_PaO_OPO_km_clust_detail_info.dat \
out WT_protein_system_cX_PaO_OPO_km_clustnum_v_time.dat \
summary WT_protein_system_cX_PaO_OPO_km_clust_summary.dat \
avgout WT_protein_system_cX_PaO_OPO_km_clust avgfmt pdb \
cpopvtime WT_protein_system_cX_PaO_OPO_popvtime.dat
```
{: .source}

> ## Important Note!
> When you go to run this subclustering, DO NOT use the MPI version of cpptraj!
> It has issues with reading and writing individual frames!
{: .callout}

One advantage to subclustering is that you can use different criteria for the
first and second clustering passes.
For instance, if you have an angle and a distance that are likely important for
the first reaction step, you can do the inital clustering on that and the
subclustering using the angle and distance of a secondary reaction step.

## Writing a Frame for QM/MM

After performing subclustering, you can once again use `cpptraj` to write PDBs
of the selected frames.
These frames are the snapshots that will be used for QM/MM.

```
trajin WT_protein_system_wat_clustnum_X.nc

autoimage

## Center the non-solvent residues -- VERY IMPORTANT FOR TINKER
center :1-477 origin mass

## Write out the specific PDBs identified with clustering
trajout WT_protein_system_sc_c0_frame_22.pdb pdb onlyframes 81
trajout WT_protein_system_sc_c1_frame_96.pdb pdb onlyframes 96
trajout WT_protein_system_sc_c2_frame_57.pdb pdb onlyframes 57
trajout WT_protein_system_sc_c3_frame_14.pdb pdb onlyframes 14
trajout WT_protein_system_sc_c4_frame_110.pdb pdb onlyframes 110
trajout WT_protein_system_sc_c5_frame_141.pdb pdb onlyframes 141
```
{: .source}

> ## Important Note!
> When you run this input, DO NOT use the MPI version of cpptraj!
> It does not support writing individual frames!
{: .callout}

All of the autoimaging and centering between subclustering and writing a
frame may be redundant, but it's better to be redundant than forget it.

{% include links.md %}
