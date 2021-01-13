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
- "Run k-means after DBscan. You can select the centroids from k-means as the
random snapshots."
---

**Highlights**
* TOC
{:toc}

## Molecular Dynamics Background

This where words would go, if I had them.

## Analyzing Molecular Dynamics Simulations

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

```
trajin ../5Y2S_wat_imaged_1-50.nc
#trajin ../path/to/other/stripped/trajectories.nc

autoimage

## O-C Distance
distance d1 :259@O :261@C out 5Y2S_OC_dist.dat

## O-C-O Angle
## 2 kinds of OCO are possible!
angle a1 :259@O :261@C :261@O1 out 5Y2S_OCO1_ang.dat
angle a2 :259@O :261@C :261@O2 out 5Y2S_OCO2_ang.dat

## kdist test
cluster C0 dbscan kdist 9 data d1,a1,a2 sieve 10

## DBSCAN based on PA-O distance and OPO angle
cluster c1 dbscan minpoints 25 epsilon 2.2 data d1,a1,a2 \
 pairdist 5Y2S_OC_OCO_db_clust_pairs.dat \
 loadpairdist \
 info 5Y2S_OC_OCO_db_clust_detail_info.dat \
 out 5Y2S_OC_OCO_db_clustnum_v_time.dat \
 summary 5Y2S_OC_OCO_db_clust_summary.dat \
 avgout 5Y2S_OC_OCO_db_clust avgfmt pdb \
 cpopvtime 5Y2S_OC_OCO_db_popvtime.dat
```
{: .source}

Because clustering is memory-intensive, we use stripped trajectories.
Typically, you will use whatever stripped trajectory you've created as part
of your analysis workflow.
In this case, we use the associated `strip.5Y2S_wat_fix.prmtop` file for the
topology information.

Since the overall goal of clustering is to sort frames into clusters based on
critical criteria for the reaction, we can perform clustering on all of our
replicates at once.
Thus, you will generally use absolute paths to call all of the replicates
one after another.

## Subclustering

DBScan is is used to group frames into like categories, so
there should be one or two clusters that arise with the proper
characteristics for the proposed reaction.
In order to more randomly select frames that match these criteria, you can
recluster those within the matching cluster (what we refer to as subclustering).

To do this more easily, we use a `grep` command to separate out the list of
frames into their individual clusters.

```bash
$ grep " 0" 5Y2S_OC_OCO_db_clustnum_v_time.dat > clust_num_0.txt
$ grep " 1" 5Y2S_OC_OCO_db_clustnum_v_time.dat > clust_num_1.txt
$ grep " 2" 5Y2S_OC_OCO_db_clustnum_v_time.dat > clust_num_2.txt
...
```

From the `clust_num_X.txt` file with the closest match, you can use Python to
print a list of frames to then use in your `subclust.in` file with `cpptraj`.
You don't have to use Python, it's just one way we've chosen to highlight.

```python
import pandas as pd
clust = pd.read_csv("clust_num_0.txt", header=None, delim_whitespace=True)
test = clust.loc[:,0].values.tolist()
test2 = ''.join(str(i)+"," for i in test)
f=open("out_clust0.txt","w+")
f.write(test2)
f.close()
```

You can then write a `subclust.in` file that writes out the relevant frames
to a new trajectory file, and cluster those frames separately.
Use the `strip.5Y2S_wat_fix.prmtop` file for the topology information.

```
trajin ../../5Y2S_wat_imaged_1-50.nc

autoimage

## Center the non-solvent residues -- VERY IMPORTANT FOR TINKER
center :1-261 origin mass

## Write out the frames of a single cluster, as identified through the
## `grep` command
trajout 5Y2S_wat_subclust_num_0.nc netcdf onlyframes \
 232,524,1041,1259,1788,1841,1847,1856,1869,1873,2280,2293,2600,2663,2710,\
 2716,2717,2877,2977,2999,3130,3343,3360,3374,3606,3616,3618,3629,5036,5044,\
 5240,5244,5415,5427,5473,5998,6168,6334,6704,6750,7056,7863,7899,8685,8756,\
 8810,9344,9802,9911,9916,9941,9980,10747,10756,10804,10806,10818

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
trajin 5Y2S_wat_subclust_num_0.nc

## O-C Distance
distance d1 :259@O :261@C

## O-C-O Angle
## 2 kinds of OCO are possible!
angle a1 :259@O :261@C :261@O1
angle a2 :259@O :261@C :261@O2

## ## 10 clusters of k-means based on OC distance and OCO angles
cluster coco kmeans clusters 10 data d1,a1,a2 \
info 5Y2S_OC_OCO_sub_km_clust_detail_info.dat \
out 5Y2S_OC_OCO_sub_km_clustnum_v_time.dat \
summary 5Y2S_OC_OCO_sub_km_clust_summary.dat \
avgout 5Y2S_OC_OCO_sub_km_clust avgfmt pdb \
cpopvtime 5Y2S_OC_OCO_sub_popvtime.dat
```
{: .source}

> ## MPI Warning
> When you go to run this subclustering, DO NOT use the MPI version of cpptraj!
> It has issues with reading and writing individual frames!
{: .warning}

One advantage to subclustering is that you can use different criteria for the
first and second clustering passes.
For instance, if you have an angle and a distance that are likely important for
the first reaction step, you can do the inital clustering on that and the
subclustering using the angle and distance of a secondary reaction step.

## Writing a Frame for QM/MM

After performing subclustering, you can once again use `cpptraj` to write PDBs
of the selected frames.
These frames are the snapshots that will be used for QM/MM.
However, it's not straightfoward, because you perform clustering on stripped
trajectories, and you need unstripped trajectories for QM/MM.

Thus, we'll use a script similar to how we did subclustering in order
to write the frames we need.

First, load in all of the frames you used when you created the unstripped
trajectory.

After subclustering, you can use the `5Y2S_OC_OCO_sub_km_clust_summary.dat`
file to write the centroids of the 10 clusters identified through k-means.
If you rerun the analysis, these centroids should be different, as k-means
is random!

Use the `5Y2S_wat_fix.prmtop` file for the topology information.

```
## Read in all the files used for to create the stripped trajectories in the
## exact same order to ensure that the proper frames are pulled
trajin /absolute/path/to/the/file/for/replicate1/5Y2S_wat_md1.mdcrd
trajin /absolute/path/to/the/file/for/replicate1/5Y2S_wat_md2.mdcrd

trajin /absolute/path/to/the/file/for/replicate2/5Y2S_wat_md1.mdcrd
trajin /absolute/path/to/the/file/for/replicate2/5Y2S_wat_md2.mdcrd

trajin /absolute/path/to/the/file/for/replicate3/5Y2S_wat_md1.mdcrd
trajin /absolute/path/to/the/file/for/replicate3/5Y2S_wat_md2.mdcrd
## ... continue

autoimage

## Center the non-solvent residues -- VERY IMPORTANT FOR TINKER
center :1-261 origin mass

## Write out the frames of a single cluster, as identified through the
## `grep` command
trajout 5Y2S_wat_subclust_num_0_unstripped.nc netcdf onlyframes \
 232,524,1041,1259,1788,1841,1847,1856,1869,1873,2280,2293,2600,2663,2710,\
 2716,2717,2877,2977,2999,3130,3343,3360,3374,3606,3616,3618,3629,5036,5044,\
 5240,5244,5415,5427,5473,5998,6168,6334,6704,6750,7056,7863,7899,8685,8756,\
 8810,9344,9802,9911,9916,9941,9980,10747,10756,10804,10806,10818

## Use this next block to create that file and clear any currently loaded
## files
#####################
go

clear trajin

go
####################

## Read in the file with specific frames from the cluster from the file you
## just wrote
trajin 5Y2S_wat_subclust_num_0_unstripped.nc

autoimage

## Center the non-solvent residues -- VERY IMPORTANT FOR TINKER
center :1-261 origin mass

## Write out the specific PDBs identified with clustering
trajout 5Y2S_subclust_c0_frame_6.pdb pdb onlyframes 6
trajout 5Y2S_subclust_c1_frame_29.pdb pdb onlyframes 29
trajout 5Y2S_subclust_c2_frame_34.pdb pdb onlyframes 34
trajout 5Y2S_subclust_c3_frame_56.pdb pdb onlyframes 56
trajout 5Y2S_subclust_c4_frame_10.pdb pdb onlyframes 10
trajout 5Y2S_subclust_c5_frame_45.pdb pdb onlyframes 45
trajout 5Y2S_subclust_c6_frame_4.pdb pdb onlyframes 4
trajout 5Y2S_subclust_c7_frame_5.pdb pdb onlyframes 5
trajout 5Y2S_subclust_c8_frame_47.pdb pdb onlyframes 47
trajout 5Y2S_subclust_c9_frame_5.pdb pdb onlyframes 5
```
{: .source}

> ## MPI Warning
> When you run this input, DO NOT use the MPI version of cpptraj!
> It does not support writing individual frames!
{: .warning}

Start by reading in the file with the subclusters, make sure that they are all
autoimaged and recentered, and then write them out with the `onlyframes`
argument.

{% include links.md %}
