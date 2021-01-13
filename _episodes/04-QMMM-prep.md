---
title: "QM/MM Structure Preparation"
teaching: 60
exercises: 0
questions:
- ""
objectives:
- "???"
keypoints:
- "File numbering is a really common trouble spot, so pay close attention to
what you're working with and where."
---

**Highlights**
* TOC
{:toc}

1. Cluster the trajectory (see [phase 3](03-MD-cluster/index.html))
2. Write out an unstripped frame centered on the origin
    - Build a parameter file!
3. Convert to TINKER XYZ
4. Build the `regions.inp`, `connect.inp`, `tinker.key`, and `BASIS` file.
5. Run an SP
6. Do a DFP for the reactant ([phase 5](05-DFP-prods/index.html))
7. Build the product from the optimized reactant ([phase 5](05-DFP-prods/index.html))
8. Compare the product and optimized reactant ([phase 5](05-DFP-prods/index.html))
9. QSM ([phase 6](06-QSM/index.html))

## Running a Single Point Energy Calculation

```

```
{: .source}

{% include links.md %}
