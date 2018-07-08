---

title: A personal cheatsheet
---

Things I do regularly and wast time on by not remembering how I did previously:

## Python

- Filter a data frame such that only rows with all non-zero elements are retained (e.g. to get core orthogroups from the gene count output of OrthoFinder):

```
core = counts.loc[(counts != 0).all(axis=1)]
```
