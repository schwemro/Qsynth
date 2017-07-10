Streamflow Generator (Qsynth)
=============================

**Streamflow generator** synthesises an artificial streamflow time series with daily values. The generation merely based on a observed streamflow time series. The implemented approach takes roughly a bearing described in Salas (1993)[^1]

License
---
This software can be distributed freely under the GPL v2 license. Please read the LICENSE for further information.

Copyright (c) 2017, Robin Schwemmle (<rschwemmle@yahoo.de>)

Structure
---
* example
* source_code
* GUI

Example
---
An example file how to use the generator is provided:

- *example/test.csv*: File which contains observed streamflow time series with daily values used as input for the streamflow generator. May also be used what your input file must look like. Otherwise streamflow will not be able to read the file. A blueprint is shown in the table below.

> **Date**           | **Q**
> -------------------|------
> dd-mm-YYYY       | x1
> dd-mm-YYYY       | x2
> dd-mm-YYYY       | x3
>  ...                        | ...

 [^1]: Salas, J.D., 1993. Analysis and modeling of hydrologic time series. In: Maidment, D.R. (Ed.), Handbook of Hydrology. McGraw-Hill, pp. 1-72.
