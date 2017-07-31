Streamflow Generator (Qsynth)
=============================

**Streamflow generator** synthesises an artificial streamflow time series with daily values. The generation merely based on an observed streamflow time series. The implemented approach takes roughly a bearing described in Salas (1993)[^1]

License
---
This software can be distributed freely under the GPL v2 license. Please read the LICENSE for further information.

© 2017, Robin Schwemmle (<rschwemmle@yahoo.de>)

Structure
---
* Directions for Use
* Example

Directions for Use
---
![](/Users/robinschwemmle/Desktop/MATLAB/Stream_Scale_Model_for_SHP_Efficiency/Qsynth/GUI/GUI.png "Figure 1")


>**Input**

>`dd/MM/yyyy`: Indicate date format of your time series

>`NA`: Indicate format of the missing values in your time series

>`Load .csv`: Imports with observed streamflow values into the application

> `Daily Streamflow`: Displays observed streamflow time series

> `Seasonality`: Daily and monthly runoff regime as well as parde coefficient are illustrated

---

>**Time Period**

>`Start Date`:

>`End Date`:

---

>**Generate Streamflow**

>`Run`:

---


>**Output**

>`Synthetic Streamflow`:

>`Filter`:

>`Autocorrelation`:

> `Histogram`:

> `IHA`:

> `Simple Test Statistic`:

> `Volume`:

> `Export`:

---

```flow
sub1=>subroutine: Input
sub2=>subroutine: Time Period
sub3=>subroutine: Generate Streamflow
sub4=>subroutine: Output

sub1->sub2->sub3->sub4->

op3=>operation: Export
```

```flow
sub1=>subroutine: Input
i1=>inputoutput: Date format
i2=>inputoutput: NA format
op1=>operation: Load .csv
sub1(right)->i1(right)->i2(right)->op1
```

```flow
sub2=>subroutine: Time Period
i3=>inputoutput: Start date
i4=>inputoutput: End date
sub2(right)->i3(right)->i4(right)
```

```flow
sub3=>subroutine: Generate Streamflow
op2=>operation: Run
sub3(right)->op2
```

```flow
sub4=>subroutine: Output
op3=>operation: Export
sub4(right)->op3
```

Example
---
An example file how to use the generator is provided:

- *example/test.csv*: File which contains an observed streamflow time series with daily values used as input for the streamflow generator. May also be used what your input file has to look like. Otherwise the application will not be able to read the file. A blueprint is depicted in the table below.

> **Date**           | **Q**
> -------------------|------
> dd-mm-YYYY       | x1
> dd-mm-YYYY       | x2
> dd-mm-YYYY       | x3
>  ...                        | ...

 [^1]: Salas, J.D., 1993. Analysis and modeling of hydrologic time series. In: Maidment, D.R. (Ed.), Handbook of Hydrology. McGraw-Hill, pp. 1-72.
