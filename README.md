# ![](images/icon_48x48.png "Icon") Streamflow Generator (Qsynth)

**Streamflow generator** synthesises an artificial streamflow time series with daily values based on an autoregressive model. The generation merely based on an observed streamflow time series. The implemented approach takes roughly a bearing described in Salas (1993)[^1]

## License

This software can be distributed freely under the GPL v2 license. Please read the LICENSE for further information.

© 2017, Robin Schwemmle (<rschwemmle@yahoo.de>)

## Structure

* Directions for Use
* Example

## Directions for Use

![](images/GUI.png "GUI")



> **Input**

>`dd/MM/yyyy`: Indicate here the date format of your time series.

>`NA`: Indicate here the format of the missing values in your time series.

>`Load .csv`: Imports with observed streamflow time series into the application. A new window pops up where the path to the .csv-file can be either entered manually or selecting it by browsing.
>![](images/loadcsv.png "Load .csv")
>`Daily Streamflow`: Displays observed streamflow time series.

>`Seasonality`: Daily and monthly runoff regime as well as parde coefficient are illustrated.

---

>**Time Period**

>`Start Date`: Indicate here the start date of the synthetic streamflow time series which will be generated.

>`End Date`: Indicate here the end date of the synthetic streamflow time series which will be generated.

---

>**Generate Streamflow**

>`Run`: Generates the synthetic streamflow time series.

---

>**AR(p)**
>
> `p`: Autoregressive order used to generate the streamflow time series. Value gets updated after streamflow has been succesfully generated.

---

>**Output**

>`Synthetic Streamflow`: Displays generated streamflow time series with and without observed values.

>`Filter`: Displays range of moving average with moving window size.

>`Autocorrelation`: Displays autocorrelation for generated synthetic and observed streamflow time series entirely, for each month and for the residuals of the autoregressive model.

>`Histogram`: Displays histogram for generated synthetic and observed streamflow time series entirely. Also shown the histogram the oserved time series after transformation to normal, after transformation to normal and standardization and the raw output of the autoregressive model. Additionally histogram of the residuals is plotted.

>`IHA`: Displays the IHA indicators of Group 1, Group 2, Group 4,Group 5-1 and Group 5-2 for the synthetic and observed streamflow.

>`Simple Test Statistic`: Displays simple test statisic which compares the minimum, maximum, mean, standard deviation and skewness coefficient between synthetic and observed time series.

>`Volume`: Displays the annually and monthly cumulated volume.

>`Export`: Exporting all the plots as .pdf and generated synthetic streamflow time series to .csv. A new window pops up where the path to the folder can be either entered manually or selecting it by browsing. All files will be exported to the indicated folder.
>![](images/export.png "Export")


---
>**Menu Bar**

> Tools to navigate inside the plots.

>![](images/pan.png "Pan"): Pan

>![](images/zoomin.png "Zoom in"): Zoom in

>![](images/zoomout.png "Zoom out"): Zoom out

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

## Example

An example file how to use the generator is provided:

- *example/example.csv*: File which contains an observed streamflow time series with daily values used as input for the streamflow generator. May also be used what your input file has to look like. Otherwise the application will not be able to read the file. A blueprint is depicted in the table below.

> **Date**           | **Q**
> -------------------|------
> dd-mm-YYYY       | x1
> dd-mm-YYYY       | x2
> dd-mm-YYYY       | x3
>  ...                        | ...

 [^1]: Salas, J.D., 1993. Analysis and modeling of hydrologic time series. In: Maidment, D.R. (Ed.), Handbook of Hydrology. McGraw-Hill, pp. 1-72.
