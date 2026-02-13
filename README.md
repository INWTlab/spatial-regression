## PoC: Extension of Multivariate Kernel Density Estimation with Anonymized Geo-coordinates

* **Key Reference Paper:** [https://academic.oup.com/jrsssa/article/180/1/161/7068203](https://academic.oup.com/jrsssa/article/180/1/161/7068203).
* **The Extension:** Incorporating a regression model into the density estimation of coordinates.
* **Comparison Method:** Implemented as the `dbivr` function in the **KernelHeaping** package: [https://cran.r-project.org/web/packages/Kernelheaping/index.html](https://cran.r-project.org/web/packages/Kernelheaping/index.html).

---

### Usage Instructions
* The script `evaluation.R` (found in the `inst/scripts` folder) contains the code to run simulation studies comparing the extension with the `dbivr` method.
* The script utilizes the following functionalities, which are provided as functions in separate scripts and called via the `source` command:
    * **Generation** of simulated geo-coordinate datasets.
    * **Estimation** of coordinate density using both the extension and the `dbivr` method.
    * Parameterized function for conducting simulation studies and calculating **comparison** metrics.
