# LB\_MissPLB — Reproducibility materials for "Scalable Logistic Biplots for Large Binary Matrices via Projected-Data Coordinate Descent"

Repository: [https://github.com/jgbabativam/LB_MissPLB](https://github.com/jgbabativam/LB_MissPLB)

This repository reproduces the data analysis, figures, and computational
benchmark reported in:

> Babativa-Marquez, J.G. \& Vicente-Villardon, J.L. "Scalable Logistic
> Biplots for Large Binary Matrices via Projected-Data Coordinate Descent."


The algorithm proposed in the paper ("PDLB": block coordinate descent with
data projection) is implemented in the CRAN package
[`BiplotML`](https://cran.r-project.org/package=BiplotML)
([source](https://github.com/jgbabativam/BiplotML)). Everything here uses
the published CRAN release of that package — no unpublished or development
code is required.

## Getting started

```sh
git clone https://github.com/jgbabativam/LB_MissPLB.git
cd LB_MissPLB
```

## Requirements

* R >= 4.1.0
* Packages: `BiplotML`, `dplyr`, `ggplot2`, `vegan` (install with
`install.packages(c("BiplotML", "dplyr", "ggplot2", "vegan"))`)

Exact package versions used to produce the results shipped in `output/` are
recorded in [`sessionInfo.txt`](sessionInfo.txt).

## License

Code and derived data in this repository are released under the MIT license
(see [LICENSE](LICENSE)), matching the license of the `BiplotML` package.
The underlying event data are publicly available from the OMC as described
above.

## Authors

* Jose Giovany Babativa-Marquez ([ORCID 0000-0002-4989-7459](https://orcid.org/0000-0002-4989-7459)) — jgbabativam@unal.edu.co
* Jose Luis Vicente-Villardon ([ORCID 0000-0001-7061-5271](https://orcid.org/0000-0001-7061-5271)) — villardon@usal.es

