<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Last-changedate](https://img.shields.io/badge/last%20change-2017--05--24-brightgreen.svg)](https://github.com/bocinsky/asian_niche/commits/master) [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.0-brightgreen.svg)](https://cran.r-project.org/) [![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/) [![DOI](https://zenodo.org/badge/52899692.svg)](https://zenodo.org/badge/latestdoi/52899692)

Research compendium for d'Alpoim Guedes and Bocinsky *in review*
----------------------------------------------------------------

**When using the code included in this research compendium, please cite *all* of the following:**

d'Alpoim Guedes, Jade and R. Kyle Bocinsky. Climate change stimulated agricultural innovation and exchange across Asia. In review.

d'Alpoim Guedes, Jade and R. Kyle Bocinsky. Research compendium for: *Climate change stimulated agricultural innovation and exchange across Asia*, in review. Version 0.9.0. Zenodo. <http://doi.org/10.5281/zenodo.583157>

d'Alpoim Guedes, Jade and R. Kyle Bocinsky. Data output for: *Climate change stimulated agricultural innovation and exchange across Asia*, in review. Version 0.9.0. Zenodo. <http://doi.org/10.5281/zenodo.583154>

### Compendium DOI:

[![DOI](https://zenodo.org/badge/52899692.svg)](https://zenodo.org/badge/latestdoi/52899692)

The files at the URL above will generate the results as found in the publication. The files hosted at <https://github.com/bocinsky/asian_niche> are the development versions and may have changed since this compendium was released.

### Authors of this repository:

-   R. Kyle Bocinsky (<bocinsky@gmail.com>)
-   Jade d'Alpoim Guedes (<jadeguedes@gmail.com>)

### Overview of contents

This repository is a research compendium for d'Alpoim Guedes and Bocinsky (in review). The compendium contains all code associated with the analyses described and presented in the publication, as well as a Docker environment (described in the `Dockerfile`) for running the code. The primary files contained in this reposity are:

-   **`README.Rmd`**: An RMarkdown source file, and its knitted Markdown output `README.md`
-   **`asian_niche.Rproj`**: An RStudio project configuration file for this study
-   **`asian_niche.R`**: The R script that runs the entire analysis, and is designed to be called from the command prompt
-   **`asian_niche.sh`**: A Bash script for building the Docker container (see below), running `asian_niche.R`, and compressing the output for uploading to Zenodo
-   **`Dockerfile`**: A text document that contains the instructions for building the Docker container
-   **`LICENSE`**: A text file of the MIT license for this repository and the code within it
-   **`src/`**: A directory containing auxiliary R code for functions used in the `asian_niche.R` script
-   **`DATA/`**: A directory containing input data for the analysis not downloaded from elsewhere. This includes `crops.csv`, a table of cultivars and their estimated thermal requirements.

### The research compendium

To download this research compendium as you see it on GitHub, for offline browsing, [install git on your computer](https://git-scm.com/) and use this line at a Bash prompt ("Terminal" on macOS and Unix-alikes, "Command Prompt" on Windows):

``` bash
git clone https://github.com/bocinsky/asian_niche.git
```

**All commands described in this research compendium are meant to be run from a Bash prompt.**

### Archaeological site data

Archaeological site location data are sensitive information due to the possibility of looting, and archaeological ethics require that we restrict access to those data. Accordingly, an essential component of this analysis is **not** shipped in this open GitHub repository or archived with Zenodo. We have instead archived the site location data necessary to run this analysis with the [Digital Archaeological Record (tDAR)](https://www.tdar.org/) under restricted access. Users who want to run this analysis need to request access through tDAR, which we will provide to any researcher with a reasonable affiliation (academic or otherwise). The main purpose is to track to whom we provide access. Please [contact the authors](mailto:bocinsky@gmail.com,jadeguedes@gmail.com) if you wish to use these data in your research or share them with others.

The data are available through tDAR at the following DOI: [10.6067/XCV8MK6G05](https://doi.org/10.6067/XCV8MK6G05). Please go to the site and select "Request Access, Submit Correction, Comment" under the downloads section in the right panel. You will have to [create a tDAR user account](https://core.tdar.org/account/new) and agree to the [tDAR user agreement](https://www.tdar.org/about/policies/terms-of-use/).

Once you have access, you should be able to download the `DALPOIMGUEDES_BOCINSKY_2017.xlsx` Microsoft Excel file. Place the file in the `asian_niche/DATA/` directory, where you should also find the `crops.csv` file.

### The Docker container

[Docker](https://www.docker.com/) is a virtual computing environment that facilitates reproducible research---it allows for research results to be produced independent of the machine on which they are computed. Docker users describe computing environments in a text format called a "Dockerfile", which when read by the Docker software builds a virtual machine, or "container". Other users can then load the container on their own computers. Users can upload container images to [Docker Hub](https://hub.docker.com/), and the image for this research is available at <https://hub.docker.com/r/bocinsky/asian_niche/>.

We have included a Dockerfile which builds a Docker container for running the analyses described in the paper. It uses [`rocker/geospatial:3.3.3`](https://hub.docker.com/r/rocker/geospatial/), which provides R, [RStudio Server](https://www.rstudio.com/products/rstudio/download-server/), the [tidyverse](http://tidyverse.org/) of R packages as its base image and adds several geospatial software packages ([GDAL](http://www.gdal.org/), [GEOS](https://trac.osgeo.org/geos/), and [proj.4](http://proj4.org/). The Dockerimage (1) adds ffmpeg, (2) updates the R packages to the [2017-05-15 MRAN snapshot](https://mran.microsoft.com/timemachine/), and (3) installs the R software packages required by the script (based on the same MRAN snapshot).

#### Downloading and running the Docker container image

The commands below demonstrate three ways to run the docker container. In each, we use the `-w` argument to set the working directory to `/asian_niche`. See this [Docker cheat sheet](https://github.com/wsargent/docker-cheat-sheet) for other arguments. Using the ":0.9.0" tag will ensure you are running the version of the code that generates the d'Alpoim Guedes and Bocinsky (2017) results---the first time you run the Docker image, it will download it from the Docker Hub.

**Be sure to copy the archaeological site data into the `asian_niche/DATA/` directory prior to running the `asian_niche.R` script!** One way you can do so by adding the archaeological site data to your local `asian_niche/DATA/` directory, commenting the `DALPOIMGUEDES_BOCINSKY_2017.xlsx` entry in the `.dockerignore` file, and rebuilding the docker image using the command below.

##### Run the analysis directly

To run the analyses directly, call the `asian_niche.R` script at the end of the run command:

``` bash
docker run -w /asian_niche bocinsky/asian_niche:0.9.0 Rscript asian_niche.R
```

##### Run the analysis interactively from the terminal

Alternatively, you can run the container in interactive mode and load the script yourself:

``` bash
docker run -w /asian_niche -it bocinsky/asian_niche:0.9.0 bash
```

The `asian_niche.R` script has been designed to be run from the shell using Python-style argument parsing. To run, simply enter `Rscript --vanilla asian_niche.R` at the shell prompt. Passing the `--vanilla` option runs the script in a "fresh" R environment. Run `Rscript asian_niche.R --help` to see all available options.

You can use the `exit` command to stop the container.

##### Run the analysis from within a Dockerized RStudio IDE

Finally, you can host RStudio Server locally to use the RStudio browser-based IDE. Run:

``` bash
docker run -p 8787:8787 bocinsky/asian_niche:0.9.0
```

Then, open a browser (we find [Chrome](https://www.google.com/chrome/) works best) and navigate to "localhost:8787" or or run `docker-machine ip default` in the shell to find the correct IP address, and log in with **rstudio**/**rstudio** as the user name and password. In the explorer (lower right pane in RStudio), navigate to the `asian_niche` directory, and click the `asian_niche.Rproj` to open the project.

#### Building the Docker container from scratch

If you wish to build the Docker container locally for this project from scratch, simply `cd` into the `asian_niche/` directory and run:

``` bash
docker build -t bocinsky/asian_niche .
```

The `-t` argument gives the resulting container image a name. You can then run the container as described above.

### One script to rule them all

Alternatively to everything above, you can run the **`asian_niche.sh`** script, which builds the Docker container, runs the analysis, exports the output for Zenodo archiving, and prepares all of the supplemental information included in the d'Alpoim Guedes and Bocinsky (2017) paper. **This is the master script we ran for d'Alpoim Guedes and Bocinsky (2017).**

Be sure to set the `ARCH_SITES` variable at the top of `asian_niche.sh` to point to the location of the archaeological site data on your local file system.

Run it all using bash:

``` bash
bash asian_niche.sh
```

### A note on run time

This analysis has been designed to take advantage of modern multi-core or multi-CPU computer architectures. By default, it will run on two cores—i.e., sections of the code will run in parallel approximately twice as fast as on a single core. The analysis also consumes quite a bit of memory. On two (relatively high-speed) cores, run-time of the entire analysis is **approximately 12 hours**. This can be sortened dramatically by running with a higher number of cores/processors and amount of memory, if available. During testing, we were able to run the entire analysis in under three hours using a 20 core cluster.

### Output

The GitHub repository for this project does not contain the output generated by the script---3.25 GB of uncompressed data. All output data is available as a separate Zenodo archive at:

<http://doi.org/10.5281/zenodo.583154>

The `OUTPUT/` directory contains all data generated by the `asian_niche.R` script:

-   `session_info.txt` describes the computational infrastructure within which the script was run
-   `packages.bib` is a BibTeX bibliography file including all packages loaded in the `asian_niche.R` script
-   `DATA/` contains data downloaded from web sources for this analysis
-   `MODELS/` contains R data objects describing the Kriging interpolation models across the study area
-   `RECONS/` contains NetCDF format raster bricks of the model output (i.e., the reconstructed crop niches)
-   `SITE_DENSITIES/` contains figures of the estimated chronometric probability density for each site in our database
-   `FIGURES/` contains all figures output by the script, including videos of how each crop niche changes over time
-   `TABLES/` contains tables of the raw site chronometric data without locational information, and the modeled chronometric probability and niche information for each site.

### Licenses

Code: [MIT](http://opensource.org/licenses/MIT) year: 2017<br> Copyright holders: R. Kyle Bocinsky and Jade d'Alpoim Guedes

### Contact

R. Kyle Bocinsky, PhD, RPA<br> Research Associate<br> Crow Canyon Archaeological Center<br> 23390 Road K, Cortez, CO 81321<br> 970.564.4384 – Office<br> 770.362.6659 – Mobile<br> <bocinsky@gmail.com> – Email<br> [bocinsky.io](http://www.bocinsky.io/) – Web

### R package citations

The following R packages were used for this analysis:

Analytics, Revolution, and Steve Weston. 2015a. *DoParallel: Foreach Parallel Adaptor for the ’Parallel’ Package*. <https://CRAN.R-project.org/package=doParallel>.

———. 2015b. *Foreach: Provides Foreach Looping Construct for R*. <https://CRAN.R-project.org/package=foreach>.

Arnold, Jeffrey B. 2017. *Ggthemes: Extra Themes, Scales and Geoms for ’Ggplot2’*. <https://CRAN.R-project.org/package=ggthemes>.

Bache, Stefan Milton, and Hadley Wickham. 2014. *Magrittr: A Forward-Pipe Operator for R*. <https://CRAN.R-project.org/package=magrittr>.

Bengtsson, Henrik. 2016. *R.utils: Various Programming Utilities*. <https://CRAN.R-project.org/package=R.utils>.

Bivand, Roger, and Nicholas Lewin-Koh. 2017. *Maptools: Tools for Reading and Handling Spatial Objects*. <https://CRAN.R-project.org/package=maptools>.

Bivand, Roger, Tim Keitt, and Barry Rowlingson. 2017. *Rgdal: Bindings for the Geospatial Data Abstraction Library*. <https://CRAN.R-project.org/package=rgdal>.

Bocinsky, R. Kyle. 2017. *FedData: Functions to Automate Downloading Geospatial Data Available from Several Federated Data Sources*. <https://CRAN.R-project.org/package=FedData>.

Boettiger, Carl. 2015. *Knitcitations: Citations for ’Knitr’ Markdown Files*. <https://CRAN.R-project.org/package=knitcitations>.

Chamberlain, Scott. 2017. *Rgbif: Interface to the Global ’Biodiversity’ Information Facility ’Api’*. <https://CRAN.R-project.org/package=rgbif>.

Douglas Nychka, Reinhard Furrer, John Paige, and Stephan Sain. 2015. “Fields: Tools for Spatial Data.” Boulder, CO, USA: University Corporation for Atmospheric Research. doi:[10.5065/D6W957CT](https://doi.org/10.5065/D6W957CT).

for R, Doug McIlroy. Packaged, Thomas P Minka, and transition to Plan 9. 2015. *Mapproj: Map Projections*. <https://CRAN.R-project.org/package=mapproj>.

Fraley, Chris, and Adrian E. Raftery. 2002. “Model-Based Clustering, Discriminant Analysis and Density Estimation.” *Journal of the American Statistical Association* 97: 611–31.

Fraley, Chris, Adrian E. Raftery, Thomas Brendan Murphy, and Luca Scrucca. 2012. *Mclust Version 4 for R: Normal Mixture Modeling for Model-Based Clustering, Classification, and Density Estimation*. Technical Report. Department of Statistics, University of Washington.

Francois, Romain. 2014. *Bibtex: Bibtex Parser*. <https://CRAN.R-project.org/package=bibtex>.

Harrell Jr, Frank E, with contributions from Charles Dupont, and many others. 2017. *Hmisc: Harrell Miscellaneous*. <https://CRAN.R-project.org/package=Hmisc>.

Hijmans, Robert J. 2016. *Raster: Geographic Data Analysis and Modeling*. <https://CRAN.R-project.org/package=raster>.

Lees, Jonathan M. 2012. *Geomapdata: Data for Topographic and Geologic Mapping*. <https://CRAN.R-project.org/package=geomapdata>.

Neuwirth, Erich. 2014. *RColorBrewer: ColorBrewer Palettes*. <https://CRAN.R-project.org/package=RColorBrewer>.

Parnell, Andrew. 2016. *Bchron: Radiocarbon Dating, Age-Depth Modelling, Relative Sea Level Rate Estimation, and Non-Parametric Phase Modelling*. <https://CRAN.R-project.org/package=Bchron>.

Pebesma, Edzer. 2017. *Sf: Simple Features for R*. <https://CRAN.R-project.org/package=sf>.

Pierce, David. 2017. *Ncdf4: Interface to Unidata netCDF (Version 4 or Earlier) Format Data Files*. <https://CRAN.R-project.org/package=ncdf4>.

Plate, Tony, and Richard Heiberger. 2016. *Abind: Combine Multidimensional Arrays*. <https://CRAN.R-project.org/package=abind>.

Sievert, Carson, Chris Parmer, Toby Hocking, Scott Chamberlain, Karthik Ram, Marianne Corvellec, and Pedro Despouy. 2017. *Plotly: Create Interactive Web Graphics via ’Plotly.js’*. <https://CRAN.R-project.org/package=plotly>.

Vaidyanathan, Ramnath, Yihui Xie, JJ Allaire, Joe Cheng, and Kenton Russell. 2016. *Htmlwidgets: HTML Widgets for R*. <https://CRAN.R-project.org/package=htmlwidgets>.

Wickham, Hadley. 2017a. *Purrrlyr: Tools at the Intersection of ’Purrr’ and ’Dplyr’*. <https://CRAN.R-project.org/package=purrrlyr>.

———. 2017b. *Tidyverse: Easily Install and Load ’Tidyverse’ Packages*. <https://CRAN.R-project.org/package=tidyverse>.

Wood, S. N. 2003. “Thin-Plate Regression Splines.” *Journal of the Royal Statistical Society (B)* 65 (1): 95–114.

———. 2004. “Stable and Efficient Multiple Smoothing Parameter Estimation for Generalized Additive Models.” *Journal of the American Statistical Association* 99 (467): 673–86.

———. 2011. “Fast Stable Restricted Maximum Likelihood and Marginal Likelihood Estimation of Semiparametric Generalized Linear Models.” *Journal of the Royal Statistical Society (B)* 73 (1): 3–36.

———. 2016. “Just Another Gibbs Additive Modeler: Interfacing JAGS and mgcv.” *Journal of Statistical Software* 75 (7): 1–15. doi:[10.18637/jss.v075.i07](https://doi.org/10.18637/jss.v075.i07).

Wood, S.N. 2006. *Generalized Additive Models: An Introduction with R*. Chapman; Hall/CRC.

Zeileis, Achim, and Gabor Grothendieck. 2005. “Zoo: S3 Infrastructure for Regular and Irregular Time Series.” *Journal of Statistical Software* 14 (6): 1–27. doi:[10.18637/jss.v014.i06](https://doi.org/10.18637/jss.v014.i06).
