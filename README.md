fusionModel
================
Kevin Ummel (<ummel@berkeley.edu>)

-   [Overview](#overview)
-   [Motivation](#motivation)
-   [Methodology](#methodology)
-   [Installation](#installation)
-   [Quick start example](#quick-start-example)
-   [Advanced examples](#advanced-examples)
-   [Analysis](#analysis)
-   [Validation](#validation)

# Overview

**fusionModel** enables variables unique to a “donor” dataset to be
statistically simulated for (i.e. *fused to*) a “recipient” dataset.
Variables common to both the donor and recipient are used to model and
simulate the fused variables. The package provides a simple and
efficient interface for general data fusion in *R*, leveraging
state-of-the-art machine learning algorithms from Microsoft’s
[LightGBM](https://lightgbm.readthedocs.io/en/latest/) framework. It
also provides tools for analyzing synthetic/simulated data, calculating
uncertainty, and validating fusion output.

fusionModel was developed to allow statistical integration of microdata
from disparate social surveys. It is the data fusion workhorse
underpinning the larger fusionACS data platform under development at the
[Socio-Spatial Climate
Collaborative](https://sc2.berkeley.edu/fusionacs-people/). In this
context, fusionModel is used to fuse variables from a range of social
surveys onto microdata from the American Community Survey, allowing for
analysis and spatial resolution otherwise impossible.

# Motivation

The desire to “fuse” or otherwise integrate independent datasets has a
long history, dating to at least the early 1970’s ([Ruggles and Ruggles
1974](https://www.nber.org/system/files/chapters/c10115/c10115.pdf);
[Alter
1974](https://www.nber.org/system/files/chapters/c10116/c10116.pdf)).
Social scientists have long recognized that large amounts of unconnected
data are “out there” – usually concerning the characteristics of
households and individuals (i.e. microdata) – which we would, ideally,
like to integrate and analyze as a whole. This process is also generally
known as “Statistical Data Integration” (SDI).

The most prominent examples of data fusion have involved administrative
record linkage. This consists of exact matching or probabilistic linking
of independent datasets, using observable information like social
security numbers, names, or birth dates of individuals. Record linkage
is the gold standard and can yield incredibly important insights and
high levels of statistical confidence, as evidenced by the pioneering
work of Raj Chetty and colleagues.

However, record linkage is rarely feasible for the kinds of microdata
that most researchers use day-to-day (nevermind the difficulty of
accessing administrative data). While the explosion of online tracking
and social network data will undoubtedly offer new lines of analysis,
for the time being, at least, social survey microdata remain
indispensable. The challenge and promise recognized 50 years ago by
Nancy and Richard Ruggles remains true today:

> Unfortunately, no single microdata set contains all of the different
> kinds of information required for the problems which the economist
> wishes to analyze. Different microdata sets contain different kinds of
> information…A great deal of information is collected on a sample
> basis. Where two samples are involved the probability of the same
> individual appearing in both may be very small, so that exact matching
> is impossible. Other methods of combining the types of information
> contained in the two different samples into one microdata set will be
> required. (Ruggles and Ruggles 1974; 353-354)

Practitioners regularly impute or otherwise predict a variable or two
from one dataset on to another. Piecemeal, *ad hoc* data fusion is a
common necessity of quantitative research. Proper data fusion, on the
other hand, seeks to systematically combine “two different samples into
one microdata set”.

Say we have microdata from two independent household surveys, A and B.
We specify that A is the “recipient” dataset and B is the “donor”.
Conceptually, our goal is to generate a new dataset, C, that has “the
best of both worlds”: the original survey responses of A plus a
realistic representation of how each respondent in A might have answered
the questionnaire of survey B. To help do this, we identify a set of
common variables, X, that both surveys solicit; in practice, these are
often things like household size, income, or respondent age. We then
attempt to fuse the variables unique to B (the “fusion variables”) to
the original microdata of A, conditional on X.

Note that this requires the simulation of outcomes across multiple
variables, with realistic variance, and a realistic representation of
the (typically unknown) joint distribution among those variables. No
easy task!

# Methodology

The fusionModel package was designed with the general needs of the
fusionACS project in mind. The goal was to create a data fusion tool
that meets the following requirements:

1.  Handle both categorical and continuous variables
2.  Accommodate semi-continuous (zero-inflated) variables
3.  Fuse variables “one-by-one” and/or in “blocks”
4.  Ensure realistic values for fused variables
5.  Scale efficiently for larger datasets
6.  Employ a data modeling approach that:
    -   Makes no distributional assumptions (i.e. non-parametric)
    -   Automatically detects non-linear and interaction effects
    -   Automatically selects predictor variables
    -   Avoids overfitting to noise
7.  Determine suitable fusion order/sequence
8.  Calculate uncertainty of analyses using fused data
9.  Assess the overall quality of fused data

In an effort to meet requirements (1-6), fusionModel combines gradient
boosting machine learning with approximate nearest neighbor matching.
This falls under the heading of “mixed methods” in statistical matching
(SM); see Section 3.1.3 in [Lewaa et
al. (2021)](https://content.iospress.com/articles/statistical-journal-of-the-iaos/sji210835)
for a review. Mixed methods, including that employed by fusionModel,
generally have the following structure:

-   Step 1: A statistical model is fit to donor data to predict fusion
    variable
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    conditional on shared variables
    ![X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X "X").
-   Step 2: A prediction of
    ![Y\|X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%7CX "Y|X")
    is made for each record in the donor
    (![\\hat{Y}\_{d}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Bd%7D "\hat{Y}_{d}"))
    and recipient
    (![\\hat{Y}\_{r}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Br%7D "\hat{Y}_{r}")).
-   Step 3: A “real” donor value from
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    is selected for each recipient record, based on the similarity of
    ![\\hat{Y}\_{d}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Bd%7D "\hat{Y}_{d}")
    and
    ![\\hat{Y}\_{r}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Br%7D "\hat{Y}_{r}").

fusionModel’s approach can be viewed as a logical extension of existing
mixed SM techniques that utilize either stochastic or predictive mean
matching in the modeling step. It builds upon them in the following
ways:

-   Microsoft’s [LightGBM](https://lightgbm.readthedocs.io/en/latest/)
    gradient boosting framework is used for Step 1. This allows the
    extensive and flexible modeling requirements specified in (6) to be
    met.
-   When
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    is continuous, the conditional mean is modeled as well as
    conditional quantiles, providing a non-parametric description of the
    conditional distribution. Step 3 then pairs records with similar
    conditional distributions.
-   When
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    is categorical, LightGBM produces conditional class probabilities
    that are easily sampled to create the fused values (i.e. Step 3 can
    be ignored).
-   When
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    is a “block” of variables (possibly a mix of categorical and
    continuous), Step 3 utilizes all of the conditional distribution
    and/or class probability information to pair similar records.
-   When
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    is semi-continuous (zero-inflated), a LightGBM classification model
    first simulates zero values. Non-zero values are then fused as
    above, conditional on
    ![Y\\neq0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%5Cneq0 "Y\neq0").

In short, the novel aspects of fusionModel are its use of
state-of-the-art gradient boosting techniques within a SM framework
*and* its use of predicted quantiles in the continuous case to describe
each record’s conditional distribution for the purposes of record
matching. Together, they obviate the need for parametric assumptions,
while automating selection of variables and detection of non-linear and
interaction effects.

The training of LightGBM models in Step 1 uses *k*-fold cross-validation
to help select a suitable set of hyperparameters to “tune” the model for
out-of-sample prediction (i.e. reduce overfitting). Since LightGBM uses
[OpenMP](https://en.wikipedia.org/wiki/OpenMP) multithreading “under the
hood”, this key piece of the processing chain can easily make use of
multiple computing cores.

Similarly, data manipulation tasks (via
[data.table](https://github.com/Rdatatable/data.table)) and input-output
operations (via [fst](https://github.com/fstpackage/fst)) are also
OpenMP-enabled for speed.

fusionModel uses the [ANN library](https://www.cs.umd.edu/~mount/ANN/)
(implemented via [RANN](https://github.com/jefferislab/RANN)) for an
approximate nearest neighbor hot deck in Step 3. This includes the
ability to sample donor records from a fixed k nearest matches or,
alternatively, to sample all records within a specified distance.
Euclidean distance is calculated using all of the variables from Step 2
describing the conditional distributions, after applying a scaling
procedure to ensure all variables are of similar magnitude.

*Under development*: The package includes an experimental chain()
function to determine a pseudo-optimal sequencing (or “chaining”) of
fusion variables. Since fusion variables early in the sequence become
available as predictors for those later on, the sequence clearly matters
to the quality of the output. But there is no consensus in the
literature on how to go about it. The nascent chain() function uses
fast-fitting, linear LASSO models to first fit a model for each Y where
all other Y’s are available as predictors. It then fits a model for each
Y using only X predictors. The R-squared of the latter is compared to
the former, and the variable with the highest ratio is selected as the
initial fusion variable. This process proceeds greedily until the full
sequence is constructed.

TO DO: analyze()…

TO DO: compare()…

# Installation

``` r
devtools::install_github("ummel/fusionModel")
library(fusionModel)
```

# Quick start example

The package includes example microdata from the 2015 Residential Energy
Consumption Survey (see `?recs` for details). For real-world use cases,
the donor and recipient data are typically independent. For illustrative
purposes, we will use the `recs` dataset to create both our “donor” and
“recipient”.

``` r
# Create donor and recipient datasets
donor <- select(recs, 2:11, electricity, natural_gas, aircon)
recipient <- select(recs, 2:11)

# Specify fusion and shared/common predictor variables
fusion.vars <- c("electricity", "natural_gas", "aircon")
predictor.vars <- names(recipient)

# Build the fusion model (see ?train)
model <- train(data = donor, y = fusion.vars, x = predictor.vars)
```

    3 fusion variables
    10 initial predictor variables
    5686 observations
    Building LightGBM models...
    Fusion model saved to: /home/kevin/Documents/Projects/fusionModel/fusion_model.fsn

By default, `train()` writes a fusion model object (“fusion_model.fsn”)
to the current working directory and returns the file path. We can use
it to fuse variables to the recipient:

``` r
# Fuse 'fusion.vars' to the recipient (see ?fuse)
sim <- fuse(data = recipient, file = model)
```

    Fusing donor variables to recipient...

``` r
head(sim)
```

       electricity natural_gas                                     aircon
    1:        5270       300.0            Central air conditioning system
    2:       12170         0.0            Central air conditioning system
    3:       19660       294.4 Both a central system and individual units
    4:        9850       992.0            Central air conditioning system
    5:        7900       395.0   Individual window/wall or portable units
    6:        2777       315.4   Individual window/wall or portable units

# Advanced examples

# Analysis

# Validation
