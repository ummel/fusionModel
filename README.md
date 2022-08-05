fusionModel
================
Kevin Ummel (<ummel@berkeley.edu>)

-   <a href="#overview" id="toc-overview">Overview</a>
-   <a href="#motivation" id="toc-motivation">Motivation</a>
-   <a href="#methodology" id="toc-methodology">Methodology</a>
-   <a href="#installation" id="toc-installation">Installation</a>
-   <a href="#simple-fusion" id="toc-simple-fusion">Simple fusion</a>
-   <a href="#advanced-fusion" id="toc-advanced-fusion">Advanced fusion</a>
-   <a href="#analyzing-fused-data" id="toc-analyzing-fused-data">Analyzing
    fused data</a>

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
like to integrate and analyze as a whole. This aim falls under the
general heading of “Statistical Data Integration” (SDI) ([Lewaa et
al. 2021](https://content.iospress.com/articles/statistical-journal-of-the-iaos/sji210835)).

The most prominent examples of data fusion have involved administrative
record linkage. This consists of exact matching or probabilistic linking
of independent datasets, using observable information like social
security numbers, names, or birth dates of individuals. Record linkage
is the gold standard and can yield incredibly important insights and
high levels of statistical confidence, as evidenced by the pioneering
work of [Raj Chetty](https://opportunityinsights.org/team/raj-chetty/)
and colleagues.

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

The size and nature of the samples involved and the intended analyses
strongly influence the choice of data integration technique and the
structure of the output. This has led to the relevant literature being
both diverse and convoluted, as practitioners take on different data
“setups” and objectives. In the context of fusionACS, we are interested
in the following problem:

We have microdata from two independent surveys, A and B, that sample the
same underlying population and time period (e.g. occupied U.S.
households nationwide in 2018). We specify that A is the “recipient”
dataset and B is the “donor”. Survey A is the American Community Survey
and invariably has a larger sample size than B. The goal is to generate
a new dataset, C, that has the original survey responses of A plus a
realistic representation of how each respondent in A might have answered
the questionnaire of survey B. To do this, we identify a set of common,
harmonized variables X that both surveys solicit; in practice, these are
often things like household size, income, or respondent age. We then
attempt to fuse a set of variables unique to B – call them Z, the
“fusion variables” – onto the original microdata of A, conditional on X.

Note that this requires the simulation of outcomes across multiple
variables, with realistic variance, and a realistic representation of
the (typically unknown) joint distribution among those variables. No
easy task!

# Methodology

The fusionModel package borrows and expands upon ideas from the
statistical matching ([D’Orazio et
al. 2006](https://onlinelibrary.wiley.com/doi/book/10.1002/0470023554)),
imputation ([Little and Rubin
2019](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119482260)),
and data synthesis ([Drechsler
2011](https://link.springer.com/book/10.1007/978-1-4614-0326-5))
literatures to create a flexible, hybrid data fusion tool. The core
approach is an extension of mixed-method predictive mean matching to
take advantage of high-performance gradient boosting algorithms. The
process is sequential, fusing variables individually or in blocks, in
order to accomodate donor datasets that are small relative to the
recipient.

The fusionModel package was designed with the general needs of the
fusionACS project in mind. The goal was to create a data fusion tool
that meets the following requirements:

1.  Accommodate donor and recipient datasets with divergent sample size
2.  Handle continuous, categorical, and semi-continuous (zero-inflated)
    variable types
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
This falls under the heading of “mixed methods” in statistical matching;
see Section 3.1.3 in [Lewaa et
al. (2021)](https://content.iospress.com/articles/statistical-journal-of-the-iaos/sji210835)
for a review. Mixed methods, including that employed by fusionModel,
generally have the following structure:

1.  A statistical model is fit to donor data to predict fusion variables
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    conditional on shared variables
    ![X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X "X").
2.  A prediction (possibly stochastic) of
    ![Y\|X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%7CX "Y|X")
    is made for each record in the donor
    (![\hat{Y}\_{d}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Bd%7D "\hat{Y}_{d}"))
    and recipient
    (![\hat{Y}\_{r}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Br%7D "\hat{Y}_{r}")).
3.  A “real” donor value from
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    is selected for each recipient record, based on the similarity of
    ![\hat{Y}\_{d}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Bd%7D "\hat{Y}_{d}")
    and
    ![\hat{Y}\_{r}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Br%7D "\hat{Y}_{r}").

fusionModel’s approach can be viewed as a logical extension of existing
mixed statistical matching techniques that utilize conditional mean
models in Step 1. It builds upon them in the following ways:

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
    ![Y\neq0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%5Cneq0 "Y\neq0").

In short, the novel aspects of fusionModel are its use of
state-of-the-art gradient boosting techniques within a statistical
matching framework and its use of predicted quantiles in the continuous
case to describe each record’s conditional distribution for the purposes
of record matching. Together, they obviate the need for parametric
assumptions, while automating selection of variables and detection of
non-linear and interaction effects.

The training of LightGBM models in Step 1 uses *k*-fold cross-validation
to help select a suitable set of hyperparameters to “tune” the model for
out-of-sample prediction (i.e. reduce overfitting). Since LightGBM uses
[OpenMP](https://en.wikipedia.org/wiki/OpenMP) multithreading “under the
hood”, this key piece of the processing chain easily makes use of
multiple computing cores. Similarly, data manipulation tasks (via
[data.table](https://github.com/Rdatatable/data.table)) and input-output
operations (via [fst](https://github.com/fstpackage/fst)) are also
OpenMP-enabled for speed.

fusionModel uses the [ANN](https://www.cs.umd.edu/~mount/ANN/) library
(implemented via [RANN](https://github.com/jefferislab/RANN)) for an
approximate nearest neighbor hot deck in Step 3. This includes the
ability to sample donor records from a fixed *k* nearest matches or,
alternatively, to sample all records within a specified distance.
Euclidean distance is calculated using all of the variables from Step 2
describing the conditional distributions, after applying a scaling
procedure to ensure all variables are of similar magnitude.

The package includes an experimental `blockchain()` function to
determine a plausible fusion sequence (i.e. chain) and, optionally, a
strategy for assigning variables to blocks. It is not always obvious
when or why or which variables to fuse as blocks. Even in the absence of
blocks (i.e. one-by-one fusion), the literature provides little guidance
on how best to sequence the fusion process. `blockchain()` attempts to
provide a data-driven answer to these questions. See the [Advanced
fusion](#advanced-fusion) section below for more information.

# Installation

``` r
devtools::install_github("ummel/fusionModel")
```

``` r
library(fusionModel)
```

    fusionModel v2.1.0 | https://github.com/ummel/fusionModel

# Simple fusion

The package includes example microdata from the [2015 Residential Energy
Consumption
Survey](https://www.eia.gov/consumption/residential/data/2015/) (see
`?recs` for details). For real-world use cases, the donor and recipient
data are typically independent and vary in sample size. For illustrative
purposes, we will randomly split the `recs` microdata into separate
“donor” and “recipient” datasets with an equal number of observations.

``` r
# Rows to use for donor dataset
d <- seq(from = 1, to = nrow(recs), by = 2)

# Create donor and recipient datasets
donor <- recs[d, c(2:16, 20:22)]
recipient <- recs[-d, 2:14]

# Specify fusion and shared/common predictor variables
predictor.vars <- names(recipient)
fusion.vars <- setdiff(names(donor), predictor.vars)
```

The `recipient` dataset contains 13 variables that are shared with
`donor`. These shared “predictor” variables provide a statistical link
between the two datasets. fusionModel exploits the information in these
shared variables.

``` r
predictor.vars
```

     [1] "income"      "age"         "race"        "education"   "employment" 
     [6] "hh_size"     "division"    "urban_rural" "climate"     "renter"     
    [11] "home_type"   "year_built"  "heat_type"  

There are 5 “fusion variables” unique to `donor`. These are the
variables that will be fused to `recipient`. This includes a mix of
continuous and categorical (factor) variables.

``` r
# The variables to be fused
sapply(donor[fusion.vars], class)
```

    $insulation
    [1] "ordered" "factor" 

    $aircon
    [1] "factor"

    $square_feet
    [1] "integer"

    $electricity
    [1] "numeric"

    $natural_gas
    [1] "numeric"

We create a fusion model using the `train()` function. The minimal usage
is shown below. See `?train` for additional function arguments and
options. This results in a “.fsn” (fusion) object being saved to
“test_model.fsn” in the current working directory.

``` r
# Train a fusion model
fsn.model <- train(data = donor, 
                   y = fusion.vars, 
                   x = predictor.vars,
                   file = "test_model.fsn")
```

    5 fusion variables
    13 initial predictor variables
    2843 observations
    Training step 1 of 5: insulation
    Training step 2 of 5: aircon
    Training step 3 of 5: square_feet
    Training step 4 of 5: electricity
    Training step 5 of 5: natural_gas
    Fusion model saved to: /home/kevin/Documents/Projects/fusionModel/test_model.fsn 
    Total processing time: 5 secs 

To fuse variables to `recipient`, we simply pass the recipient data and
path of the .fsn model to the `fuse()` function. Each variable specified
in `fusion.vars` is fused in the order provided. The default setting
`k = 5` means that the five nearest neighbors in the donor are
identified for each recipient observation, and one of them is randomly
selected to provide the fused value.

``` r
# Fuse 'fusion.vars' to the recipient
sim <- fuse(data = recipient, 
            file = fsn.model,
            k = 5)
```

    5 fusion variables
    13 initial predictor variables
    2843 observations
    Generating 1 implicate 
    Fusion step 1 of 5: insulation
    -- Predicting LightGBM models
    -- Simulating fused values
    Fusion step 2 of 5: aircon
    -- Predicting LightGBM models
    -- Simulating fused values
    Fusion step 3 of 5: square_feet
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Fusion step 4 of 5: electricity
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Fusion step 5 of 5: natural_gas
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Total processing time: 3.98 secs 

Notice that the status messages for the first two fusion variables –
“insulation” and “aircon” – do not include the nearest neighbors step.
This is because they are categorical variables being fused on their own.
In this case, fusionModel uses the LightGBM model’s predicted class
probabilities to simulate the fused value, rather than *k*-nearest
neighbors (as is done for the continuous variables).

Let’s look at the the recipient dataset’s fused/simulated variables.
Note that your results will look different, because each call to
`fuse()` generates a unique, probabilistic set of outcomes.

``` r
head(sim)
```

                 insulation                                   aircon square_feet
    1:       Well insulated          Central air conditioning system        1848
    2: Adequately insulated Individual window/wall or portable units        2100
    3:       Well insulated Individual window/wall or portable units         481
    4: Adequately insulated          Central air conditioning system        2215
    5:       Well insulated          Central air conditioning system        3146
    6:     Poorly insulated Individual window/wall or portable units         840
       electricity natural_gas
    1:       11330           0
    2:       14180           0
    3:        2396           0
    4:        7440        1152
    5:       19100           0
    6:       10360         172

We can do some quick sanity checks to compare the distribution of the
fusion variables in `donor` with those in `sim`. This, at least,
confirms that the fusion output is not obviously wrong.

``` r
sim <- data.frame(sim)

# Compare means of the continuous variables
cbind(donor = colMeans(donor[fusion.vars[3:5]]), sim = colMeans(sim[fusion.vars[3:5]]))
```

                    donor        sim
    square_feet  2070.784  2063.4478
    electricity 10994.517 11067.6012
    natural_gas   338.154   326.6724

``` r
# Compare frequencies of categorical variable classes
cbind(donor = table(donor$insulation), sim = table(sim$insulation))
```

                         donor  sim
    Not insulated           40   30
    Poorly insulated       459  471
    Adequately insulated  1401 1404
    Well insulated         943  938

``` r
cbind(donor = table(donor$aircon), sim = table(sim$aircon))
```

                                               donor  sim
    Central air conditioning system             1788 1822
    Individual window/wall or portable units     545  519
    Both a central system and individual units   125  112
    No air conditioning                          385  390

And we can look at kernel density plots of the non-zero values for the
continuous variables to see if the univariate distributions in `donor`
are generally similar in `sim`.

![](man/figures/README-unnamed-chunk-10-1.png)<!-- -->

Since each call to `fuse()` returns a single probabilistic version of
the fusion variables, we generally want to create multiple versions –
called *implicates* – in order to reduce bias in point estimates and
calculate associated uncertainty. We can do this using the `M` argument
within `fuse()`. Here we generate 10 implicates.

``` r
# Fuse multiple implicates to the recipient
sim10 <- fuse(data = recipient, 
              file = fsn.model,
              k = 5,
              M = 10)
```

    5 fusion variables
    13 initial predictor variables
    2843 observations
    Generating 10 implicates 
    Fusion step 1 of 5: insulation
    -- Predicting LightGBM models
    -- Simulating fused values
    Fusion step 2 of 5: aircon
    -- Predicting LightGBM models
    -- Simulating fused values
    Fusion step 3 of 5: square_feet
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Fusion step 4 of 5: electricity
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Fusion step 5 of 5: natural_gas
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Total processing time: 5.59 secs 

Note that each implicate in `sim10` is identified by the “M”
variable/column.

``` r
head(sim10)
```

       M           insulation                                   aircon square_feet
    1: 1       Well insulated          Central air conditioning system        1848
    2: 1 Adequately insulated          Central air conditioning system        2697
    3: 1     Poorly insulated Individual window/wall or portable units         552
    4: 1     Poorly insulated          Central air conditioning system        3300
    5: 1       Well insulated          Central air conditioning system        1202
    6: 1     Poorly insulated          Central air conditioning system        1066
       electricity natural_gas
    1:       15330           0
    2:        6920           0
    3:        2396           0
    4:        5940         778
    5:       19060         201
    6:        7780           0

``` r
nrow(sim10)
```

    [1] 28430

``` r
table(sim10$M)
```


       1    2    3    4    5    6    7    8    9   10 
    2843 2843 2843 2843 2843 2843 2843 2843 2843 2843 

# Advanced fusion

The fusionModel package allows variables to be fused individually and/or
in blocks. When variables are fused in a block, the fused values always
consist of records drawn directly from the donor. This means that fused
results for “blocked” variables will always consist of a “real” joint
observation. If all of the fusion variables are in a single block, then
the operation is equivalent to sampling complete records from the donor
(in terms of the fusion variables). This is generally only advisable
when the donor sample size is at least as large as the recipient.

There may be obvious reasons (usually having to do with the structure or
definition of certain variables) to fuse some variables as a block. But,
in general, it is not obvious when or why or which variables to block.
Even in the absence of blocks (i.e. one-by-one fusion), it is usually
unclear how to *chain* (sequence) the fusion process. Presumably, it
makes sense to fuse some variables ahead of others, but the literature
has very little to say about this question.

The experimental `blockchain()` function attempts to provide a plausible
fusion sequence (i.e. chain) and, optionally, a strategy for assigning
variables to blocks. The algorithm uses cross-validated LASSO models fit
via [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html).
It first builds a “complete” model for each *y* (fusion) variable using
all other variables as predictors; the cross-validated model skill is
the maximum possible. Next, a “minimal” model is fit using only the *x*
predictors; the model skill is divided by the maximum to create a
“score” metric. The *y* with the maximum score is assigned to the first
position in the fusion chain and included as a predictor in all
subsequent models. The remaining *y* are assigned to the chain in the
same, greedy fashion.

When a fusion variable (*y1*) is selected for inclusion in the chain,
its score is compared to that of the previous iteration; i.e. its score
prior to including the preceding fusion variable (*y0*) as a predictor.
If the score does not improve significantly, then *y1* is grouped into a
block with *y0*. The general logic here is that chaining makes sense
if/when it adds substantial explanatory power (i.e. when *y0* helps
predict *y1*). If chaining does not appear to do this, then the default
preference is to fuse the variables jointly as a block. See
`?blockchain` for options and additional details.

Let’s use `blockchain()` to get some guidance on how/whether to sequence
and block the `fusion.vars` from the previous example.

``` r
fchain <- blockchain(data = donor,
                     y = fusion.vars,
                     x = predictor.vars)
```

    Preparing data
    Constructing cross-validation folds
    Fitting complete models
    Determining order and blocks

``` r
fchain
```

    [[1]]
    [1] "insulation" "aircon"    

    [[2]]
    [1] "electricity" "natural_gas"

    [[3]]
    [1] "square_feet"

The resulting `fchain` list suggests that “insulation” and “aircon” be
fused in the initial block, followed by “electricity” and “natural_gas”
in a separate block, and then “square_feet” on its own. You may get
slightly different results, since `blockchain()` relies on random
cross-validation. We can then pass the suggested blocking and chaining
strategy to `train()`, whose `y` argument accepts a list for the
purposes of specifying variable blocks.

For this call to `train()`, we will also specify a set of
hyperparameters to search over when training each LightGBM gradient
boosting model (see `?train` for details). Since this requires
additional computation, the `cores` argument is used to enable parallel
processing.

``` r
# Train a fusion model with variable blocks
fsn.model <- train(data = donor, 
                   y = fchain, 
                   x = predictor.vars,
                   file = "test_model.fsn",
                   hyper = list(boosting = c("gbdt", "goss"),
                                num_leaves = c(10, 30),
                                feature_fraction = c(0.7, 0.9)),
                   cores = 2)
```

    5 fusion variables
    13 initial predictor variables
    2843 observations
    Using OpenMP multithreading within LightGBM (2 cores)
    Training step 1 of 3: insulation, aircon
    Training step 2 of 3: electricity, natural_gas
    Training step 3 of 3: square_feet
    Fusion model saved to: /home/kevin/Documents/Projects/fusionModel/test_model.fsn 
    Total processing time: 19.8 secs 

And then fuse multiple implicates, per usual. Note that the order of the
columns in `sim10` reflects the order of the variables passed in
`fchain`.

``` r
# Fuse multiple implicates to the recipient
sim10 <- fuse(data = recipient, 
              file = fsn.model,
              k = 5,
              M = 10)
```

    5 fusion variables
    13 initial predictor variables
    2843 observations
    Generating 10 implicates 
    Fusion step 1 of 3: insulation, aircon
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Fusion step 2 of 3: electricity, natural_gas
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Fusion step 3 of 3: square_feet
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Total processing time: 5.58 secs 

``` r
head(sim10)
```

       M           insulation                                   aircon electricity
    1: 1       Well insulated Individual window/wall or portable units       13920
    2: 1     Poorly insulated          Central air conditioning system       17400
    3: 1 Adequately insulated Individual window/wall or portable units        1480
    4: 1 Adequately insulated          Central air conditioning system        6690
    5: 1 Adequately insulated          Central air conditioning system       21300
    6: 1     Poorly insulated Individual window/wall or portable units       16600
       natural_gas square_feet
    1:         0.0        1549
    2:         0.0        1938
    3:        99.5         413
    4:       534.0        2081
    5:         0.0        2605
    6:         0.0        1872

Note that this time the *k*-NN step is used to select simulated values
for “insulation” and “aircon”. This is because *k*-NN is always used
when simulating a block of variables, regardless of type. In this
particular case, the conditional class probabilities of both variables
are used to determine the nearest neighbors in the donor data.

The same is true for “electricity” and “natural_gas”, though they use
conditional means and quantiles to identify neighbors in the donor. We
can confirm that the simulated values consist only of observed
combinations drawn from the donor (i.e. sampled from the observed joint
distribution).

``` r
# This will return TRUE, confirming simulated outcomes are drawn from the donor
v <- c("electricity", "natural_gas")
check <- data.frame(sim10)
check <- unique(check[v])
check <- apply(check, 1, paste, collapse = "_")
all(check %in% apply(donor[v], 1, paste, collapse = "_"))
```

    [1] TRUE

# Analyzing fused data

The `analyze()` function provides a convenient way to calculate point
estimates and associated standard errors and confidence intervals when
using multiple implicates. It can return means, proportions, and linear
regression coefficients (see `?analyze` for details). The standard
errors are “pooled” ([Rubin
1987](https://onlinelibrary.wiley.com/doi/book/10.1002/9780470316696))
and account for variation both within and across implicates. In general,
more implicates give more reliable point estimates and standard errors.

For example, to calculate the mean value of the “electricity” variable
for the recipient dataset using the multiple implicates in `sim10`, we
can do the following.

``` r
analyze(electricity ~ 1,
        implicates = sim10, 
        donor.N = nrow(donor))
```

    Assuming uniform 'sample_weights`

    # A tibble: 1 × 10
      response    metric estimate std_e…¹ lower…² upper…³ stati…⁴ pvalue  degf  nobs
      <chr>       <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl> <int> <int>
    1 electricity mean     11022.    141.  10745.  11300.    78.3      0   231  2843
    # … with abbreviated variable names ¹​std_error, ²​lower_ci, ³​upper_ci,
    #   ⁴​statistic

When the response variable is categorical, `analyze()` returns the
proportions associated with each factor level.

``` r
analyze(aircon ~ 1,
        implicates = sim10, 
        donor.N = nrow(donor))
```

    Assuming uniform 'sample_weights`

    # A tibble: 4 × 10
      response level      estim…¹ std_e…² lower…³ upper…⁴ stati…⁵ pvalue  degf  nobs
      <chr>    <fct>        <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl> <int> <int>
    1 aircon   Central a…  0.640  0.0105   0.619   0.661    61.0       0   117  1820
    2 aircon   Individua…  0.193  0.00877  0.175   0.210    22.0       0   117   548
    3 aircon   Both a ce…  0.0377 0.00457  0.0286  0.0467    8.25      0   117   107
    4 aircon   No air co…  0.130  0.00743  0.115   0.144    17.4       0   117   368
    # … with abbreviated variable names ¹​estimate, ²​std_error, ³​lower_ci,
    #   ⁴​upper_ci, ⁵​statistic

If we want to perform an analysis across subsets of the recipient
population – for example, calculate the distribution of “aircon”
outcomes by household “income” – we can use the `by` and `static`
arguments.

``` r
analyze(aircon ~ 1,
        implicates = sim10, 
        donor.N = nrow(donor),
        static = recipient,
        by = "income")
```

    Assuming uniform 'sample_weights`

    # A tibble: 32 × 11
       income    respo…¹ level estim…² std_e…³ lower…⁴ upper…⁵ stati…⁶  pvalue  degf
       <ord>     <chr>   <fct>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl> <int>
     1 Less tha… aircon  Cent…  0.444  0.0251  0.394    0.493    17.7  0         112
     2 Less tha… aircon  Indi…  0.330  0.0237  0.282    0.377    13.9  0         112
     3 Less tha… aircon  Both…  0.0285 0.00947 0.00968  0.0472    3.00 0.00329   112
     4 Less tha… aircon  No a…  0.199  0.0218  0.155    0.242     9.09 0         112
     5 $20,000 … aircon  Cent…  0.588  0.0240  0.541    0.636    24.6  0         112
     6 $20,000 … aircon  Indi…  0.227  0.0187  0.190    0.264    12.1  0         112
     7 $20,000 … aircon  Both…  0.0406 0.0112  0.0185   0.0628    3.64 0.00042   112
     8 $20,000 … aircon  No a…  0.144  0.0173  0.110    0.178     8.35 0         112
     9 $40,000 … aircon  Cent…  0.620  0.0278  0.565    0.676    22.3  0         112
    10 $40,000 … aircon  Indi…  0.208  0.0237  0.161    0.255     8.79 0         112
    # … with 22 more rows, 1 more variable: nobs <int>, and abbreviated variable
    #   names ¹​response, ²​estimate, ³​std_error, ⁴​lower_ci, ⁵​upper_ci, ⁶​statistic
    # ℹ Use `print(n = ...)` to see more rows, and `colnames()` to see all variable names

Finally, any linear regression model can be specified using the formula
interface – just as in
[lm()](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/lm).

``` r
analyze(electricity ~ square_feet + hh_size,
        implicates = sim10, 
        donor.N = nrow(donor),
        static = recipient)
```

    Assuming uniform 'sample_weights`

    # A tibble: 3 × 11
      response    metric term   estim…¹ std_e…² lower…³ upper…⁴ stati…⁵ pvalue  degf
      <chr>       <chr>  <chr>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl> <int>
    1 electricity coef   (Inte… 4216.   321.    3580.   4852.      13.1      0   120
    2 electricity coef   hh_si… 1094.   106.     884.   1304.      10.3      0   120
    3 electricity coef   squar…    1.90   0.109    1.69    2.12    17.5      0   120
    # … with 1 more variable: nobs <int>, and abbreviated variable names ¹​estimate,
    #   ²​std_error, ³​lower_ci, ⁴​upper_ci, ⁵​statistic
    # ℹ Use `colnames()` to see all variable names

Happy fusing!
