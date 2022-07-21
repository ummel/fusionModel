fusionModel
================
Kevin Ummel (<ummel@berkeley.edu>)

-   [Overview](#overview)
-   [Motivation](#motivation)
-   [Methodology](#methodology)
-   [Installation](#installation)
-   [Simple fusion](#simple-fusion)
-   [Advanced fusion](#advanced-fusion)
-   [Analyzing fused data](#analyzing-fused-data)

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
    (![\\hat{Y}\_{d}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Bd%7D "\hat{Y}_{d}"))
    and recipient
    (![\\hat{Y}\_{r}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Br%7D "\hat{Y}_{r}")).
3.  A “real” donor value from
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    is selected for each recipient record, based on the similarity of
    ![\\hat{Y}\_{d}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Bd%7D "\hat{Y}_{d}")
    and
    ![\\hat{Y}\_{r}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BY%7D_%7Br%7D "\hat{Y}_{r}").

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
    ![Y\\neq0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%5Cneq0 "Y\neq0").

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
provide a data-driven answer to these question. See the [Advanced
fusion](#advanced-fusion) section for more information.

# Installation

``` r
devtools::install_github("ummel/fusionModel")
library(fusionModel)
```

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
d <- seq(1, nrow(recs), by = 2)

# Create donor and recipient datasets
donor <- select(recs[d, ], 2:13, square_feet, electricity, natural_gas, insulation, aircon)
recipient <- select(recs[-d, ], 2:13)

# Specify fusion and shared/common predictor variables
predictor.vars <- names(recipient)
fusion.vars <- setdiff(names(donor), predictor.vars)
```

The `recipient` dataset contains 12 variables that are shared with
`donor`. These shared “predictor” variables provide a statistical link
between the two datasets. fusionModel exploits the information in these
shared variables.

``` r
predictor.vars
```

     [1] "income"      "age"         "race"        "education"   "employment" 
     [6] "hh_size"     "division"    "urban_rural" "climate"     "renter"     
    [11] "home_type"   "year_built" 

There are 5 “fusion variables” unique to `donor`. These are the
variables that will be fused to `recipient`. This includes a mix of
continuous and categorical (factor) variables.

``` r
# The variables to be fused
sapply(donor[fusion.vars], class)
```

    $square_feet
    [1] "integer"

    $electricity
    [1] "numeric"

    $natural_gas
    [1] "numeric"

    $insulation
    [1] "ordered" "factor" 

    $aircon
    [1] "factor"

We build our fusion model using the `train()` function. The minimal
usage is shown below. See `?train` for additional function arguments and
options. This results in a “.fsn” (fusion) object being saved to
“test_model.fsn” in the current working directory.

``` r
# Build the fusion model (see ?train)
fsn.model <- train(data = donor, 
                   y = fusion.vars, 
                   x = predictor.vars,
                   file = "test_model.fsn")
```

    5 fusion variables
    12 initial predictor variables
    2843 observations
    Building LightGBM models...
    Fusion model saved to: /home/kevin/Documents/Projects/fusionModel/test_model.fsn 

To fuse variables to `recipient`, we simply pass the recipient data and
path of the .fsn model to the `fuse()` function. Each variable specified
in `fusion.vars` is fused in the order provided. The default setting
`k = 5` means that the five nearest neighbors in the donor are
identified for each recipient observation, and one of them is randomly
selected to provide the fused value.

``` r
# Fuse 'fusion.vars' to the recipient (see ?fuse)
sim <- fuse(data = recipient, 
            file = fsn.model,
            k = 5)
```

    Fusion 1 of 5: square_feet
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Fusion 2 of 5: electricity
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Fusion 3 of 5: natural_gas
    -- Predicting LightGBM models
    -- Finding nearest neighbors
    -- Simulating fused values
    Fusion 4 of 5: insulation
    -- Predicting LightGBM models
    -- Simulating fused values
    Fusion 5 of 5: aircon
    -- Predicting LightGBM models
    -- Simulating fused values

Notice that the status messages for the final two fusion variables –
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

       square_feet electricity natural_gas           insulation
    1:        2904       21940           0       Well insulated
    2:        2721       14270           0       Well insulated
    3:        4100        3560         914     Poorly insulated
    4:        4100       14060         976       Well insulated
    5:         853        2040         427 Adequately insulated
    6:        3288       13700        1016       Well insulated
                                           aircon
    1: Both a central system and individual units
    2:            Central air conditioning system
    3:   Individual window/wall or portable units
    4: Both a central system and individual units
    5:   Individual window/wall or portable units
    6:            Central air conditioning system

We can do some quick sanity checks to compare the distribution of the
fusion variables in `donor` with those in `sim`. This, at least,
confirms that the fusion output is not obviously wrong.

``` r
sim <- data.frame(sim)

# Compare means of the continuous variables
cbind(donor = colMeans(donor[fusion.vars[1:3]]), sim = colMeans(sim[fusion.vars[1:3]]))
```

                     donor        sim
    square_feet  2088.2807  2077.5202
    electricity 11104.1047 11153.1611
    natural_gas   335.0414   331.7955

``` r
# Compare frequencies of categorical variable classes
cbind(donor = table(donor$insulation), sim = table(sim$insulation))
```

                         donor  sim
    Not insulated           41   32
    Poorly insulated       474  471
    Adequately insulated  1366 1347
    Well insulated         962  993

``` r
cbind(donor = table(donor$aircon), sim = table(sim$aircon))
```

                                               donor  sim
    Central air conditioning system             1772 1826
    Individual window/wall or portable units     576  502
    Both a central system and individual units   128  130
    No air conditioning                          367  385

And we can look at kernel density plots of the non-zero values for the
continuous variables to see if the univariate distributions in `donor`
are generally preserved in `sim`.

![](man/figures/README-unnamed-chunk-9-1.png)<!-- -->

Since each call to `fuse()` returns a single probabilistic version of
the fusion variables, we generally want to create multiple versions –
called *implicates* – in order to reduce bias in point estimates and
calculate associated uncertainty. The `fuseM()` function is used to
generate multiple implicates. Here we generate `M = 10` implicates.

``` r
# Fuse multiple implicates to the recipient (see ?fuseM)
sim10 <- fuseM(data = recipient, 
               file = fsn.model,
               k = 5,
               M = 10)
```

Note that each implicate in `sim10` is identified by the “M”
variable/column.

``` r
nrow(sim10)
```

    [1] 28430

``` r
table(sim10$M)
```


       1    2    3    4    5    6    7    8    9   10 
    2843 2843 2843 2843 2843 2843 2843 2843 2843 2843 

``` r
head(sim10)
```

       M square_feet electricity natural_gas           insulation
    1: 1        1240       23270           0       Well insulated
    2: 1        2560       10230         358       Well insulated
    3: 1        2078        5310        1054 Adequately insulated
    4: 1        1124        7940         911 Adequately insulated
    5: 1         841        4120         279       Well insulated
    6: 1        4304       18100        1247       Well insulated
                                         aircon
    1:          Central air conditioning system
    2: Individual window/wall or portable units
    3:          Central air conditioning system
    4:          Central air conditioning system
    5:          Central air conditioning system
    6:          Central air conditioning system

# Advanced fusion

The fusionModel package allow for variables to be fused individually
and/or in blocks. When variables are fused in a block, the fused values
always consist of records drawn directly from the donor. This means that
fused results for “blocked” variables will always consist of a “real”
joint observation. If *all* of the fusion variables are in a single
block, then the operation is equivalent to sampling “complete” records
from the donor (in terms of the fusion variables). This is generally
only advisable when the donor sample size is at least as large as the
recipient.

There may be obvious reasons (usually having to do with the structure or
definition of certain variables) to fuse some variables as a block. But,
in general, it is not obvious when or why or which variables to block.
Even in the absence of blocks (i.e. one-by-one fusion), it is usually
unclear how to *sequence* the fusion process. Presumably, it makes sense
to fuse some variables ahead of others, but the literature has very
little to say about this question.

The experimental `blockchain()` function attempts to provide a plausible
fusion sequence (i.e. chain) and, optionally, a strategy for assigning
variables to blocks. The algorithm uses cross-validated LASSO models fit
via [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html).
It first builds a “complete” model for each *y* (fusion) variable using
*all* other variables as predictors; the cross-validated model skill is
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
                     x = predictor.vars,
                     delta = 0.01,
                     criterion = "min")
```

``` r
fchain
```

    [[1]]
    [1] "aircon"     "insulation"

    [[2]]
    [1] "electricity"

    [[3]]
    [1] "square_feet"

    [[4]]
    [1] "natural_gas"

The resulting `fchain` list suggests that “aircon” and “insulation”
should be fused jointly in the initial block followed by the remaining
fusion variables in a preferred order. You may get slightly different
results, since `blockchain()` relies on random cross-validation. We can
then pass the suggested blocking and chaining strategy to `train()`,
whose `y` argument accepts a list for the purposes of specifying
variable blocks.

For this call to `train()`, we will also specify a set of
hyperparameters to search over when training each LightGBM gradient
boosting model (see `?train` for details). Since this requires
additional computation, the `cores` argument is used to enable parallel
processing.

``` r
# Build a fusion model with variable blocks
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
    12 initial predictor variables
    2843 observations
    Building LightGBM models...
    Fusion model saved to: /home/kevin/Documents/Projects/fusionModel/test_model.fsn 

And then fuse multiple implicates, per usual. Note that the order of the
columns in `sim10` reflects the order of the variables passed in
`fchain`.

``` r
# Fuse multiple implicates to the recipient (see ?fuseM)
sim10 <- fuseM(data = recipient, 
               file = fsn.model,
               k = 5,
               M = 10)
```

``` r
head(sim10)
```

       M                          aircon           insulation electricity
    1: 1 Central air conditioning system       Well insulated       31730
    2: 1 Central air conditioning system       Well insulated        9290
    3: 1 Central air conditioning system Adequately insulated        9350
    4: 1 Central air conditioning system Adequately insulated        8410
    5: 1 Central air conditioning system     Poorly insulated        6180
    6: 1 Central air conditioning system       Well insulated        7690
       square_feet natural_gas
    1:        2118         644
    2:        1550         856
    3:        1008         572
    4:        3468        1020
    5:         880         378
    6:        4761        1210

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

    # A tibble: 1 × 10
      response    metric estimate std_error lower_ci upper_ci statistic pvalue  degf
      <chr>       <chr>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>  <dbl> <int>
    1 electricity mean     11096.      147.   10806.   11387.      75.6      0   137
    # … with 1 more variable: nobs <int>

When the response variable is categorical, `analyze()` automatically
returns the proportions associated with each factor level.

``` r
analyze(aircon ~ 1,
        implicates = sim10, 
        donor.N = nrow(donor))
```

    # A tibble: 4 × 10
      response level     estimate std_error lower_ci upper_ci statistic pvalue  degf
      <chr>    <fct>        <dbl>     <dbl>    <dbl>    <dbl>     <dbl>  <dbl> <int>
    1 aircon   Central …   0.652    0.0105    0.631    0.673      62.1       0   109
    2 aircon   Individu…   0.196    0.00917   0.178    0.215      21.4       0    64
    3 aircon   Both a c…   0.0343   0.00454   0.0249   0.0436      7.54      0    25
    4 aircon   No air c…   0.118    0.00724   0.103    0.132      16.2       0    68
    # … with 1 more variable: nobs <int>

If we want to perform an analysis across subsets of the recipient
population – for example, calculate the mean of “electricity” by
household “income” – we can use the `by` and `static` arguments.

``` r
analyze(aircon ~ 1,
        implicates = sim10, 
        donor.N = nrow(donor),
        static = recipient,
        by = "income")
```

    # A tibble: 32 × 11
       income   response level estimate std_error lower_ci upper_ci statistic pvalue
       <ord>    <chr>    <fct>    <dbl>     <dbl>    <dbl>    <dbl>     <dbl>  <dbl>
     1 Less th… aircon   Cent…   0.636    0.0247   0.587     0.685      25.8  0     
     2 Less th… aircon   Indi…   0.200    0.0251   0.148     0.252       7.99 0     
     3 Less th… aircon   Both…   0.0329   0.0116   0.00443   0.0614      2.84 0.0302
     4 Less th… aircon   No a…   0.131    0.0206   0.0877    0.173       6.35 0     
     5 $20,000… aircon   Cent…   0.641    0.0213   0.599     0.683      30.1  0     
     6 $20,000… aircon   Indi…   0.207    0.0195   0.168     0.247      10.6  0     
     7 $20,000… aircon   Both…   0.0326   0.00976  0.00997   0.0551      3.34 0.0106
     8 $20,000… aircon   No a…   0.119    0.0163   0.0856    0.153       7.30 0     
     9 $40,000… aircon   Cent…   0.65     0.0255   0.599     0.701      25.5  0     
    10 $40,000… aircon   Indi…   0.202    0.0228   0.155     0.248       8.84 0     
    # … with 22 more rows, and 2 more variables: degf <int>, nobs <int>

Finally, any linear regression model can be specified using the formula
interface – just as in
[lm()](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/lm).

``` r
analyze(electricity ~ square_feet + hh_size,
        implicates = sim10, 
        donor.N = nrow(donor),
        static = recipient)
```

    # A tibble: 3 × 11
      response    metric term  estimate std_error lower_ci upper_ci statistic pvalue
      <chr>       <chr>  <chr>    <dbl>     <dbl>    <dbl>    <dbl>     <dbl>  <dbl>
    1 electricity coef   (Int…  6519.     359.     5808.    7229.      18.2    0    
    2 electricity coef   hh_s…   -35.7    106.     -248.     177.      -0.336  0.738
    3 electricity coef   squa…     2.23     0.119     1.99     2.47    18.7    0    
    # … with 2 more variables: degf <int>, nobs <int>

Happy fusing!
