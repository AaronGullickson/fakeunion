# fakeunion - An R package for generating counterfactual unions

[![DOI](https://zenodo.org/badge/88232957.svg)](https://zenodo.org/badge/latestdoi/88232957)

This R package contains functions that will generate a set of "fake" unions (e.g. marriages, cohabitations), given a set of real unions and alternate partners. The key function `generateCouples` will sample from the list of alternate partners to create a dataset that contains both real and fake unions. This dataset can then be used to run a fixed effect conditional logit model to predict what partner characteristics are likely to be associated with a real match.

## Installation

The package can be installed in R from GitHub with the `devtools` library:

```
library(devtools)
install_github("AaronGullickson/fakeunion")
```

## Usage

The examples section of the help file for `generateCouples` contains a simple example showing how to use the program, but I provide some further details here.

There are three core data frame objects that are required:

- **actual**: this data frame contains information about the real unions and the characteristics of spouses in those unions. This dataset should contain a variable that identifies the geographic unit that unions should be sampled within (e.g. state, city, county). It should also contain unique id variables that identify the two partners and end with "h" and "w" (e.g. "idh" and "idw"). Further data can be included here on the couple or the partner characteristics. All partner characteristic variables should end with either an "h" or "w" to identify which partner they belong to.
- **men**: This data frame contains information about potential alternate men that could have been selected as partners. This dataset should contain the same geographic unit variable and an unique id that is labeled the same way as `actual` (e.g. "idh"). It should also include characteristics that match the "h" partner characteristics in `actual`. This dataset can also contain a variable for weights that can be used to sample individuals with different probability. The weight variable should not be identified with an "h" ending.
- **women**: This data frame contains information about potential alternate women that could have been selected as partners. It should have the same format as `men` but variables should end in "w".

Here are examples of what the format should look like based on the example acs data supplied in the package:

**Actual**:

```
   state              idw   agew racew             idh ageh raceh
  New York  2011.1.856963.2   53     W 2011.1.856963.1   69     W
  New York  2013.1.815005.1   30     W 2013.1.815005.2   33     W
  New York  2014.1.877381.2   38     W 2014.1.877381.1   39     W
California  2015.1.100895.1   52     W 2015.1.100895.2   46     W
California  2011.1.199641.2   34     W 2011.1.199641.1   31     W
  New York  2011.1.870767.1   48     W 2011.1.870767.2   39     W
```

**Men**:

```
    state perwt             idh ageh raceh
California    58 2012.1.193121.1   41     B
California    73 2015.1.160865.3   23     W
  New York    72 2013.1.872930.1   70     W
  New York   112 2011.1.859134.3   30     W
California    57 2011.1.199178.1   47     W
  New York    56 2012.1.873801.3   33     W
```

**Women:**

```
     state perwt             idw agew racew
California    40 2015.1.192456.3   18     W
  New York   466 2015.1.850982.2   26     W
  New York    79 2013.1.879650.1   61     W
  New York   103 2014.1.816232.1   63     B
California   110 2011.1.140627.1   71     W
  New York    79 2014.1.809120.1   85     W
```

The main function `generateCouples` can then be called by specifying these datasets as well as the geographic cluster for sampling and the number of partners to be sampled. Here is an example using the ACS data provided in the package

```
markets <- generateCouples(3,acs.couples,acs.men,acs.women,
                           geo="state",weight="perwt")
```

This program will randomly choose either the husband or wife of the actual couples and sample three random partners from within the state without replacement. The sampling will be done using the weights provided by the `perwt` variable.

### Processing time

I have tested out the processing time of the `generateCouples` under a variety of scenarios. The American Community Survey data I use has 73,536 actual couples, 2.35 million alternate male partners, and 2.76 million alternate female partners. I test clustering at three different levels: metro area, state, and country. Since all respondents were resident in the USA, the "country" level is effectively a single cluster.

I also created a smaller dataset by sampling 10% of the data from my actual couples and alternates. Finally, I tested the effect of taking probability samples by weighting and using flat sampling. The results are presented in the table below. I record processing time in seconds. The computations were run from a high-end linux computer.

|size   |	cluster	| weights	   | N3	    | N5	  | N10	  | N30	  | N50    |
|:------|:--------|:-----------|-------:|------:|------:|------:|-------:|
|Small	| metro	  | unweighted | 17.36	| 16.76 |	26.10	| 64.17	| 102.68 |
|Small	| metro	  | weighted	 | 18.40	| 16.89	| 26.39	| 64.52	| 101.88 |
|Small	| state	  | unweighted |	4.36	|  4.28	|  6.68	| 17.12	|  25.81 |
|Small	| state	  | weighted	 |  4.36	|  4.27	|  6.67	| 17.18	|  25.48 |
|Small	| country	| unweighted |	6.57	|  6.64	|  6.83	|  7.51	|   8.22 |
|Small	| country	| weighted	 |  6.57	|  6.64	|  6.83	|  7.51	|   8.22 |
|Large	| metro	  | unweighted |173.44	|175.10	|270.74	|666.56	|1091.80 |
|Large	| metro	  | weighted	 |167.32	|171.18	|273.09	|659.00 |1088.92 |
|Large	| state	  | unweighted | 58.70	| 61.36	| 87.54	|183.15	| 287.52 |
|Large	| state	  | weighted	 | 58.58	| 61.36	| 86.14	|185.42	| 285.77 |
|Large	| country	| unweighted |625.25	|624.38	|626.86	|638.11	| 648.72 |
|Large	| country	| weighted	 |625.25	|624.38	|626.86	|638.11	| 648.72 |


### Conditional Logit Models

The data that is produced by `generateCouples` can be used directly in conditional logit models where the `group` variable is used as the fixed effect. These models can be estimated in R using the `clogit` command in the `survival` library and specifying the `group` by the `strata` function, like so:

```
library(survival)
clogit(choice~I(ageh-agew)+I((ageh-agew)^2)+I(raceh!=racew)
              +strata(group), data=markets)
```

This model can also be estimated in Stata using `xtlogit`.

### Convenience functions

This package also contains two convenience functions for calculating the conditional logit models. First, the function `poolChoiceModel` can be used to run the same conditional logit model on a list of sampled datasets. Because the dataset generated by `generateCouples` is a sample that will vary each time it is conducted, it may be useful to use a technique similar to multiple imputation and conduct *m* parallel analyses that can then be pooled together to account for additional between sample variation. This function will do that for the user based on a formula and a list of generated datasets. An example is given in the help file.

The function `createExogamyTerms` will create a set of interaction terms between partners characteristics but will also create an "endogamous" category when partners share the same characteristic and set this category as the reference. The user can also specify that the exogamy terms should be symmetric such that it does not matter which spouse had which characteristic.
