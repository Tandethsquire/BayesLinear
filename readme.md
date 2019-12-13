# Bayes Linear Statistics Package #
---
### A set of functions for doing Bayes Linear calculations ###

## Introduction ##
Typically in statistics, when we want to model a physical situation, we start by stipulating probability distributions and from there, determine expectation and variance. This can be cumbersome and computationally inefficient, due to the complexity of a general probability distribution. If we want to adjust our model due to observed information, this requires a redefinition of the probability distribution in order to refine the observable quantities.

In Bayes Linear Statistics, we instead define our belief specification purely in terms of the expectation (presumed quantities) and variance (uncertainty in the statement of the quantities). A probability distribution can be created from this, by thinking of P(X) as an indicator function applied to E(X). However, in a general setting, we only wish to modify the model based on observed quantities, and so the probability distribution is not needed; the observations motivate a change in the expectation and variance directly. Similarly, we can use the expectations and variances of some set of 'data' (random variables for which we expect to obtain observations) to adjust beliefs, even in the absence of data.

Such an approach allows us to consider which models will most affect the understanding of our beliefs, by considering the amount of variance resolved. Given a set of beliefs, B, and some data quantities, D, we can obtain an adjustment of the expectation and variance of B through some straightforward matrix algebra, and extract information like total variance reduction and canonical quantities of the beliefs (those combinations of elements of B that are most affected by D). We can compare multiple different belief specifications by looking at the reductions over B. If we then obtain data, we can perform diagnostics on the adjustments, and use this to motivate our choice of model.

## Structure ##
The code is in multiple parts (which will continue to be modified and re-jigged as time goes on). Files described earlier will typically (but not always) be dependencies of the later files.
- *boilerplate.js* : This contains a set of generic, all-purpose functions for doing generally boring things (e.g. finding the sign of a number, printing out a matrix in user-friendly format)
- *matrixAlgebra.js* : This does a large part of the computational heavy lifting. There are the usual matrix operations (addition, multiplication, trace etc), functions to find eigenvalues and eigenvectors, inversion of matrices, and so on. Most of the important functionality is based about QR Decomposition, which reduces a matrix to upper-triangular form.
- *varianceComparison.js* : This deals with belief comparisons. Functions here calculate canonical quantities of some chosen belief(s), find observed adjustments, work out residuals for specification diagnostics, and package this information up in an easily-exportable way.
- *RV.js* : A small file for defining and dealing with random variables. Creates objects of type `rV`, with expectation, variance, and covariances as properties.
- *uncertaintyResolution.js* : This has a similar scope to varianceComparison.js, but it centres around node resolutions and node diagnostics. Functions operating on `rV` objects are used to generate (potentially partial) transformation matrices for B given D, calculate uncertainty resolution from one random variable to another, find bearings and sizes of adjustments, find the heart of a transform, and create a tree of dependencies from the result.
- *<filename>.html*: Typically used for displaying the results of the above packages. D3.js is used heavily to create visualisations of these adjustments, dependencies, and resolutions.

## Usage ##
Any pure JavaScript files are written with command-line usage in mind, and particularly using node.js. Outputs can be written to file (in json or csv format) to be pushed to the .html visualisation.

## To do ##
- There are some aspects of `uncertaintyResolution.js` that would benefit from more testing (for example, the `path_correlation` function and related), and where it could be made more user-friendly. It would be preferable to have the user input the node structure (eg combined nodes) as a command-line argument, rather than changing the code directly.
- Data visualisation: any alternative ways to view the data is always useful. In particular, the plots from `beliefComparison.html` suffer from redundancy in displayed information
- More functionality! Suggestions welcome.

## Main File Summaries ##
### varianceComparison.js ###
Take a system of two variables, X1 and X2. Define the variance in two different specifications H1 and H2 as
- **H1**: Var(X1)=2, Var(X2)=4, Cov(X1,X2)=1, E(X1)=0, E(X2)=0;
- **H2**: Var(X1)=2, Var(X2)=2, Cov(X1,X2)=1, E(X1)=0, E(X2)=1.

Assume that observation of X1 and X2 gives sample expectations x1=3 and x2=1.
``` javascript
var U1 = [[2,1],[1,4]];
var U2 = [[2,1],[1,2]];
var E1 = [0,0];
var E2 = [0,0];
var X = [3,1];
var canonical = arrangeCanonical(canonicalQuantities(U1,U2),U1,U2,E1,E2)
console.log(canonical);
/**
 * Outputs a dictionary of quantities Z1 = -0.267261 X1 + 0.534522 X2 and Z2 = 0.707107 X1
 * with variances 1,0.4285714 and 1,1 under H1 and H2, respectively.
 */
console.log(standardised(canonical, U1, U2, E1, E2, X));
/**
 * Outputs the observed result as a dictionary, with expectations and variances
 * under both H1 and H2 for comparison:
 * observed: -0.1428571
 * exp1: 0
 * exp2: 0.2857143
 * var1 = 0.2857143
 * var2 = 0.122449
 */
```

### uncertaintyResolution.js ###
Suppose we have a variable Y which is given by a linear relation Y=a+2b+e. The variances of a,b,e are assumed to be

Var(a)=4, Var(b)=3, Var(e)=0.5, Cov(a,b)=-1, Cov(a,e)=0=Cov(b,e)

with expectations

E(a)=1, E(b)=2, E(e)=0.

e is an error or 'noise' term. Then Y is fully determined by its parents.
``` javascript
var a = new RV.rV('a',1,4), b = new rV('b',2,3);
a.setCov(b,-1);
var e = new RV.rV('e',0,0.5);
e.setCov(a,0);
e.setCov(b,0);
var Y = new RV.rV('Y', {'rvs': [a,b,e], 'coeffs': [1,2,1]});
console.log(Y);
/**
 * Output is a rV with fields
 * nm: 'Y', exp: 5, cov: {Y: 12.5, a: 2, b: 5, e: 0.5}
 */
// If parent data is not known, instead use derive_parents([a,b,e,Y])
var parentData = [
  ['a', a, []],
  ['b', b, [1]],
  ['e', e, []],
  ['Y', Y, [1,2,3]]
];
console.log(buildNodeList(parentData, [a,b,e,Y]));
/**
 * Output is a dictionary of nodes and links:
 * { nodes:
 *  [ { id: 1, name: 'a', parents: [] },
 *   { id: 2, name: 'b', parents: [Array] },
 *   { id: 3, name: 'e', parents: [] },
 *   { id: 4, name: 'Y', parents: [Array] } ],
 *links:
 * [ { source: 1, target: 2, leaving: 0.0833, arriving: 0.0833 },
 *   { source: 1, target: 4, leaving: 0.08, arriving: 0.2933 },
 *   { source: 2, target: 4, leaving: 0.6667, arriving: 0.88 },
 *   { source: 3, target: 4, leaving: 0.04, arriving: 0.04 } ] }
 */
```

It is also possible to import variance data, instead of building random variables by hand, using the `read_from_csv` function (which, at present, is a bit unstable) and the `RV.rvs_from_matrix` function. Some working code using `sample.csv` is included in the file.
