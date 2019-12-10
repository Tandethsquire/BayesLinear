// Imports from other packages
const base = require('./boilerplate');
const mat = require('./matrixAlgebra');

/**
 * The class definition for a random variable.
 * In the absence of proper function overloading in JS, the definition takes
 * two object: a name and a dictionary.
 * The dictionary is either the expectation and variance, or it's three arrays:
 * the dependent variables, the coefficients, and any non-dependent variables.
 * For instance, a variable X1 with expectation 1 and variance 2 would be created
 * as rV('X1', {'mu': 1, 'sigma': 2}); if X1 has been thus defined, then we could
 * create X2 = 4*X1 as rV('X2', {'rvs': [X1], 'coeffs': [4]}).
 * @param  {String} name The name of the random variable
 * @param  {Dict} args A dictionary of arguments to pass (see above for details)
 * @return {rV}      The output random variable object.
 */
var rV = function(name, args) {
  if (args.mu !== undefined && args.sigma !== undefined) {
    if (isNaN(args.sigma) || isNaN(args.mu))
      throw new Error('Either the expectation or variance are non-numeric');
    this.nm = name;
    this.exp = args.mu;
    // This creates a dictionary of covariances between the rV and others; it starts
    // by setting Var(rV)=sigma.
    this.cov = {};
    this.cov[name] = args.sigma;
  }
  else if (args.rvs !== undefined && args.coeffs !== undefined) {
    var rvs = args.rvs, coeffs = args.coeffs;
    var exp = base.mRound(rvs.reduce((a,b,i) => a + coeffs[i]*b.exp, 0),7);
    var sigma = base.mRound(rvs.reduce((a,b,i) => a + coeffs[i]*calculate_covariance(b,rvs,coeffs), 0),7);
    this.nm = name;
    this.exp = exp;
    this.cov = {};
    this.cov[name] = base.mRound(sigma,7);
    for (var i=0; i<rvs.length; i++) {
      var s = rvs[i];
      this.setCov(s,base.mRound(calculate_covariance(s,rvs,coeffs),7));
    }
    if (args.independent !== undefined)
    {
      var ind = args.independent;
      for (var i=0; i<ind.length; i++)
        this.setCov(ind[i],base.mRound(calculate_covariance(ind[i],rvs,coeffs),7));
    }
  }
  else {
    throw new Error('Random variable declaration incomplete or incorrect.');
  }
}

/**
 * Sets the covariance between two random variables. If oneway=true, the covariance
 * is only added to 'that' dictionary. This is useful for ephemeral variables where
 * we don't need it beyond an intermediate calculation, and so don't want it in
 * 'proper' variable dictionaries.
 * @param  {rV}  that           The other random variable
 * @param  {Number}  cov            The covariance
 * @param  {Boolean} [oneway=false] Should we only set the covariance in that?
 */
rV.prototype.setCov = function(that, cov, oneway = false) {
  if (!(that instanceof rV))
    throw new Error('Cannot set a covariance between a random variable and a non-random one');
  this.cov[that.nm] = cov;
  if (!oneway)
    that.cov[this.nm] = cov;
}

/**
 * Standardises a random variable. If an rV X has expectation M and variance S,
 * then S(X)=(X-M)/sqrt(S).
 * @param  {Array[rV]} others Any random variables that we want to calculate covariance wrt.
 * @return {rV}        The standardised rV.
 */
rV.prototype.standardise = function(others) {
  var constRV = new rV('temp_variable', {'mu': -1*this.exp, 'sigma': 0});
  var totalRV = [this].concat(others);
  totalRV.forEach( function (s) {
    constRV.setCov(s,0);
  });
  var sigma = this.cov[this.nm];
  var result = new rV(`S(${this.nm})`, {'rvs': [this, constRV], 'coeffs': [1/Math.sqrt(sigma), 1/Math.sqrt(sigma)], 'independent': others});
  ([result].concat(totalRV)).forEach( function(s) {
    delete s.cov['temp_variable'];
  })
  return result;
}

/**
 * Calculates covariance between a random variable and a linear combination thereof.
 * Given rVs X1, X2, X3 and an rV Y = X1+2*X2-4*X3, then Cov(X1,Y) is calculated as
 * calculate_covariance(X1,[X1,X2,X3],[1,2,-4]).
 * @param  {rV} rv     The primitive random variable
 * @param  {Array[rV]} rvs    The collection of rVs that make up the linear combination.
 * @param  {Array[Number]} coeffs The corresponding coefficients
 * @return {Number}        The covariance.
 */
var calculate_covariance = function(rv, rvs, coeffs) {
  return rvs.reduce((a,b,i) => a + coeffs[i]*rv.cov[b.nm], 0);
}

var rvs_from_matrix = function(names, expectations, variances) {
  var rvArr = [];
  if (names.length != expectations.length || names.length != variances.length)
    throw new Error("The number of named random variables does not match the specification.")
  for (var i=0; i<names.length; i++)
  {
    var tempRV = new rV(names[i], {mu: expectations[i], sigma: variances[i][i]});
    rvArr.push(tempRV);
  }
  for (var i=0; i<rvArr.length; i++)
  {
    for (var j=0; j<i; j++)
    {
      rvArr[i].setCov(rvArr[j], variances[i][j]);
    }
  }
  return rvArr;
}

/**
 * Creates a variance (or covariance) matrix from two collections of rVs.
 * @param  {Array[rV]} args1 The first set (creating the rows of the matrix)
 * @param  {Array[rV]} args2 The second set (creating the columns)
 * @return {Array[Array[Number]]}       The (co)variance matrix.
 */
var build_variance_matrix = function (args1, args2) {
  var outarr = new Array(args1.length).fill(0).map(s => new Array(args2.length).fill(0));
  return outarr.map((s,i) => s.map((elem,j) => args1[i].cov[args2[j].nm]));
}

module.exports = {rV, build_variance_matrix, rvs_from_matrix};
