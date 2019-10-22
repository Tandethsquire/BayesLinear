/**
 * 
 * A standalone file for expectation and variance adjustement of beliefs.
 * Uses functions from boilerplate.js, matrixAlgebra.js, and RV.js.
 *
 */


// From boilerplate

/**
 * Rounds a number to a given precision, up to some constraints on closeness
 * to -1, 0 or 1.
 * @param  {Number} n    The number to round
 * @param  {Number} prec The required precision (as decimal places)
 * @return {Number}      The rounded number
 */
function mRound(n,prec)
{
  if (Math.abs(n)<0.00001) return 0;
  if (Math.abs(n-1)<0.00001) return 1;
  if (Math.abs(n+1)<0.00001) return -1;
  return parseFloat(parseFloat(n).toFixed(prec));
}
/**
 * Finds the sign of a number (-1 if negative, 1 if positive, 0 else)
 * @param  {Number} x The number
 * @return {Number}   The sign: either -1, 0, or 1
 */
var sign = function(x) {
  return ((x==0)? 0 :((x < 0) ? -1: 1));
}

// From matrixAlgebra
// Global precision for calculation
const epsilon = 0.00000001;
/**
 * Creates an identity matrix of a given size.
 * @param  {Number} size An integer indicating the size of the identity
 * @return {Array[Array[Number]]}      The size*size identity matrix.
 */
var identity = function (size) {
  var idArr = new Array(size).fill(0).map(s => new Array(size).fill(0));
  return idArr.map((s,i) => s.map((el,j) => (i==j)? 1: 0));
}

/**
 * Creates a diagonal matrix whose diagonal entries are given by a list.
 * @param  {Array[Number]} elems The diagonal elements
 * @return {Array[Array[Number]]}       The diagonal matrix
 */
var diagonal = function (elems) {
  var result = new Array(elems.length).fill(0).map(s => new Array(elems.length).fill(0))
  return result.map((s,i) => s.map((elem,j) => (i==j) ? elems[i]: 0));
}
/**
 * Tests if a matrix is upper triangular; i.e. if every entry below the leading diagonal is zero.
 * @param  {Array[Array[Number]]} mat The matrix to check
 * @return {Boolean}
 */
var upperTri = function (mat) {
  for (var col=0; col<mat[0].length; col++)
  {
    for (var row=col+1; row<mat.length; row++)
    {
      if (Math.abs(mat[row][col])>epsilon)
        return false;
    }
  }
  return true;
}
/**
 * Calculates the norm of a vector w.r.t. the standard inner product.
 * @param  {Array[Number]} vec The vector
 * @return {Number}     The vector's length
 */
var vNorm = function (vec) {
  return Math.sqrt(vec.reduce((a,b) => a+b*b, 0));
}
/**
 * Multiplies two matrices.
 * The number of columns of A must be the same as the number of rows of B.
 * @param  {Array[Array[Number]]} A The first matrix
 * @param  {Array[Array[Number]]} B The second matrix
 * @return {Array[Array[Number]]}   AB
 */
var mult = function (A, B) {
  if (A[0].length != B.length)
    return TypeError('Matrices are not compatible for multiplication.');
  var result = new Array(A.length).fill(0).map(row => new Array(B[0].length).fill(0));
  return result.map((row,i) => row.map((val,j) => A[i].reduce((sum,elem,k) => sum + elem*B[k][j], 0)));
}
/**
 * Transposes a matrix: exchanges rows for columns.
 * @param  {Array[Array[Number]]} A The matrix
 * @return {Array[Array[Number]]}   A^T
 */
var transpose = function (A)
{
  var result = new Array(A[0].length).fill(0).map(s => new Array(A.length).fill(0));
  return result.map((s,i) => s.map((el,j) => A[j][i]));
}
/**
 * Creates a Householder matrix from a given vector.
 * @param  {Array[Number]} vec The starting vector
 * @return {Array[Array[Number]]}     The output transformation matrix
 */
var make_householder = function (vec) {
  var tempv = vec.map(x => [x/(vec[0]+sign(vec[0])*vNorm(vec))]), H = identity(vec.length);
  tempv[0] = [1];
  var correction = mult(tempv, transpose(tempv)).map(x => x.map(y => 2*y/Math.pow(vNorm(tempv.map(s => s[0])),2)));
  return H.map((s,i) => s.map((el,j) => el - correction[i][j]));
}
/**
 * Performs a single Householder transformation on a matrix.
 * This transforms a column of the matrix A to be upper triangular.
 * If A is symmetric, then it also does the same to the equivalent row.
 * A is modified, and we store the transformation for later use.
 * @param  {Array[Array[Number]]} A The matrix tp apply householder to
 * @return {[Array[Array[Number]],Array[Array[Number]]]}   The transformation matrix, and the transformed A.
 */
var householder = function (A) {
  var m = A.length, n = A[0].length, Q = identity(m);
  for (var i=0; i<n-(m==n); i++)
  {
    var H = identity(m), tempvec = [];
    for (var j=i; j<m; j++)
      tempvec.push(A[j][i]);
    var hh = make_householder(tempvec);
    for (var row=i; row<m; row++)
    {
      for (var col=i; col<m; col++)
        H[col][row] = hh[row-i][col-i];
    }
    Q = mult(Q,H);
    A = mult(H,A,H);
  }
  return [Q,A];
}
/**
 * Performs QR decomposition on a matrix.
 * At each stage, we convert A to a matrix QR, where Q is upper triangular and R is orthogonal.
 * Then we calculate RQ and repeat until A is upper triangular.
 * If A is symmetric, then this process generates a diagonal matrix and an orthogonal matrix;
 * these are the eigenvalues and eigenvectors of A.
 * There is a fallback of 200 iterations, in case we cannot quite achieve the needed precision.
 * @param  {Array[Array[Number]]} A The starting matrix
 * @return {[Array[Array[Number]],Array[Array[Number]]]}   The transformation matrix, and the result
 */
var QRDecompose = function (A) {
  var hfirst = householder(A);
  var Q = hfirst[0], R = hfirst[1], wQ = identity(Q.length), mat = A, fallback = 200, index = 0;
  while (!upperTri(A) && index<fallback)
  {
    var hnext = householder(mat);
    Q = hnext[0], R = hnext[1];
    wQ = mult(wQ, Q);
    mat = mult(R, Q);
    index++;
  }
  return [wQ,mat];
}
/**
 * Inverts a matrix, based on QR decomposition.
 * Will only work properly for symmetric matrices!
 * @param  {Array[Array[Number]]} A A (hopefully symmetric) matrix
 * @return {Array[Array[Number]]}   The inverse of A.
 */
var inverse = function (A) {
  var qrd = QRDecompose(A);
  var evcs = transpose(qrd[0]);
  var det = qrd[1].map((s,i) => s[i]).reduce((a,b) => a*b, 1);
  if (mRound(det)==0)
  {
    var evecs = [], evals = [];
    for (var i=0; i<qrd[1].length; i++)
    {
      if (mRound(qrd[1][i][i],6) != 0)
      {
        evecs.push(evcs[i]);
        evals.push(qrd[1][i][i]);
      }
    }
    Q = transpose(evecs), D = diagonal(evals.map(s => 1/s));
  }
  else
  {
    var Q = qrd[0], D = qrd[1].map((s,i) => s.map((el,j) => (i==j)? 1/el: 0));
  }
  return mult(Q, mult(D, transpose(Q)));
}

// Stuff from RV
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
    var exp = mRound(rvs.reduce((a,b,i) => a + coeffs[i]*b.exp, 0),7);
    var sigma = mRound(rvs.reduce((a,b,i) => a + coeffs[i]*calculate_covariance(b,rvs,coeffs), 0),7);
    this.nm = name;
    this.exp = exp;
    this.cov = {};
    this.cov[name] = mRound(sigma,7);
    for (var i=0; i<rvs.length; i++) {
      var s = rvs[i];
      this.setCov(s,mRound(calculate_covariance(s,rvs,coeffs),7));
    }
    if (args.independent !== undefined)
    {
      var ind = args.independent;
      for (var i=0; i<ind.length; i++)
        this.setCov(ind[i],mRound(calculate_covariance(ind[i],rvs,coeffs),7));
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
  var constRV = new rV('temp', {'mu': -1*this.exp, 'sigma': 0});
  var totalRV = [this].concat(others);
  totalRV.forEach( function (s) {
    constRV.setCov(s,0);
  });
  var sigma = this.cov[this.nm];
  var result = new rV(`S(${this.nm})`, {'rvs': [this, constRV], 'coeffs': [1/Math.sqrt(sigma), 1/Math.sqrt(sigma)], 'independent': others});
  ([result].concat(totalRV)).forEach( function(s) {
    delete s.cov['temp'];
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

/**
 * Calculates the expectation correction for given beliefs and data
 * @param  {Array[rV]} beliefs     The set of belief quantities B
 * @param  {Array[rV]} data        The set of data quantities D
 * @param  {Array[Number]} observation The set of observations for D
 * @param {Array[Array[Number]]} [invmatrix=null] An optional inverse matrix (if multiple such adjustements needed, then this saves computational time)
 * @return {Array[Number]}             The corrected values.
 */
var expectationCorrection = function (beliefs,data,observation,invmatrix=null) {
  var bexp = beliefs.map(s=>[s.exp]), dexp = data.map(s=>[s.exp]);
  var bdcov = build_variance_matrix(beliefs,data);
  if (invmatrix === null)
    dvarinv = inverse(build_variance_matrix(data,data));
  else
    dvarinv = invmatrix;
  var correction = mult(bdcov, mult(dvarinv, observation.map((s,i) => [s-dexp[i][0]])));
  return bexp.map((s,i) => mRound(s[0]+correction[i][0],4));
}

/**
 * [description]
 * @param  {Array[rV]} beliefs          The set of belief quantities B
 * @param  {Array[rV]} data             The set of data quantities D
 * @param  {Array[Array[Number]]} [invmatrix=null] An optional inverse matrix (if multiple such adjustements needed, then this saves computational time)
 * @return {Array[Array[Number]]}                  The corrected variances and covariances.
 */
var varianceCorrection = function (beliefs, data, invmatrix=null) {
  var bvar = build_variance_matrix(beliefs,beliefs), bdcov = build_variance_matrix(beliefs,data);
  if (invmatrix === null)
    dvarinv = inverse(build_variance_matrix(data,data));
  else
    dvarinv = invmatrix;
  var correction = mult(bdcov, mult(dvarinv,transpose(bdcov)));
  return bvar.map((s,i) => s.map((t,j) => mRound(t - correction[i][j],4)));
}

//TESTING
var g0 = new rV('G0', {'mu': 4.16, 'sigma': 1.12});
var g2 = new rV('G2', {'mu': 6.25, 'sigma': 2.43});
var d0 = new rV('D0', {'mu': 4.16, 'sigma': 1.12});
var d2 = new rV('D2', {'mu': 6.25, 'sigma': 2.43});
g0.setCov(g2,0.72);
d0.setCov(d2,0.72);
g0.setCov(d0,0.62);
g0.setCov(d2,0.3);
g2.setCov(d0,0.3);
g2.setCov(d2,0.43);
var observations = [5.4,9.8];
console.log(expectationCorrection([g0,g2],[d0,d2],observations));
console.log(varianceCorrection([g0,g2],[d0,d2]));
