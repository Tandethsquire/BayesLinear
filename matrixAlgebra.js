// Import the base functionality
const base = require('./boilerplate');

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
 * Scales a matrix by a constant. I can't remember if this gets used anywhere.
 * @param  {Array[Array[Number]]} A The matrix to scale
 * @param  {Number} l The scaling factor
 * @return {Array[Array[Number]]}   The scaled matrix
 */
var scale = function (A, l) {
  return A.map(s => s.map(elem => l*elem));
}

/**
 * Concatenates two matrices by appending the columns of B to A.
 * @param  {Array[Array[Number]]} A The matrix to append to
 * @param  {Array[Array[Number]]} B The matrix to append
 * @return {Array[Array[Number]]}   The concatenated matrices
 */
var concatenate = function(A, B) {
  if (A.length == 0)
    return B;
  if (B.length == 0)
    return A;
  if (A.length != B.length)
    throw new Error("Cannot concatenate: matrices do not have the same number of rows.");
  var outmat = [];
  for (var i=0; i<A.length; i++)
  {
    var temparr = [];
    for (var j=0; j<A[i].length; j++)
      temparr.push(A[i][j]);
    for (var j=0; j<B[i].length; j++)
      temparr.push(B[i][j]);
    outmat.push(temparr);
  }
  return outmat;
}

/**
 * Finds the trace of a matrix; that is, the sum of its diagonal entries.
 * The matrix must be square.
 * @param  {Array[Array[Number]]} mat The matrix
 * @return {Number}     The trace
 */
var trace = function (mat) {
  if (mat.length != mat[0].length)
    return TypeError('Matrix trace is not defined on non-square matrices.');
  return mat.map((s,i) => s[i]).reduce((a,b) => a+b,0);
}

/**
 * Checks if a matrix is symmetric. By definition, only works on square matrices.
 * @param  {Array[Array[Number]]} mat The matrix
 * @return {Boolean}
 */
var sym = function (mat) {
  if (mat.length != mat[0].length)
    return false;
  return !(mat.map((s,i) => s.map((el,j) => (Math.abs(el-mat[j][i])>epsilon)).some(x=>x)).some(x=>x));
}

/**
 * Checks if two matrices are the same. Usually used to check if a matrix is the zero matrix using isZero (below).
 * @param  {Array[Array[Number]]} mat1 The first matrix
 * @param  {Array[Array[Number]]} mat2 The second matrix
 * @return {Boolean}
 */
var equal = function (mat1, mat2)
{
  var A = round(mat1), B = round(mat2);
  if (A.length != B.length)
    return false;
  for (var i=0; i<A.length; i++)
  {
    if (A[i].length != B[i].length)
      return false;
    for (var j=0; j<A[i].length; j++)
      if (Math.abs(A[i][j]-B[i][j]) > 0.001)
        return false;
  }
  return true;
}

/**
 * Checks if a matrix equals the zero matrix
 * @param  {Array[Array[Number]]} A The matrix to check
 * @return {Boolean}
 */
var isZero = function (A) {
  var zeroMat = new Array(A.length).fill(0).map(s => new Array(A[0].length).fill(0));
  return equal(A,zeroMat);
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
 * Creates a zero vector of a given length.
 * @param  {Number} size The required vector length
 * @return {Array[Number]}      The output zero vector
 */
var zeroV = function (size) {
  return new Array(size).fill(0).map(s => [0]);
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
 * Normalises a vector: divides every element of the vector by the vector norm
 * @param  {Array[Number]} vec The vector
 * @return {Array[Number]}     The rescaled vector
 */
var vNormalise = function (vec) {
  return vec.map(s => s/vNorm(vec));
}

/**
 * Adds two matrices together. The matrices must have the same dimensions.
 * @param  {Array[Array[Number]]} A The first matrix
 * @param  {Array[Array[Number]]} B The second matrix
 * @return {Array[Array[Number]]}   A+B
 */
var add = function (A, B) {
  if (A.length != B.length || A[0].length != B[0].length)
    return TypeError('Matrices are not the same size.');
  var result = new Array(A.length).fill(0).map(row => new Array(A[0].length).fill(0));
  return result.map((row,i) => row.map((val,j) => A[i][j] + B[i][j]));
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
 * Rounds all elements of a matrix to 7 decimal places.
 * 7dp was chosen because this seems to strike the right balance between accuracy
 * and computational efficiency.
 * @param  {Array[Array[Number]]} A The matrix
 * @return {Array[Array[Number]]}   The rounded matrix
 */
var round = function (A) {
  return A.map(row => row.map(entry => base.mRound(entry,7)));
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
 * Finds the norm of a vector with respect to a non-standard inner product.
 * Typically, the inner product describes a variance matrix, so this effectively
 * finds the variance of a linear combination of quantities under a particular
 * variance specification.
 * @param  {Array[Array[Number]]} v1 The vector (as a matrix)
 * @param  {Array[Array[Number]]} A The inner product
 * @param  {Array[Array[Number]]} v2=null Optional: second vector. If null, assume v1.v1
 * @return {Number}   The length/norm/variance
 */
var innerProd = function (v1, A, v2=null) {
  if (v2 === null)
    return base.mRound(mult(transpose(v1),mult(A,v1))[0][0],7);
  else
    return base.mRound(mult(transpose(v1),mult(A,v2))[0][0],7);
}

/**
 * Creates a Householder matrix from a given vector.
 * @param  {Array[Number]} vec The starting vector
 * @return {Array[Array[Number]]}     The output transformation matrix
 */
var make_householder = function (vec) {
  var tempv = vec.map(x => [x/(vec[0]+base.sign(vec[0])*vNorm(vec))]), H = identity(vec.length);
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
 * Sorts eigenvalues and eigenvectors of a matrix. Usually these were generated by QRDecompose.
 * The two matrices must be organised in the standard way: column i of evecs has eigenvalue given by
 * the i-th diagonal entry of evals.
 * @param  {Array[Array[Number]]} evecs A matrix of eigenvectors
 * @param  {Array[Array[Number]]} evals A diagonal matrix of eigenvalues.
 * @return {Dict}       A collected dictionary, sorted in descending order of eigenvalue.
 */
var eigens = function (evecs, evals) {
  var evcs = evecs, evls = round(evals), dict = [];
  for (var i=0; i<evals.length; i++)
  {
    dict.push({'eval': evls[i][i], 'evec': new Array(evcs[i].length).fill(0).map((s,j) => evcs[i][j])});
  }
  return dict.sort((a,b) => (a.eval < b.eval) ? 1 : -1);
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
  if (base.mRound(det)==0)
  {
    var evecs = [], evals = [];
    for (var i=0; i<qrd[1].length; i++)
    {
      if (base.mRound(qrd[1][i][i],6) != 0)
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

/**
 * Calculates the rank of a matrix: finds the eigenvalues of the matrix, removes
 * repeats (this is not necesarily accurate...) and returns the length of the
 * resulting array, less one if 0 is an element of the array.
 * @param  {Array[Array[Number]]} A The matrix whose rank we wish to calculate
 * @return {Number}   The matrix rank.
 */
var rank = function (A) {
  var evals = QRDecompose(A)[1].map((s,i) => base.mRound(s[i],7));
  return evals.reduce((a,b) => (b!=0)? a+1: a, 0);
}

// Exports for use elsewhere.
module.exports = {add,concatenate,diagonal,eigens,innerProd,inverse,mult,QRDecompose,round,trace,transpose,vNorm,vNormalise,sym,equal,isZero,scale, rank};
