// Import required dependencies
const base = require('./boilerplate');
const mat = require('./matrixAlgebra');
const fs = require('fs');

/**
 * Scales a vector to unit variance under a choice of two variance matrices A and B.
 * If the vector has zero variance under A, we scale under B. If this is also 0,
 * then we scale the vector to have unit length.
 * @param  {Array[Number]} v The vector to normalise
 * @param  {Array[Array[Number]]} A The first variance matrix
 * @param  {Array[Array[Number]]} B The second variance matrix
 * @return {Array[Number]}   The output, scaled, vector.
 */
var varianceScale = function (v, A, B) {
  if (mat.innerProd(v,A) != 0)
    return v.map(s => [s[0]/Math.sqrt(mat.innerProd(v,A))]);
  if (mat.innerProd(v,B) != 0)
    return v.map(s => [s[0]/Math.sqrt(mat.innerProd(v,B))]);
  return mat.vNormalise(v).map(s => [s]);
}

/**
 * Finds the correction to a vector such that its expectation under model one is 0.
 * If the resulting expectation under model two is negative, we then rescale the vector
 * to have positive expectation.
 * @param  {Array[Number]} v  The vector to correct
 * @param  {Array[Number]} E1 The first expectation specification
 * @param  {Array[Number]} E2 The second specification
 * @return {[Array[Number], Number]}    The vector (possibly rescaled), and the correction.
 */
var expectationCorrection = function (v, E1, E2) {
  var expc = -1*v.map((s,i) => s[0]*E1[i]).reduce((a,b) => a+b, 0);
  var exp2 = v.map((s,i) => s[0]*E2[i]).reduce((a,b) => a+b, 0) + expc;
  if (exp2 < 0)
    return [v.map(s => [-1*s[0]]), -1*expc];
  return [v, expc];
}

/**
 * Finds the canonical quantities for a two-specification system (A and B).
 * Applies QR decomposition to A+B; zero-eval eigenvectors are scaled and stored.
 * Then takes the remaining space w.r.t. B and applies QR decomposition.
 * The results are concatenated, giving a full canonical quantity description.
 * @param  {Array[Array[Number]]} A The matrix of the first specification.
 * @param  {Array[Array[Number]]} B The matrix of the second specification.
 * @return {Array[Array[Number]]}   The canonical quantities, expressed as columns of a matrix
 */
var canonicalQuantities = function(A, B) {
  if (A == null)
  {
    A = new Array(B.length).fill(0).map(s => new Array(B.length).fill(0));
  }
  var msum = mat.add(A, B);
  var qR = mat.QRDecompose(msum);
  var ev = mat.eigens(mat.transpose(qR[0]),qR[1]);
  var positive = ev.filter(s => s.eval>0), zero = ev.filter(s => s.eval==0), negative = ev.filter(s => s.eval<0);
  if (negative.length > 0)
    console.log("The matrix system has some negative eigenvalues (their canonical quantities have been omitted): treat the analysis with caution!");
  var W = [], R = [], diag1 = new Array(positive.length).fill(0).map(s => new Array(positive.length).fill(0));
  for (var i=0; i<positive.length; i++) {
    W.push(positive[i].evec);
    diag1[i][i] = positive[i].eval;
  }
  if (zero.length != 0) {
    for (var i=0; i<zero.length; i++)
      R.push(mat.vNormalise(zero[i].evec));
    R = mat.transpose(R);
  }
  W = mat.transpose(W);
  var Theta = diag1.map((row,i) => row.map((elem,j) => (i==j) ? 1/Math.sqrt(elem): 0));
  var G = mat.mult(mat.mult(Theta,mat.mult(mat.transpose(W),mat.mult(B,W))),Theta);
  var gqr = mat.QRDecompose(G);
  var gev = mat.eigens(mat.transpose(gqr[0]),gqr[1]);
  var gnonzero = gev.filter(s => s.eval > 0), gzero = gev.filter(s => s.eval == 0);
  var Y = [], V = [], diag2 = new Array(gnonzero.length).fill(0).map(s => new Array(gnonzero.length).fill(0));
  for (var i=0; i<gnonzero.length; i++) {
    Y.push(gnonzero[i].evec);
    diag2[i][i] = gnonzero[i].eval;
  }
  if (gzero.length != 0) {
    for (var i=0; i<gzero.length; i++)
      V.push(gzero[i].evec);
    V = mat.transpose(V);
  }
  Y = mat.transpose(Y);
  var Delta = diag2.map((row,i) => row.map((elem, j) => (i==j) ? 1/Math.sqrt(elem) : 0));
  var comp1 = mat.mult(W, mat.mult(Theta, mat.mult(Y,Delta)));
  if (V.length != 0)
    var comp2 = mat.mult(W, mat.mult(Theta, V));
  else
    var comp2 = [];
  var comp3 = R;
  return mat.round(mat.concatenate(mat.concatenate(comp1,comp2),comp3));
}

/**
 * Takes a canonical quantity, and its index, and formats it in a human-readable way.
 * @param  {Array[Number]}  v                The canonical quantity
 * @param  {Number}  i                The index (i.e. which canonical quantity is it)
 * @param  {Number}  [ecorr=0]        Any expectation correction; 0 by default
 * @param  {Boolean} [wantLabel=true] Should a label be prefixed to the string; true by default
 * @return {String}                   The output string
 */
var formatCQ = function (v,i,ecorr=0,wantLabel=true) {
  var prestr = "";
  if (wantLabel)
    prestr = `Z${i+1}=`;
  return `${prestr}${v.map((elem, j) => (Math.abs(elem[0])<0.000001)? "": `${base.mRound(elem[0],6)}*X${j+1}`).join('+')}+${(Math.abs(ecorr)>0)? ("+" + ecorr): ""}`.replace(/\++/g,"+").replace(/\+\-/g,"-").replace(/\=\+/,"=").replace(/(?:^\+|\+$)/g,"");
}

/**
 * Creates a dictionary of canonical quantities.
 * Dictionary entries are: human readable quantity; coefficients of the CQ in vector form;
 * expectation correction; variances under the two specifications' expectations under the two specifications.
 * The entries are sorted as follows:
 * Quantities with non-zero variance under H1 (which is thus scaled to 1) are first, sorted in descending order
 * of their variance under H2;
 * Quantities with zero variance under H1 (and so variance 1 under H2) are next, in the order provided;
 * Quantities with zero variance under both specifications are last.
 * @param  {Array[Array[Number]]} cMat The matrix of canonical quantities
 * @param  {Array[Array[Number]]} A    The first variance specification
 * @param  {Array[Array[Number]]} B    The second variance specification
 * @param  {Array[Number]} E1   The first expectation specification
 * @param  {Array[Number]} E2   The second expectation specification
 * @return {Dict}      A dictionary of canonical quantities, with vital statistics.
 */
var arrangeCanonical = function (cMat, A, B, E1, E2) {
  var dict = [], wMat = mat.transpose(cMat);
  for (var i=0; i<wMat.length; i++)
  {
    var cq = varianceScale(wMat[i].map(s => [s]), A, B);
    var corr = expectationCorrection(cq, E1, E2);
    cq = corr[0];
    var ecorr = corr[1];
    var V1 = mat.round(mat.mult(mat.transpose(cq), mat.mult(A, cq)))[0][0];
    var V2 = mat.round(mat.mult(mat.transpose(cq), mat.mult(B, cq)))[0][0];
    dict.push({'quantity': formatCQ(cq,i,ecorr,false), 'vectorcoeffs': cq.map(s => s[0]), 'correction': ecorr, 'var1': V1, 'var2': V2, 'exp1': 0, 'exp2': base.mRound(cq.map((s,i) => s[0]*E2[i]).reduce((a,b) => a+b, 0) + ecorr,7)});
  }
  dict.sort( function (a,b) {
    if (a.var1 ==1 && b.var1 == 1)
      return (a.var2 - b.var2);
    if (a.var1 == 0 && b.var1 == 0)
      return (b.var2 - a.var2);
    if (a.var1 == 0 || b.var1 == 0)
      return ((a.var1 == 0)? 1: -1)
  });
  for (var i=0; i<dict.length; i++)
    dict[i].label = `Z${i+1}`;
  return dict;
}

/**
 * A measure of the predictive difference between two specifications.
 * Canonical quantities are collected into three sets: positive variance under H1 and H2,
 * positive under H1 and zero under H2; and vice versa. Then the predictive difference is
 * a particular ratio of these values. Closer to 1 indicates little difference in predictive power.
 * @param  {Dict} dict A dictionary of canonical quantities.
 * @param  {Array[Array[Number]]} A    The variance specification H1
 * @param  {Array[Array[Number]]} B    The variance specification H2
 * @return {Number}      The predictive measure.
 */
var totalVarianceDifference = function (dict, A, B) {
  var a0 = 0; aplus = 0; aplusl = 0; b0 = 0; Y = new Array(A.length).fill(1);
  for (var i=0; i<dict.length; i++)
  {
    var c1 = mat.mult(mat.transpose(dict[i].vectorcoeffs.map(s => [s])),mat.mult(A,Y.map(s=>[s])))[0][0];
    var c2 = mat.mult(mat.transpose(dict[i].vectorcoeffs.map(s => [s])),mat.mult(B,Y.map(s=>[s])))[0][0];
    if (dict[i].var1 != 0)
    {
      if (dict[i].var2 != 0)
      {
        aplus += c1*c1;
        aplusl += dict[i].var2*c1*c1;
      }
      else
        a0 = c1*c1;
    }
    else if (dict[i].var2 != 0)
      b0 += c2*c2;
  }
  return (aplusl + b0)/(aplus + a0);
}

/**
 * Generates residuals of a system given an observation.
 * @param  {Dict} dict A dictionary of canonical quantities.
 * @param  {Array[Number]} obs  The observation of the quantities
 * @return {Dict}      The dictionary, with the residuals added as entries
 */
var observedResiduals = function (dict, obs) {
  for (var i=0; i<dict.length; i++)
  {
    var obz = dict[i].vectorcoeffs.map((s,i) => obs[i]*s).reduce((a,b) => a+b, 0) + dict[i].correction;
    dict[i].obs = base.mRound(obz,5);
    dict[i].r1 = base.mRound(Math.abs(obz),5);
    (dict[i].var2 == 0) ? dict[i].r2 = 0: dict[i].r2 = base.mRound(Math.abs((obz-dict[i].exp2)/Math.sqrt(dict[i].var2)),5);
  }
  return dict;
}

/**
 * Generates the vital stats for the observed comparison G21
 * @param  {Dict} dict The dictionary of canonical quantities
 * @param  {Array[Array[Number]]} A    The first variance specification
 * @param  {Array[Array[Number]]} B    The second varaince specification
 * @param  {Array[Number]} E1   The first expectation specification
 * @param  {Array[Number]} E2   The second expectation specification
 * @param  {Array[Number]} obs  The observed quantities
 * @return {Dict}      A dictionary of the vital statistics for the observed comparison
 */
var standardised = function (dict, A, B, E1, E2, obs) {
  var gVec = new Array(A.length).fill(0), corr = 0;
  for (i=0; i<dict.length; i++)
  {
    var qu = dict[i];
    if (qu.var1 != 0 && qu.var2 != 0)
    {
      gVec = gVec.map((s,j) => s + qu.exp2*qu.vectorcoeffs[j]);
      corr += qu.exp2*qu.correction;
    }
  }
  gObs = gVec.map((s,i) => s*obs[i]).reduce((a,b) => a+b, 0) + corr;
  g1 = gVec.map((s,i) => s*E1[i]).reduce((a,b) => a+b, 0) + corr;
  g2 = gVec.map((s,i) => s*E2[i]).reduce((a,b) => a+b, 0) + corr;
  gVar1 = mat.innerProd(gVec.map(s=>[s]),A);
  gVar2 = mat.innerProd(gVec.map(s=>[s]),B);
  return {'observed': gObs, 'exp1': g1, 'exp2': g2, 'var1': gVar1, 'var2': gVar2};
}

//TESTING: Canonical Quantities
// var U1 = [[2,1,5,1,0,0,0],[1,4,6,-3,0,0,0],[5,6,16,-1,0,0,0],[1,-3,-1,4,0,0,0],[0,0,0,0,6,-4,7],[0,0,0,0,-4,6,-5],[0,0,0,0,7,-5,9]];
// var U2 = [[2,1,3,1,0,0,0],[1,2,3,-1,0,0,0],[3,3,6,0,0,0,0],[1,-1,0,2,0,0,0],[0,0,0,0,6,0,7],[0,0,0,0,0,6,-1],[0,0,0,0,7,-1,9]];
// var E1 = [0,0,0,0,-1,0,0];
// var E2 = [0,0,0,0,1,0,0];
// var X = [3,1,4,2,13,-5,14];
// var can = arrangeCanonical(canonicalQuantities(U1,U2),U1,U2,E1,E2);
// console.log(can);
// console.log(standardised(can,U1,U2,E1,E2,X));
