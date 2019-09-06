const base = require('./boilerplate');
const mat = require('./matrixAlgebra');
const RV = require('./RV.js');
const fs = require('fs');

/**
 * Checks if a set Y is Bayes linear sufficient for the separation of sets X and Z.
 * This condition holds if Var(Z,X)=Cov(Z,Y)Var(Y)^(-1)Cov(Y,X).
 * @param  {Array[rV]} X A collection of quantities
 * @param  {Array[rV]} Y The collection of quantities that may separate X and Z
 * @param  {Array[rV]} Z A collection of quantities
 * @return {Boolean}
 */
function separation(X,Y,Z) {
  var lhs = RV.build_variance_matrix(Z,X);
  var rhs = mat.mult(RV.build_variance_matrix(Z,Y), mat.mult(mat.inverse(RV.build_variance_matrix(Y,Y)), RV.build_variance_matrix(Y,X)));
  return mat.equal(lhs,rhs);
}

/**
 * Calculates the resolution transform matrix for given beliefs and data.
 * T = Var(B)^(-1)Cov(B,D)Var(D)^(-1)Cov(D,B).
 * If we want a symmetrised version (perhaps for finding evecs), then set
 * sym=true. Be aware that we then have to transform back at the end of whatever
 * we're doing!
 * @param  {Array[rV]}  B           Beliefs
 * @param  {Array[rV]}  D           Data quantities
 * @param  {Boolean} [sym=false] Do we want a symmetrised version?
 * @return {Array[Array[Number]]}              The resolution transform
 */
var resolution_transform = function (B, D, sym=false) {
  var bVar = RV.build_variance_matrix(B,B), bdVar = RV.build_variance_matrix(B,D), dVar = RV.build_variance_matrix(D,D);
  var T = mat.mult(mat.inverse(bVar), mat.mult(bdVar, mat.mult(mat.inverse(dVar), mat.transpose(bdVar))));
  if (sym) {
    var qrb = mat.QRDecompose(bVar);
    var evecs = qrb[0], evals = qrb[1];
    var d1 = evals.map((s,i) => s.map((el, j) => (i==j)? Math.sqrt(el): 0));
    var d2 = evals.map((s,i) => s.map((el, j) => (i==j && el!=0)? 1/Math.sqrt(el): 0));
    T = mat.mult(d1,mat.mult(mat.transpose(evecs), mat.mult(T, mat.mult(evecs,d2))));
  }
  return T;
}

/**
 * Similar to above, but partial resolution. Effectively B adjusted by D2 in the absence of D1.
 * @param  {Array[rV]}  B           Beliefs
 * @param  {Array[rV]}  D1          The first set of data
 * @param  {Array[rV]}  D2          The second set of data
 * @param  {Boolean} [sym=false] Symmetrised?
 * @return {Arra[Array[Number]]}              The partial res trefo.
 */
var partial_resolution_transform = function(B, D1, D2, sym=false) {
  var bVar = RV.build_variance_matrix(B,B), done = RV.resolution_transform(B, D1), dboth = RV.resolution_transform(B, D1.concat(D2));
  var T = mat.mult(mat.inverse(bVar), mat.mult(bVar, mat.add(dboth, mat.scale(done,-1))));
  if (sym) {
    var qrb = mat.QRDecompose(bVar);
    var evecs = qrb[0], evals = qrb[1];
    var d1 = evals.map((s,i) => s.map((el, j) => (i==j)? Math.sqrt(el): 0));
    var d2 = evals.map((s,i) => s.map((el, j) => (i==j && el!=0)? 1/Math.sqrt(el): 0));
    T = mat.mult(d1,mat.mult(mat.transpose(evecs), mat.mult(T, mat.mult(evecs,d2))));
  }
  return T;
}

/**
 * Finds the heart of the transform for a given resolution transform.
 * Given dependent quantities, we split them into 'meaningful' quantities D
 * and 'less meaningful' quantities E (think coefficients of a linear
 * regression vs. noise terms). Then there is a subspace of B that is most
 * affected by the meaningful parameters.
 * @param  {Array[rV]} B Beliefs
 * @param  {Array[rV]} D The meaningful data quantities
 * @param  {Array[rV]} E Less meaningful quantities
 * @return {Array[Array[Array[Number]]]}   The meaningful space, and its orth compl.
 */
var heart_of_transform = function (B,D,E) {
  var varB = RV.build_variance_matrix(B, B), bEspace = mat.QRDecompose(varB);
  // This is why we need a symmetrised resolution transform.
  var tSymEspace = mat.QRDecompose(resolution_transform(B,D,true));
  var bPsi = mat.diagonal(bEspace[1].map((s,i) => (s[i]!=0) ? 1/Math.sqrt(s[i]) : 0));
  var tEspace = mat.mult(bEspace[0], mat.mult(bPsi, tSymEspace[0]));
  var wPlus = [], wZero = [];
  for (var i=0; i<tEspace.length; i++)
  {
    var coeffs = mat.transpose(tEspace)[i];
    var rv = new RV.rV(`W${i+1}`, {'rvs': B, 'coeffs': coeffs, 'independent': D.concat(E)});
    (mat.round(tSymEspace[1])[i][i] != 0) ? wPlus.push(rv): wZero.push(rv);
  }
  return [wPlus, wZero];
}

/**
 * Checks how much uncertainty is resolved by a given transform: namely the
 * trace of T (above) divided by the rank of Var(B).
 * The rank calculation might go squiffy at times, because of mat.rank().
 * @param  {Array[rV]} B Beliefs
 * @param  {Array[rV]} D Data quantities
 * @return {Number}   Resolution
 */
var uncertainty_resolution = function (B,D) {
  var rk = mat.rank(RV.build_variance_matrix(B,B));
  return base.mRound(mat.trace(resolution_transform(B,D))/rk,4);
}

/**
 * Given beliefs, data quantities, and observation, construct the adjustment bearing:
 * this is the linear combination of elements of B that experiences the largest
 * variance reduction.
 * @param  {Array[rV]} B Beliefs
 * @param  {Array[rV]} D Data quantities
 * @param  {Array[rV]} d Data observation
 * @return {rV}   The bearing Z.
 */
var adjustment_bearing = function (B, D, d) {
  var bMu = B.map(s => s.exp), dMu = D.map(s => s.exp);
  var edb = mat.mult(RV.build_variance_matrix(B,D), mat.mult(mat.inverse(RV.build_variance_matrix(D,D)),d.map((s,i) => [s - dMu[i]])));
  var opMat = mat.round(mat.mult(mat.transpose(edb), mat.inverse(RV.build_variance_matrix(B,B))));
  var bearingCoeffs = opMat[0], bearingAdjustment = -1*mat.mult(opMat, bMu.map(s=>[s]))[0][0];
  var bearingName = B.reduce((a,b,i) => a + ((base.mRound(bearingCoeffs[i],2)==0)? "" : (((a=="")? "": "+") + base.mRound(bearingCoeffs[i],2) + "*" + b.nm)), "").replace(/\+\-/g,"-");
  return new RV.rV(bearingName, {'rvs': B, 'coeffs': bearingCoeffs, 'independent': D.concat(B)});
}

/**
 * Finds path correlation between two adjustments. It returns a number in the range
 * [-1,1]: the closer to 1 (-1), the more complementary (contradictory) the two data
 * sources are. Very much a work in progress.
 * @param  {Array[rV]} dependent The beliefs we're adjusting
 * @param  {Array[rV]} D1        The first set of data quantities
 * @param  {Array[Number]} d1        The corresponding observations of D1
 * @param  {Array[rV]} D2        The second set of data quantities
 * @param  {Array[Number]} d2        The corresponding observations of D2
 * @return {Number}           The correlation
 */
function pathCorrelation (dependent, D1, d1, D2, d2) {
  /**
   * A helper function to unpack coefficients from a name of an adjustment bearing.
   * I'm not super-keen on this, as it's dependent on regEx and any change in the
   * specification of adjustment_bearing screws this up.
   * @param  {rV} rV The adjustment bearing
   * @return {Dict}    The dictionary of coefficients
   */
  var coeffsFromName = function (rV) {
    var coeffsdict = {};
    Object.keys(rV.cov).forEach(function(s,i){
      if (i != 0) {
        var regEx = new RegExp("(\\+?\\-?[0-9\.]*)\\*" + s.replace(/\(/,"\\(").replace(/\)/,"\\)"));
        var mtch = rV.nm.match(regEx);
        if (mtch !== null)
          coeffsdict[s] = parseFloat(mtch[1]);
      }
    });
    return coeffsdict;
  }

  var z1 = adjustment_bearing(dependent, D1, d1), z2 = adjustment_bearing(dependent, D1.concat(D2), d1.concat(d2));
  var a1 = coeffsFromName(z1), a2 = coeffsFromName(z2);
  var v1 = [], v2 =[];
  for (var i=0; i<dependent.length; i++)
  {
    var varnm = dependent[i].nm;
    (a1[varnm] !== undefined)? v1.push([a1[varnm]]): v1.push([0]);
    (a2[varnm] !== undefined)? v2.push([a2[varnm]]): v2.push([0]);
  }
  var varMat = RV.build_variance_matrix(dependent,dependent);
  var vp = v2.map((s,i) => [s[0] - v1[i][0]]);
  return base.mRound(mat.innerProd(v1,varMat,vp)/Math.sqrt(mat.innerProd(v1,varMat,v1)*mat.innerProd(vp,varMat,vp)),4);
}

/**
 * Combines random variables into one list, either by pushing all elements into
 * one big list, or concatenating lists.
 * @param  {rV|Array[rV]} rvlist The list of random variables.
 * @return {Array[rV]}        The output array
 */
function flatten(rvlist) {
  if (rvlist instanceof RV.rV)
    return [rvlist];
  var outarr = [];
  for (var i=0; i<rvlist.length; i++)
    (rvlist[i] instanceof RV.rV) ? outarr.push(rvlist[i]) : outarr = outarr.concat(rvlist[i]);
  return outarr;
}

/**
* Creates a dictionary with all the information needed to build a node diagram.
* Given information about the ordering and parentage of the variables, adds a parent
* entry to each dictionary entry, and adds a child entry to the corresponding parent.
* It then works out the uncertainty resolution each parent provides its child, and works out
* the information that leaves P to C and the information that arrives at C from P.
* These two need not be the same!
* For each object, we get a node with a list of uncertainty resolutions from parents;
* for each parent-child relationship, we get a link with the info arriving and leaving.
 * @param  {Array[String,rV,Array[Number]]} parentData A name, object, and parent index list for each element
 * @param  {Array[rV]} rvs        The random variables in question. The ordering of these is the same as that in parentData.
 * @return {Dict}            A dictionary of nodes and links.
 */
function build_node_list (parentData, rvs) {
  var data = {'nodes': [], 'links': []};
  for (var i=0; i<parentData.length; i++) {
    var child = parentData[i];
    var cObj = flatten(child[1]);
    data.nodes.push({'id': i+1, 'name': child[0], 'parents': []});
    for (var j=0; j<child[2].length; j++) {
      var oldRes = data.nodes[i].parents.reduce((a,b) => a+b.w, 0);
      var parents = flatten(child[2].slice(0,j+1).map(s => rvs[s-1]));
      var res = base.mRound(uncertainty_resolution(cObj,parents)-oldRes, 4);
      data.nodes[i].parents.push({'id': child[2][j], 'w': res});
      var leaving = base.mRound(uncertainty_resolution(cObj, flatten(rvs[child[2][j]-1])),4);
      if (child[2].length ==1) var arriving = leaving;
      else {
        var totalParents = flatten(child[2].map(s => rvs[s-1]));
        var parentsLessOne = flatten(child[2].slice(0,j).concat(child[2].slice(j+1)).map(s => rvs[s-1]));
        var arriving = base.mRound(uncertainty_resolution(cObj,totalParents)-uncertainty_resolution(cObj,parentsLessOne),4);
      }
      data.links.push({'source': child[2][j], 'target': data.nodes[i].id, 'leaving': leaving, 'arriving': arriving});
    }
  }
  return data;
}

// TESTING
// var xvals = [-0.1656, -0.1386, -0.1216, -0.0776, -0.0396, -0.0166, -0.0026, 0.0384, 0.0774, 0.1114, 0.1564, 0.1784];
// var yvals = [63.7, 59.5, 67.9, 68.8, 66.1, 70.4, 70, 73.7, 74.1, 79.6, 77.1, 82.8];
// var zvals = [20.3, 24.2, 18, 20.5, 20.1, 17.5, 18.2, 15.4, 17.8, 13.3, 16.7, 14.8];
// var a = new RV.rV('a', {'mu': 75, 'sigma': 4});
// var b = new RV.rV('b', {'mu': 40, 'sigma': 225});
// var c = new RV.rV('c', {'mu': 20, 'sigma': 1});
// var d = new RV.rV('d', {'mu': -30, 'sigma': 144});
// a.setCov(b,-6);
// a.setCov(c,-1);
// a.setCov(d,0);
// b.setCov(c,0);
// b.setCov(d,-90);
// c.setCov(d,-2.4);
// var varArr = [a,b,c,d];
// var e = [], f = [];
// for (var i=1; i<13; i++)
// {
//   e.push(new RV.rV(`e${i}`, {'mu': 0, 'sigma': 6.25}));
//   f.push(new RV.rV(`f${i}`, {'mu': 0, 'sigma': 4}));
// }
// for (var i=0; i<e.length; i++)
// {
//   varArr.forEach( function(s) {
//     s.setCov(e[i],0);
//   });
//   for (var j=0; j<f.length; j++)
//   {
//     e[j].setCov(e[i], (i==j)? 6.25: 0);
//     f[j].setCov(e[i], (i==j)? 2.5: 0);
//   }
// }
// for (var i=0; i<f.length; i++)
// {
//   varArr.forEach( function(s) {
//     s.setCov(f[i],0);
//   });
//   for (var j=0; j<e.length; j++)
//   {
//     f[j].setCov(f[i], (i==j)? 4: 0);
//     e[j].setCov(f[i], (i==j)? 2.5: 0);
//   }
// }
// var y = [], z = [];
// for (var i=0; i<12; i++)
// {
//   y.push(new RV.rV(`Y${i+1}`, {'rvs': [a,b,e[i]], 'coeffs': [1,xvals[i],1], 'independent': varArr.concat(e).concat(f).concat(y)}));
// }
// for (var i=0; i<12; i++)
// {
//   z.push(new RV.rV(`Z${i+1}`, {'rvs': [c,d,f[i]], 'coeffs': [1,xvals[i],1], 'independent': varArr.concat(e).concat(f).concat(y).concat(z)}));
// }
//
// // var stArr = []
// // var l = varArr.length;
// // for (var i=0; i<l; i++)
// // {
// //   var temp = varArr[i].standardise(varArr.concat(e).concat(f).concat(y).concat(z));
// //   varArr.push(temp);
// //   stArr.push(temp);
// // }
//
// var parentData = [
//   ['E', e, []],
//   ['F', f, [1]],
//   ['a', a, []],
//   ['b', b, [3]],
//   ['c', c, [3,4]],
//   ['d', d, [3,4,5]],
//   ['Y', y, [3,4,1]],
//   ['Z', z, [5,6,2]]
// ]
// console.log(build_node_list(parentData,[e,f,a,b,c,d,y,z]));
