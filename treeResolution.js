// Imports from other packages
const base = require('./boilerplate');
const mat = require('./matrixAlgebra');
const fs = require('fs');

// A global array to store random variables in.
// This is so that when a new variable is generated, we can find covariances
// between all previously defined quantities.
globalArr = [];

/**
 * The class definition for a random variable.
 * @param  {String} name        The name of the random variable
 * @param  {Number} mu          The expectation of the RV
 * @param  {Number} sigma       The variance of the RV (NOT the standard deviation)
 * @param  {String|null} [vstr=null] The 'variance string'. See the function comments.
 * @return {rV}             The random variable object.
 */
var rV = function (name, mu, sigma, vstr=null) {
  if (isNaN(sigma) || isNaN(mu))
    throw new Error('Either the expectation or variance are non-numeric');
  this.nm = name;
  // If we want to create a new rV by translating an old one, we don't want the
  // constant to affect the variance. We use vstr to keep only the parts of the
  // quantity that affect the variance calculations.
  (vstr !== null) ? this.vs = vstr : this.vs = name;
  this.exp = mu;
  // This creates a dictionary of covariances between the rV and others; it starts
  // by setting Var(rV)=sigma.
  this.cov = {};
  this.cov[name] = sigma;
}
/**
 * Sets the covariance between two rVs, by appending entries to thier cov dicts.
 * @param  {rV}  that           The other rV
 * @param  {Number}  cov            The covariance between the two
 * @param  {Boolean} [oneway=false] If true, we only add the covariance to this.cov
 */
rV.prototype.setCov = function(that, cov, oneway = false) {
  if (!(that instanceof rV))
    throw new Error('Cannot set a covariance between a random variable and a non-random one');
  this.cov[that.nm] = cov;
  if (!oneway)
    that.cov[this.nm] = cov;
}
/**
 * Adds to a rV. If we're adding another rV to it, then both the expectation and
 * the variance are recalculated; if we're adding a number, then we change the expectation
 * and use vstr for the purpose to which it is designed.
 * @param  {rV|Number}  that              The quantity to add
 * @param  {String}  [label=undefined] An optional label for naming the result
 * @param  {Boolean} [ephemeral=true]  If this is an intermediate calculation, we use the one-way setCov.
 */
rV.prototype.add = function (that, label = undefined, ephemeral = true) {
  if (that instanceof rV) {
    var expectation = this.exp + that.exp;
    var sigma = this.cov[this.nm] + 2*that.cov[this.nm] + that.cov[that.nm];
    var name = `(${this.nm}+${that.nm})`, vstr = name;
  }
  else if (!isNaN(that))
  {
    var expectation = this.exp + that;
    var sigma = this.cov[this.nm];
    var name = `(${this.nm}+${that.toString()})`, vstr = this.nm;
  }
  else {
    throw new Error('The object you are trying to add is not a variable.')
  }
  var ident;
  (label !== undefined) ? ident = label: ident = name;
  var res = new rV(ident, expectation, sigma, vstr);
  globalArr.forEach(function(elem) {
    res.setCov(elem, calcCov(elem, res.vs), ephemeral);
  });
  return res;
}
/**
 * Multiplies an rV by a constant.
 * @param  {Number}  k                 The constant
 * @param  {String}  [label=undefined] An optional label for naming the result
 * @param  {Boolean} [ephemeral=true]  If this is an intermediate calculation, we use the one-way setCov.
 */
rV.prototype.multiply = function (k, label = undefined, ephemeral = true) {
  var expectation = k * this.exp;
  var sigma = k * k * this.cov[this.nm];
  var name = `${k}*${this.nm}`, ident;
  (label !== undefined) ? ident = label: ident = name;
  var res = new rV(ident, expectation, sigma, k+"*"+this.vs);
  globalArr.forEach(function(elem) {
    res.setCov(elem,calcCov(elem,res.vs), ephemeral);
  });
  return res;
}
/**
 * Creates a standardised version of a rV via S(X)=(X-E(X))/sqrt(Var(X))
 * @return {rV} The standardised random variable.
 */
rV.prototype.standardise = function ()
{
  return this.add(-1*this.exp).multiply(1/Math.sqrt(this.cov[this.nm]),`S${this.nm}`, false);
}

/**
 * Calculates covariance between a random variable and a string representing
 * a linear combination of rVs.
 * @param  {rV} rv  The random variable
 * @param  {String} str The linear combination
 * @return {Number}     The resulting covariance
 */
var calcCov = function (rv, str) {
  if (!(rv instanceof rV) || (typeof str !== 'string'))
    throw new Error('Can\'t calculate covariance.');
  globalArr.forEach(function(elem) {
    var reg = new RegExp("(?<=^|[\\(\\+\\-\\*])"+elem.nm+"(?=[\\)\\+\\-\\*]|$)")
    //var reg = new RegExp(elem.nm+'(?=\\D|$)', 'g');
    str = str.replace(reg,rv.cov[elem.nm]);
  });
  return eval(str);
}

/**
 * Makes a variance matrix from two collections of random variables.
 * If we have A=[a,b,c] and B=[d,e,f], where a,b,...,f are all rVs,
 * then bVM(A,A)=Var(A), bVM(A,B)=Cov(A,B), etc.
 * @param  {Array[rV]} args1 The first set of rVs.
 * @param  {Array[rV]} args2 The second set of rVs.
 * @return {Array[Array[Number]]}       The variance matrix.
 */
var buildVarianceMatrix = function (args1, args2) {
  var outarr = new Array(args1.length).fill(0).map(s => new Array(args2.length).fill(0));
  return outarr.map((s,i) => s.map((elem,j) => args1[i].cov[args2[j].nm]));
}

/**
 * Finds the transformation matrix for B given D.
 * @param  {Array[rV]}  B           The collection of beliefs
 * @param  {Array[rV]}  D           The collection of data
 * @param  {Boolean} [sym=false] Do we want the transformation matrix to be symmetrised (eg for finding evals)?
 * @return {Array[Array[Number]]}              The transform.
 */
var dataTransform = function(B,D,sym = false) {
  var bvar = buildVarianceMatrix(B,B), covar = buildVarianceMatrix(B,D), dvar = buildVarianceMatrix(D,D);
  var T = mat.mult(mat.inverse(bvar), mat.mult(covar, mat.mult(mat.inverse(dvar),mat.transpose(covar))));
  if (sym)
  {
    var qrb = mat.QRDecompose(bvar);
    var evecs = qrb[0], evals = qrb[1];
    var d1 = evals.map((s,i) => s.map((el,j) => (i==j)? Math.sqrt(el): 0));
    var d2 = evals.map((s,i) => s.map((el,j) => (i==j)? 1/Math.sqrt(el): 0));
    T = mat.mult(d1,mat.mult(mat.transpose(evecs), mat.mult(T, mat.mult(evecs,d2))));
  }
  return T;
}

/**
 * Similar to the above, but for partial transforms.
 * This calculates the transform of B given D1, and then calculates the transform
 * of that given D2. Not the same as transforming B given D1 & D2.
 * @param  {Array[rV]}  B           Beliefs
 * @param  {Array[rV]}  D1          Data 1
 * @param  {Array[rV]}  D2          Data 2
 * @param  {Boolean} [sym=false] Does this want to be symmetrised?
 * @return {Array[Array[Number]]}              The transform matrix
 */
var partialDataTransform = function(B,D1,D2,sym=false)
{
  var bvar = buildVarianceMatrix(B,B), dalone = dataTransform(B,D1), dboth = dataTransform(B,D1.concat(D2));
  var T = mat.mult(mat.inverse(bvar),mat.mult(bvar, mat.add(dboth,mat.scale(dalone,-1))));
  if (sym)
  {
    var qrb = mat.QRDecompose(bvar);
    var evecs = qrb[0], evals = qrb[1];
    var d1 = evals.map((s,i) => s.map((el,j) => (i==j)? Math.sqrt(el): 0));
    var d2 = evals.map((s,i) => s.map((el,j) => (i==j)? 1/Math.sqrt(el): 0));
    T = mat.mult(d1,mat.mult(mat.transpose(evecs), mat.mult(T, mat.mult(evecs,d2))));
  }
  return T;
}

/**
 * Works out the uncertainty resolution for a given adjustment.
 * This is calculated as the trace of the transformation matrix, divided by
 * the rank of the belief variance matrix. At present, it does not calculate
 * the rank and instead assumes bVM(B,B) is full rank.
 * @param  {Array[rV]} B Beliefs
 * @param  {Array[rV]} D Data
 * @return {Number}   The resolution.
 */
var uncertaintyResolution = function(B,D) {
  // // Need to calculate rank properly! Also related to calculating inverses properly.
  return base.mRound(mat.trace(dataTransform(B,D))/B.length,4);
}

/**
 * Creates a random variable from a linear combination of other random variables.
 * The coefficients are given by a vector, where the length of the vector is the
 * same as the length of 'dependent'.
 * The 'independent' quantites are those in globalArr that do not directly contribute
 * to the random variable, but for which there may be non-zero covariance.
 * @param  {Array[Number]} vec         The vector of coefficients
 * @param  {Array[rV]} dependent   The directly dependent quantities
 * @param  {Array[rV]} independent The 'independent' quantities (defined above)
 * @param  {String} name        The name to give the result
 * @return {rV}             The random variable
 */
function buildRVFromVector(vec,dependent,independent,name)
{
  var exp = vec.map((s,i) => s*dependent[i].exp).reduce((a,b) => a+b, 0);
  var sigma = mat.innerProd(vec.map(s => [s]),buildVarianceMatrix(dependent,dependent));
  var outputRV = new rV(name,exp,sigma);
  var strform = vec.map((s,i) => s+"*"+dependent[i].nm).reduce((a,b) => a+"+"+b,"");
  independent.forEach(function(elem) {elem.setCov(outputRV,calcCov(elem,strform))});
  return outputRV;
}

/**
 * Finds the bearing of an adjustment; this is the linear combination of elements
 * of B that experiences the largest reduction in variance.
 * @param  {Array[rV]} B The beliefs.
 * @param  {Array[rV]} D The data.
 * @param  {Array[Number]} d The observed data
 * @return {[Array[Number],Number,rV]}   The coefficients of the bearing; the size of the adjustment, and the bearing as an rV
 */
function adjustmentBearing(B,D,d) {
  var beliefExp = B.map(s => s.exp), dataExp = D.map(s => s.exp);
  var edb = mat.mult(buildVarianceMatrix(B,D),mat.mult(mat.inverse(buildVarianceMatrix(D,D)),d.map((s,i) => [s - dataExp[i]])));
  var opmat = mat.round(mat.mult(mat.transpose(edb),mat.inverse(buildVarianceMatrix(B,B))));
  var bng = buildRVFromVector(opmat[0],B,B.concat(D),"Zd"), adj = -1*mat.mult(opmat,beliefExp.map(s=>[s]))[0][0];
  globalArr.push(bng);
  return [opmat[0],adj,bng.add(adj)];
}

/**
 * Finds the heart of the transform; this is the subspace of the full space that has non-zero
 * resolution with respect to the meaningful parameters D.
 * @param  {Array[rV]} B Beliefs
 * @param  {Array[rV]} D Meaningful parameters
 * @param  {Array[rV]} E 'Meaningless' parameters (eg error components)
 * @return {[Array[Array[Number]],Array[Array[Number]]]}   The meaningful space, and its orthogonal complement
 */
function heartOfTransform(B, D, E)
{
  var varB = buildVarianceMatrix(B, B), bEv = mat.QRDecompose(varB), tSymEv = mat.QRDecompose(dataTransform(B,D,true));
  var bPsi = mat.diagonal(bEv[1].map((s,i) => (s[i]!=0) ? 1/Math.sqrt(s[i]) : 0));
  var tEv = mat.mult(bEv[0], mat.mult(bPsi, tSymEv[0]));
  var Wplus = [], Wzero = [];
  for (var i=0; i<tEv.length; i++)
  {
    var coeffs = mat.transpose(tEv)[i];
    var rv = buildRVFromVector(coeffs, B, D.concat(E), `W${i+1}`);
    (mat.round(tSymEv[1])[i][i] != 0) ? Wplus.push(rv) : Wzero.push(rv);
  }
  return [Wplus, Wzero];
}

/**
 * Checks if Y is Bayes linear sufficent for the separation of X and Z.
 * @param  {Array[rV]} X A collection of quantities
 * @param  {Array[rV]} Y The collection of quantities that may separate X and Z
 * @param  {Array[rV]} Z A collection of quantities
 * @return {Boolean}
 */
function separation(X, Y, Z)
{
  var lhs = buildVarianceMatrix(Z, X);
  var rhs = mat.mult(buildVarianceMatrix(Z,Y), mat.mult(mat.inverse(buildVarianceMatrix(Y,Y)), buildVarianceMatrix(Y,X)));
  return mat.equal(lhs,rhs);
}

/**
 * A generic function for combining lists of rVs.
 * For each element of rvlist, we either push it into an array (if it's a rV),
 * or we concatenate the array with it (if it's an array of rVs.)
 * @param  {Array[Array[rV]|rV]} rvlist The elements to combine
 * @return {Array[rV]}        The resulting cleaned array
 */
function combineRV(rvlist)
{
  if (rvlist instanceof rV)
    return [rvlist];
  var outarr = [];
  for (var i=0; i<rvlist.length; i++)
    (rvlist[i] instanceof rV) ? outarr.push(rvlist[i]) : outarr = outarr.concat(rvlist[i]);
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
 * @param  {Array[rV]} objects    A list of objects, ordered in the same order as parentData
 * @return {Dict}            The dictionary of links and nodes
 */
function buildNodeList (parentData, objects) {
  var data = {'nodes': [], 'links': []}
  for (var i=0; i<parentData.length; i++)
  {
    var child = parentData[i];
    var cObj = combineRV(child[1]);
    data.nodes.push({'id': i+1, 'name': child[0], 'parents': []});
    for (var j=0; j<child[2].length; j++)
    {
      var oldRes = data.nodes[i].parents.reduce((a,b) => a + b.w, 0)
      var parents = combineRV(child[2].slice(0,j+1).map(s => objects[s-1]));
      var res = base.mRound(uncertaintyResolution(cObj, parents) - oldRes, 4);
      data.nodes[i].parents.push({'id': child[2][j], 'w': res});
      var leaving = base.mRound(uncertaintyResolution(cObj,combineRV(objects[child[2][j]-1])),4);
      if (child[2].length == 1)
        var arriving = leaving;
      else
      {
        var totalParents = combineRV(child[2].map(s => objects[s-1]));
        var parentsLessOne = combineRV(child[2].slice(0,j).concat(child[2].slice(j+1)).map(s => objects[s-1]));
        var arriving = base.mRound(uncertaintyResolution(cObj,totalParents) - uncertaintyResolution(cObj,parentsLessOne),4);
      }
      data.links.push({'source': child[2][j], 'target': data.nodes[i].id, 'leaving': leaving, 'arriving': arriving});
    }
  }
  return data;
}

// TESTING: Random Variables
// Y = a + b xi +ei
// Z = c + d xi +fi
// i runs from 1 to 12
// Defining a, b, c, d
// var xvals = [-0.1656, -0.1386, -0.1216, -0.0776, -0.0396, -0.0166, -0.0026, 0.0384, 0.0774, 0.1114, 0.1564, 0.1784];
// var yvals = [63.7, 59.5, 67.9, 68.8, 66.1, 70.4, 70, 73.7, 74.1, 79.6, 77.1, 82.8];
// var zvals = [20.3, 24.2, 18, 20.5, 20.1, 17.5, 18.2, 15.4, 17.8, 13.3, 16.7, 14.8];
// var a = new rV('a', 75, 4);
// var b = new rV('b', 40, 225);
// var c = new rV('c', 20, 1);
// var d = new rV('d', -30, 144);
// // Creating the set E = [e1,e2,...] and F = [f1,f2,...]
// var e = [], f = [];
// for (var i=1; i<13; i++)
// {
//   e.push(new rV(`e${i}`,0,6.25));
//   f.push(new rV(`f${i}`,0,4));
// }
// // Setting the 'artifical' covariances by hand
// a.setCov(b,-6);
// a.setCov(c,-1);
// a.setCov(d,0);
// b.setCov(c,0);
// b.setCov(d,-90);
// c.setCov(d,-2.4);
// var globalArr = [a,b,c,d];
// for (var i=0; i<e.length; i++)
// {
//   e[i].setCov(a,0);
//   e[i].setCov(b,0);
//   e[i].setCov(c,0);
//   e[i].setCov(d,0);
//   f[i].setCov(a,0);
//   f[i].setCov(b,0);
//   f[i].setCov(c,0);
//   f[i].setCov(d,0);
//   for (var j=0; j<i; j++)
//   {
//     e[i].setCov(e[j],0);
//     f[i].setCov(f[j],0);
//   }
//   for (var j=0; j<f.length; j++)
//   {
//     (i==j) ? e[i].setCov(f[j],2.5): e[i].setCov(f[j],0);
//   }
//   globalArr.push(e[i]);
//   globalArr.push(f[i]);
// }
// // Creating Y=[Y1,Y2,...] and Z=[Z1,Z2,...]
// var y = [], z = [];
// // The covariances of these are calculated, since Y and Z are linear in previouslt determined quantities
// for (var i=0; i<12; i++)
// {
//   tempy = a.add(e[i].add(b.multiply(xvals[i])),`Y${i+1}`,false);
//   globalArr.push(tempy);
//   tempz = c.add(f[i].add(d.multiply(xvals[i])),`Z${i+1}`,false);
//   globalArr.push(tempz);
//   y.push(tempy);
//   z.push(tempz);
// }
// // Observed values of Y and Z
// var yvals = [63.7,59.5,67.9,68.8,66.1,70.4,70,73.7,74.1,79.6,77.1,82.8], zvals = [20.3,24.2,18,20.5,20.1,17.5,18.2,15.4,17.8,13.3,16.7,14.8];

// Straightforward relations: Y and Z are children
// var parentData = [
//   ['E',e,[]],
//   ['F',f,[1]],
//   ['a',a,[]],
//   ['b',b,[3]],
//   ['c',c,[3,4]],
//   ['d',d,[3,4,5]],
//   ['Y',y,[3,4,1]],
//   ['Z',z,[5,6,2]]
// ]
//
// Reversed so that a,b,c,d,E,F are children of Y, Z, rather than parents.
// var childData = [
//   ['Y', y ,[]],
//   ['Z', z, []],
//   ['a', a, [1,2]],
//   ['b', b, [1,2]],
//   ['c', c, [1,2]],
//   ['d', d, [1,2]],
//   ['E', e, [1,2]],
//   ['F', f, [1,2]]
// ]
//
// The heart of the transform: [a,b,c,d] are the meaningful quantities
// var heart = heartOfTransform(y.concat(z), [a,b,c,d], e.concat(f));
// var heartData = [
//   ['a', a, []],
//   ['b', b, [1]],
//   ['c', c, [1,2]],
//   ['d', d, [1,2,3]],
//   ['E', e, []],
//   ['F', f, [5]],
//   ['W+', heart[0], [1,2,3,4,5,6]],
//   ['W0', heart[1], [5,6]]
// ];
//
// Combining [a,b] and [c,d], motivated by the structure of Y and Z.
// var combinedData = [
//   ['Gy', [a,b], []],
//   ['Gz', [c,d], [1]],
//   ['E', e, []],
//   ['F', f, [3]],
//   ['Y', y, [1,3]],
//   ['Z', z, [2,4]]
// ];
// var outData = buildNodeList(combinedData, [[a,b],[c,d],e,f,y,z]);
// fs.writeFileSync('combineddata.json', JSON.stringify(outData));

Standardisation of [a,b,c,d] to find the standardised bearing
var origarr = [a,b,c,d], starr = [];
for (var i=0; i<origarr.length; i++)
{
  var temp = origarr[i].standardise();
  globalArr.push(temp);
  starr.push(temp);
}

Finding the adjustment with respect to Y&Z
console.log(adjustmentBearing(starr,y.concat(z),yvals.concat(zvals)));

Given data, partial adjustment by Y and then by Z.
var tPart = partialDataTransform(starr,y,z,true);
console.log(mat.QRDecompose(tPart).map(s => mat.round(s)));
