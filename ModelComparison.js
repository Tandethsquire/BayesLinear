const base = require('./boilerplate.js');
const mat = require('./matrixAlgebra.js')
const RV = require('./RV.js');
const comp = require('./varianceComparison.js');

var n = 13;
var mus = new Array(n).fill(0).map(s => base.mRound(Math.random()*10-5,4));
rwalkvar = 0.75;
autovar = 0.1;
var alpha = 0.5;
// Random Walk
var xarr = [];
var temparr = [];
var walkv = new RV.rV('Z', {'mu': 0, 'sigma': rwalkvar});
for (var i=0; i<n; i++) {
  if (i==0) {
    xarr.push(new RV.rV(`X${i}`,{'mu': mus[i], 'sigma': 0}));
    xarr[i].setCov(walkv,0);
  }
  else {
    var temp = new RV.rV(`temp${i}`, {'mu': mus[i]-mus[i-1], 'sigma': 0});
    temp.setCov(walkv,0);
    for (var j=0; j<xarr.length; j++)
      temp.setCov(xarr[j],0);
    temparr.push(temp);
    xarr.push(new RV.rV(`X${i}`,{'rvs': [xarr[i-1], temp, walkv], 'coeffs': [1,1,1], 'independent': xarr.concat(temparr)}));
  }
}
for (var i=0; i<xarr.length; i++)
{
  delete xarr[i].cov['Z'];
  for (var j=1; j<temparr.length+1; j++)
    delete xarr[i].cov[`temp${j}`];
}
var random_walk_variance = RV.build_variance_matrix(xarr.slice(1),xarr.slice(1));

// Autocorrelation
var yarr = [];
var varr = [];
var temparr = [];
var autow = new RV.rV('W', {'mu': 0, 'sigma': autovar});
for (var i=0; i<n; i++) {
  if (i==0) {
    yarr.push(new RV.rV(`Y${i}`, {'mu': mus[i], 'sigma': 0}));
    varr.push(new RV.rV(`V${i}`, {'mu': 0, 'sigma': 0}));
    yarr[i].setCov(autow, 0);
    varr[i].setCov(autow, 0);
    varr[i].setCov(yarr[i],0);
  }
  else {
    var temp = new RV.rV(`temp${i}`, {'mu': mus[i]-mus[i-1], 'sigma': 0});
    temp.setCov(autow,0);
    for (var j=0; j<yarr.length; j++)
    {
      temp.setCov(yarr[j],0);
      temp.setCov(varr[j],0);
    }
    temparr.push(temp);
    varr.push(new RV.rV(`V${i}`, {'rvs': [varr[i-1], autow], 'coeffs': [alpha, 1], 'independent': varr.concat(yarr).concat(temparr)}));
    yarr.push(new RV.rV(`Y${i}`, {'rvs': [yarr[i-1],temp,varr[i]], 'coeffs': [1,1,1], 'independent': varr.concat(yarr).concat(temparr).concat([autow])}));
  }
}
for (i=0; i<yarr.length; i++)
{
  //delete yarr[i].cov['W'];
  for (var j=1; j<temparr.length+1; j++)
    delete yarr[i].cov[`temp${j}`];
  for (var j=0; j<varr.length; j++)
    delete yarr[i].cov[`V${j}`];
}

autocorrelation_variance = RV.build_variance_matrix(yarr.slice(1),yarr.slice(1));

//TESTING
// console.log(comp.arrangeCanonical(comp.canonicalQuantities(random_walk_variance,autocorrelation_variance),random_walk_variance,autocorrelation_variance,mus,mus,mus));
//
console.log(mat.round(mat.QRDecompose(autocorrelation_variance)[1]))
//console.log(comp.arrangeCanonical(comp.canonicalQuantities(random_walk_variance,autocorrelation_variance)));
