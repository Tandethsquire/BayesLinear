/**
 * Finds the sign of a number (-1 if negative, 1 if positive, 0 else)
 * @param  {Number} x The number
 * @return {Number}   The sign: either -1, 0, or 1
 */
var sign = function(x) {
  return ((x==0)? 0 :((x < 0) ? -1: 1));
}

/**
 * Copies a matrix and returns the copy, to avoid passing by reference.
 * @param  {Array[Array[Number]]} arr The array to copy
 * @return {Array[Array[Number]]}     A new copy of the same array
 */
var copyArray = function (arr) {
  var result = new Array(arr.length).fill(0).map(s => new Array(arr[0].length).fill(0));
  return result.map((s,i) => s.map((elem,j) => arr[i][j]));
}

/**
 * Checks if an element is in any array. Only works if the type of the element
 * has a well-defined equality measure.
 * @param  {Array[Object]} arr  The array to search in
 * @param  {Object} elem The element to find
 * @return {Boolean}
 */
var inArray = function (arr,elem) {
  return arr.some(s => s==elem);
}

/**
 * Prints a matrix in Mathematica format (mostly for debugging)
 * @param  {Array[Array[Number]]} arr A matrix
 * @return {String}     The output string {{a,b,c,..},{d,e,f,..},..}
 */
var printArray = function (arr) {
  var outstr = "{";
  for (var i=0; i<arr.length; i++)
  {
    (arr[i] instanceof Array) ? outstr += printArray(arr[i]) : outstr += arr[i];
    (i!=arr.length-1) ? outstr += ",": outstr += "";
  }
  return outstr + "}";
}

/**
 * Creates a human-readable string from a matrix.
 * @param  {Array[Array[Number]]} A The matrix
 * @return {String}   The output string: [[a,b,c,..],[d,e,f,..],...]
 */
var readableMatrix = function (A) {
  return "[" + A.map(s => `[${s.join(',')}]`).join('\n') + "]"
}

/**
 * Creates a LaTeX representation of a matrix.
 * @param  {Array[Array[Number]]} A The matrix
 * @return {String}   The output string: \begin{matrix} a & b & c & ... \\ d & e & f & .. \\ ...\end{matrix}
 */
var texMatrix = function (A) {
  return "\\begin\{pmatrix\}" + A.map(s => s.join('\&')).join('\\\\') + '\\end\{pmatrix\}';
}

/**
 * Converts a string to a matrix.
 * @param  {String} str The matrix in the form [[a,b,c,..]\n[d,e,f..]\n...]
 * @return {Array[Array[Number]]}     The output matrix as an array
 */
var readMatrix = function (str) {
  var rows = str.split('\n').map(s => s.replace(/(\[+|\]+)/,''));
  return rows.map(s => s.split(',').map(elem => parseFloat(elem)));
}

/**
 * Converts an array into a human-readable vector
 * @param  {Array[Number]} vec The vector
 * @return {String}     The vector in the form a \n b \n c \n ...
 */
var readableVector = function (vec) {
  return vec.join('\n');
}

/**
 * Converts a vector into a LaTeX table of coefficients, indexed by X_i
 * @param  {Array[Number]} vec The vector
 * @return {String}     A HTML table with LaTeX elements inside.
 */
var tableTexVector = function (vec) {
  return '<table>' + vec.map((s,i) => `<tr><td>\\(X_${i+1}\\)</td></tr>\\(${s}\\)</td></tr>`).join('') + '</table>';
}

/**
 * Converts a vector in string form into an array.
 * @param  {String} str The vector, in the form a \n b \n c \n...
 * @return {Array[Number]}     The output array
 */
var readVector = function (str) {
  return str.split('\n').map(s => parseFloat(s));
}

/**
 * Rounds a number to a required precision.
 * One subtlety is that if the value is within 10^(-4) of either 0 or 1,
 * it will treat this as being 0 or 1 resp.
 * This is due to the prevalence of 0.99998 and 0.00002 etc in variance results.
 * @param  {Number} n    The number to round
 * @param  {Number} prec The required precision (number of decimal places)
 * @return {Number}      The rounded number
 */
var mRound = function (n,prec) {
  return ((Math.abs(1-n) <= 0.0001) ? 1 : ((Math.abs(n) <= 0.0001) ? 0: ((Math.abs(n+1) <= 0.0001)? -1 : parseFloat(n.toFixed(prec)))));
}

// Exports to be used by other packages.
module.exports = {sign,copyArray,inArray,printArray,readableMatrix,texMatrix,readMatrix,readableVector,tableTexVector,readVector, mRound};
