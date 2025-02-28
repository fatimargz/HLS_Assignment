# HLS\_Assignment
HW 1 for C2-the-P2 Course @ NIU

check.ipynb : I created the fit function in python jupyter notebook first and passed the data.txt values. This was done to cross check with fit.C simulation results.

## Optimized fit function
fit.C : source code for fit function with optimization techniques. I used HLS interface pragmas, pipeline, and array partitioning. Inner loops were unable to unroll. 
fit.h : header that declares fit.C 

main.C : source code that calls in fit function. Used for testbench by calling in data from main.h. 
main.h : header that includes data as arrays. 

__Results__ : 
simulation\_with\_HLSBitType.png : screenshot of parameter results from simulation. Print statements include results for parameters and their uncertainties along with chi2. These results do not match with correct/expected values found in check.ipynb. This is probably due to precision error with ap-byte types.

synthesis\_with\_HLSBitTypes.png : screenshot result after synthesizing fit function. Latency = 1. I do get warning about unable to uroll inner loops. This is because inner loops rely on information after outerloop.

## Unoptimized fit function
fit1.C : My semi-optimized function that used float and int types.

simulation\_preHLSBitTypes.png : This is a screenshot of the simulation results before I changed the types to bit types. The parameters, uncertainties, and chi2 results are exactly as expected in comparison to my algorithm in python, check.ipynb.
 
