Full publication: https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-3-117


Installation and execution of AntEpiSeeker Software



1) Installation:

For linux system version, users may simply unzip the package into one folder and run the program by typing "./AntEpiSeeker"(gcc version 3.4.6 
or above is needed).

For windows system version, unzip all files into one folder. The GNU Scientific Library (GSL) files "libgsl.dll",
"libslcblas.dll" and "WinGsl.dll" in this package should be put in the same folder with AntEpiSeeker as well as your system folder 
(e.g., "c:\\windows\\system32"). Start MS-DOS by Start>Run>cmd and change your working directory to AntEpiSeeker. 
Run the program by typing "AntEpiSeeker.exe".

2) Input Format

The user should create a comma-delimited file (by default, "data0.txt") which contains the case-control genotype data 
as the input for the program. The first line of the input file contains the SNP IDs and the name for the response variable (last columumn). 
The following lines are the genotype data which should be coded by 0, 1 and 2 with each line corresponding to one individual. 
The last column should contain the disease status of each individual coded by 0 and 1. 
The following is a sample data file for 5 individuals (3 cases and 2 controls) each genotyped for 10 SNPs.

rs0,rs1,rs2,rs3,rs4,rs5,rs6,rs7,rs8,rs9,class
2,2,2,2,1,0,1,2,2,0,1
0,1,2,0,1,0,2,2,2,1,1
2,1,0,1,1,2,0,2,1,0,1
1,1,2,0,1,0,1,0,2,2,0
0,2,0,0,1,0,1,1,0,0,0

3) Parameter file

The parameters for running AntEpiSeeker are specified in the "parameters.txt" file. These parameters include iAntCount,
iItCountLarge, iItCountSmall, alpha, iTopModel, iTopLoci, rou, phe, largehapsize, smallhapsize, iEpiModel, pvalue, 
INPFILE, OUTFILE. 

The parameter "iEpiModel" specifies the number of SNPs in an epistatic interaction. The parameters "largehapsize",
"smallhapsize" must be greater than "iEpiModel". For two-locus interaction model, we suggest largehapsize=6, 
smallhapsize=3, iEpiMode=2; For three-locus interaction model, we suggest largehapsize=6, smallhapsize=4, iEpiModel=3.

The parameters "iItCountLarge", "iItCountSmall" should be chosen according to the number of SNPs genotyped in the data
(Denoted by L) . Typically, we suggest iItCountSmall=0.1*L and iItCountLarge=0.5*iItCountSmall.


4) Output
Two output files will be generated by AntEpiSeeker. "AntEpiSeeker.log" under the AntEpiSeeker directory records the 
intermediate results including the detected top ranked haplotypes and the loci with top pheromone levels. The result 
file specified by the user in the parameter file reports the detected epistatic interactions with significant user-defined p value threshold
(after the Bonferroni correction).

5) Support
Questions and comments should be directed to Yupeng Wang (wyp1125@gmail.com).

