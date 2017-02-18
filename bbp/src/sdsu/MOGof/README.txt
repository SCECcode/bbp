a. Input seismograms must be the same format (i.e. acceleration)
b. Input seismograms must be the same length (seconds)
c. Input seismograms must have the same units (i.e. cm/s)
d. Each seismogram must have the same number of headers preceding it.
e. The columns following the headers are in the format (time[s], x comp, y comp, z comp)
f. Filtering cut-offs are in periods (seconds)
g. Must have an 'out' directory in the folder where the algorithm is to be run.
h. BBanalyze.m file must be in the 'out' directory if the MatLab routine is to be used 
    or the user can comment out the call to the BBanalyze.m routine in the OUTPUT.m file 
    (line 331 if no NGA values are calculated) after the algorithm has been run.
i. Optimal speed in the algorithms execution: if the number of points for each seismogram
    is limited to 8000 points or less.
j. If a segmentation fault occurs shortly after exectution try:
    ulimit -s unlimited
    to increase the stack size and allow more memory for the program.
k. The out folder will contain all of the exported data from the algorithm (except for a 
    *.log file).
l. The OUTPUT.m Matlab routine will organize the output data into an excel file (*.xls) and
    a MatLab structure (*.mat).
m. Calculated values are exported for peak elastic displacement, IE ratios, SA16, spectral 
    duration (normalized by period; horizontal components), PGA (horizontal components), 
    PGV (horizontal components) and the filtered/resized seismograms.
n. For continuity, the first file listed in PARAM.dat should be the recorded data and the 
    second file listed should be the synthetic data.
o. To change the strain hardening ratio, alter the value in InElastic_V1.2.F line 42:
      'data alpha/0.02/' to the desired value.
p. To change the damping ratio, alter the value in InElastic_V1.2.F line 43:
      'data z/0.05/' to the desired value and
    alter the value in SpecResp_V1.2.F line 19:
      'beta = 0.05' to the desired value.