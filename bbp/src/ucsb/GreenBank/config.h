#!/bin/csh -f 
setenv PREPROCFLAGS 
setenv CPP "/lib/cpp -traditional " 
setenv OPT_LO "ifort -c -auto -tpp6 -mp1 -O0 " 
setenv OPT_MED "ifort -c -auto -tpp6 -mp1 -O2 " 
setenv OPT_HI "ifort -c -auto -tpp6 -mp1 -ip -O3 " 
setenv LOAD "ifort" 
setenv LOADLIB 
