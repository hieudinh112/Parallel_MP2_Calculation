# Helper Files

This directory contains parsers for converting Q-Chem output to data usable by the code. There
are some helper files to plot the scaling.


## Libraries 

python/3.9.12 is required to run this code.

## Files

- c_parse: Converts coefficient matrix output from Q-Chem outputs to readable format.

- parse: Extracts number of electrons, number of basis function, number of auxiliary basis functions and molecular orbital energies from Q-Chem output.

- overlap_parse: Converts overlap matrix outputs from Q-Chem outputs to readable format.

- eri_parse: Converts three electron intergal matrix outputs from Q-Chem outputs to readable format.
