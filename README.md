# Relaxation Assisted Seperation Toolbox

Relaxation assisted seperation (RAS) is a method for analyzing and resolving NMR data using relaxation in MATLAB.


## Requirements
- All testing has been done on MATLAB 2020a and newer
- The Parrallel Computing Toolbox in MATLAB
- Previously developed NNLS functions are included with this distribution

## Installation
- Download the 'functions' folder and add it to the MATLAB path
- Example scripts and data are provided in the other folders to demonstrate usage

## Features
- Some basic processing example scripts are provided with experimental datasets to demonstrate how to pre-process relaxation data into the appropriate matrix or 3D object including: Reading Bruker FIDs/SERs, CPMG echo coaddition, phasing, automatic phasing up to second order, denoising relaxation data with PCA
- 2D RAS functions are provided for processing T1 and T2-type datasets
- 3D RAS scripts are provided for processing T1-T2 datasets
- The kernel structure of either of the above can be simply modified to process other types of exponential decays (e.g., T1rho, dipolar T1, diffusion, etc.)

## Citing
If you use RAS please cite the following:
TBD
A.R. Altenhof, M.J. Jaroszewicz, L. Frydman, and R.W. Schurko. 3D Relaxation-Assisted NMR: Achieving site-specific resolution in ultra-wideline NMR patterns. 
https://doi.org/###

## Creators
- Michael Jaroszewicz
- Adam Altenhof

## Contact
Please contact rschurko@fsu.edu with any quesitons, feedback, or suggestions.

## Support
This software was supported (in part) by the National Science Foundation Chemical Measurement and Imaging Program, with partial co-funding from the Solid State and Materials Chemistry Program (NSF-2003854).

## License
MIT

[//]: # ()

   [dill]: <https://github.com/joemccann/dillinger>
