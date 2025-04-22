# Data-driven-thermal-conduction
A data-driven solver framework for thermal conduction applied to in-situ data. The temperature measurements (in-situ data) is given in a separate file.

## Important notes:
Currently, the simulation framework is only tested under Ubuntu 20.04. The directories are named after the sections of a scientific paper for which the source code is used. The files load.dat contain the load function for the surface. The values originate from the top sensor from the in-situ data. The values in between are linearly interpolated.

## Required Packages
For running these source codes matplotlib and numpy are required. The tested versions are: numpy 1.24.4 and matplotlib 3.7.5.

## Theory
The source is utilized in a scientific publication, which is currently under review.

## Description
The source code provided here is for one-dimensional heat conduction problems. This challenges is solved via the finite element method and the data driven solver. Everything is implemented with Python 3.8.6.

## Temperature-measurements-BASt-interpolated
This directory contains a spread sheet with the sensor data measured by BASt in the duraBASt project. For each sensor the temperature values after outlier removal and interpolation are listed.

## Identified-material-parameters
In this directory, a spread is saved in which the best material parameter set per week is saved. These are the results of an optimization framework, which cannot be included in this repository. The developments are further discussed in the scientific publication. 

## Section 3.2
Here, the weeks 91-100 are used for material parameter identification with the finite element method, while week 90 is used for validation. All simulations can be started by executing the shell script runpy.sh with a `1` as additional argument.

## Section 4.2
In these directories, the data driven solver implementation is tested on artificially generated data. The data is generated directly in the python script dd.py. The script can be run with executing the shell script runpy.sh with `1` as additional argument. The weeks 91-100 are used to adjust the simulation paramters. Week 90 is used for validation purposes.

## Section 5.1
In these directories, the data driven solver is applied to in-situ data. In `data-generation`, the script for generating the data-sets is located. The file sensor.dat contains all sensor data for weeks 91-100. This is an input file for the python script lumped-C.py, which generates the data-sets stored in the files g.dat and q.dat. Those files are required for running the dd.py script in the other directories. Weeks 91-100 are used to adjust the simulation parameters, while week 90 is used for validation. The simulations are run via the shell script run.py.  In the file g.dat the temperature gradients are stored, while q.dat contains the corresponding heat fluxes.

## Section 5.2 
In these directories an optimized measruing concept is tested. In data_generation, all scripts required to generate the data sets are located. Firstly, a finite element simulation needs to be run by the script fe.py. Afterwards, the file generategandq.py needs to be executed to process the raw simulation data contained in T.dat and Q.dat into g.dat and q.dat. T.dat contains temperatures, while Q.dat contains the heat fluxes from the simulation. Weeks 91-100 are used to adjust simulation parameters. Week 90 serves as validation. The simulations are run via the shell script run.py. The directory Figure-19 contains the simulation script with additional output. Here, the final distance for each material point is saved for further processing. The file plot-distances.py generates a three-dimensional plot of the final distances over the one-dimensional specimen and time.

## Section 5.3 
In these directories, the data driven solver is shown to handle changing material characteristics. In data_generation, the script do.sh runs for the first, second and third week a finite element simulation and, afterwards, runs the script generategandq.py to fascilitate the simulation results into input data sets for the data driven solver. The first week is week 2 from the measurements, the second week is week from the measurements and the third week is week 4 from the measurements. In weeks2-4 the data driven simulation is run by the shell script runpy.sh. In it the data sets change by simulated time.
