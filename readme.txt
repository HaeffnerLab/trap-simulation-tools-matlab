1. Introduction
---------------

This package is used to read and post-process electrostatic simulations of ion trap potentials. 

Before you use this package you must use a numerical electrostatics solver which returns the potentials 
of all the electrodes of a trap. These potentials of all electrodes are stored in a text file with format:

X1 Y1 Z1 V1
...
XN YN ZN VN

You can find example simulations under eurotrap-pt1.txt etc. These were produced using Kilian Singer's 
BEMsolver package. 

2. Preliminaries and naming conventions
---------------------------------------

We have adapted this library to use with surface electrode traps. You should be bale to use it with any 
other trap geometry. Our convention for numbering the trap electrodes is: 
  a. Each trap has N electrodes and one ground electrode (the RF electrode is included to the N)
  b. For a linear segmented trap We start counting on the lower left corner. The lower left DC electrode 
     is 1. As you move up along the axis (staying on the same side of the RF rail you count up). When you 
     reach the top, you cross to the right side of the RF rail and start counting at the bottom. When you 
     reach the top, you cross to the inside of the RF rail and count DC electrodes from bottom to top. When 
     you are done counting all the DCs, you count the RF electrodes.

3. Instructions 
---------------

Day-to-day usage will only involve messing around with analyze_trap,  project_parameters, and set_voltages. 
You will occasionally need to dive into the lower level functions.
_______________________________________________________________________________________________________________________
-----------------------------------------------------------------------------------------------------------------------

Main scripts and functions
--------------------------

analyze_trap: 
    This is the main script which coordinates all the other scripts. 
    Do not remove any of the function or script calls from here, unless you understand what you are doing. 
    You are allowed to flip the boolean and string flags which determine which of the following functions 
    plot what. 
    
project_parameters: 
    Here you define all the parameters pertaining to the trap simulation, and the system parameters 
    (paths for reading/writing data). This is all the messing around you will need to do. Any attempts 
    to change parameters in the following functions will probably lead to regret and frustration.

set_voltages:
    This is a script, not a function, and you should us this to edit voltage/multipole parameters before running ppt2. 
    
import_data:
    You should not have to look inside this function, unless you are trying something new (and possibly inappropriate). 
    import_data reads the text files which your BEM solver simulations generated, and converts them to a .mat structure (matlab)
    or a pickled data object (python). import_data *does not* make any changes to the structure 'data'. This is the job 
    of the next function. 
    
    The main thing to remember is that consecutive BEM solver simulations must have overlapping 
    first and last axial positions, i.e. [Z(last)]_simulation(n) = [Z(first)]_simulation(n+1), where Z is along the 
    trap axis.
    
get_trapping_field:
    You should not have to look inside this function, unless you are trying something new (and possibly inappropriate). 
    get_trapping_field adds the grid and electrostatic potential arrays to the structure 'data'. It looks through the data 
    structures which BEM solver generated and creates new data structure whith a grid which is centered around 
    the position where you want to trap. You define the trapping position in project_parameters. 
    
    At the end of get_trapping_field, the structure 'data' should contain fields such as data.trapConfiguration.X, 
    data.trapConfiguration.EL_DC1, etc.

expand_field:
    You should not have to look inside this function, unless you are trying something new (and possibly inappropriate). 
    expand_field does the spherical harmonic expansion of the potentials for all the used electrodes. It then stores the 
    harmonic expansion coefficients in the array data.trapConfiguration.multipoleCoefficients.
    
    As an option, expand_field will also compute all the potentials from the spherical harmonic coefficients and replace the 
    original values with the computed ones. This can be useful as a data smoothing step so that the algorithms which use 
    numerical derivatives (e.g. in ppt2 for pseudopotential calculation, find_saddle, trap_depth) work better.
    
trap_knobs:
    You should not have to look inside this function, unless you are trying something new (and possibly inappropriate). 
    trap_knobs calculates the independent multipole control parameters (see G. Littich Master's thesis for an explanation of 
    what these are). It stores the control parameters in data.trapConfiguration.multipoleControl. It also stores the kernel 
    vectors of the multipoleControl space in data.trapConfiguration.multipoleControl. 
    
    Normally, you are done at this point, you can save the so-called data.trapConfiguration.multipoleControl array to a file 
    and import it to LabRad. As of 4/4/2014, the array is saved as Multipole_Control_File.txt, under the _post_processed folder
    
set_dc:
    You should not have to look inside this function, unless you are trying something new (and possibly inappropriate). 
    set_dc takes as inputs multipole parameters or a set of Mathieu parameters, and produces an 1-D array with the voltages
    you need to apply to all of the electrodes. The entries corresponding to multipole-controlled electrodes receive some 
    values. If some electrodes are under manual control (per your definition in project_parameters) these are set to zero. 
    You can add the values to these electrodes at a higher level (as in set_voltages)
    
post_process_trap:
    You should not have to look inside this function, unless you are trying something new (and possibly inappropriate). 
    ppt2 is a function which takes in a set of electrode voltages, RF frequency, and amplitude and calculates what the trap 
    should look like for the ionb. it plots the rf potential, rf pseudopotential, dc potentail, and total trap potential. It 
    also calculates the trap frequencies, trap depth, and trapping position. The options for plotting these potentials are;
    'no plots', '1d plots', '2d plots', and '2d plots and 1d plots'.
    
    It updates the data.trapInstance structure with all these parameters that it has calculated. 
    
    Two obsolete (pre-2009) functionalities of ppt2 are: calculating micromotion compensation parameters given a stray 
    electric field, and calculating a stray electric field, given micromotion compensated voltages. You will not use these, 
    unless you have an interest in the history of science and technology.
    

Lower level scripts and functions
---------------------------------
....

_______________________________________________________________________________________________________________________
-----------------------------------------------------------------------------------------------------------------------
To Do:

* Double check that data stiching is done correctly...OK

* Double check that the electrode mapping and  multipole masking is done properly

* Homogenize the naming conventions of scripts and functions. Smooth matlab/python naming inconsistencies

* Expand this readme to include lower level functions

* Add more plotting flags to the remaining functions

* Eliminate name-space polution after post_process_trap

* Generate a set of BEM solver test data which allow transparent tutorial/debugging

* Join import_data functionality with the now defunct script for importing data from CPO

* Change the length scales to micron, but keep mm for the E and U

* Internally rescale to geometric mean for spher_harm_exp