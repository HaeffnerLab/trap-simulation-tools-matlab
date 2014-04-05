function data = project_parameters
% The electrode counting scheme is: first count electrodes on the left of the RF, bottom to top, 
% then count electrodes on the right of the RF, bottom to top, then Center electrode, then RF

% SI values for physical constants etc 
qe=1.60217646e-19; 
mp=1.67262158e-27;

platform = 'matlab';
projectName = 'eurotrap-pt';                                               %%%% used in import_data
                                                                           % Usually this is the name of the trap. 
                                                                           % A directory '/simulationDirectory/projectName_time_post-processed' will be created. 
                                                                           % This is where the post processed data files are stored.
%simulationDirectory = ...
%     '/media/Data/quantumlab/trap-simulations/G_trap_translume_v2/'; 	   %%%% used in import_data
simulationDirectory = '/home/nikos/qlab/trap-simulations/eurotrap/';       %%%% used in import_data. 
                                                                           % This is where the simulation files are stored.
dataNames= 'eurotrap-pt';                                                  %%%% used in import_data	
                                                                           % Used to open BEM txt files, named: dataNames1.txt, dataNames2.txt,...
                                                                           % Also used as the name of the .mat structure where post-processed data are saved. 
timestarted = datestr(now,'mm-dd-yyyy-HHhMM');                             %%%% used in import_data
%file_name = sprintf('%s%c%s',projectName,'_',timestarted);                %%%% used in import_data
                                                                           % This is the _time_ in '/simulationDirectory/projectName_time_post-processed
matDataPath = ...
  [simulationDirectory,projectName,'_',timestarted,'_','post-processed/']; %%%% used in import_data
                                                                           % This is where the post processed data files are stored.
NUM_AXIS  = 21;                                                            %%%% used in import_data 
                                                                           % Number of data points per axis used in import_data
NUM_ELECTRODES  = 22;                                                      %%%% used in import_data, get_trapping_field, expand_field, dcpotential, post_process_trap
                                                                           % Number of non-ground electrodes, including all DCs, center, RF
NUM_USED_ELECTRODES = 22;                                                  %%%% used in import_data
                                                                           % If the RF is not biased, set this to NUM_ELECTRODES-1, etc.

electrodeMapping = ...                                                     %%%% used in import_data, expand_field
            [1 1; 2 2; 3 3; 4 4; 5 5; ...                                  % Here define the electrode combinations
		    6 6; 7 7; 8 8; 9 9; 10 10; ...                                 % the convention is physical electrode -> functional electrode
		    11 11; 12 12; 13 13; 14 14; 15 15; ...                         % if electrodes 1 and 2 are combined into one electrode, then enter [1 1; 2 1; 3 2;...
		    16 16; 17 17; 18 18; 19 19; 20 20; ...                         % if electrodes 1 and 4 are not in use (grounded), then enter [1 0; 2 1; 3 2; 4 0...
		    21 21; 22 22];                                                 % NUM_ELECTRODES (i.e. last) is the RF electrode
manualElectrodes = ...                                                     % used in import_data, expand_field, dcpotential
                   [0 0 0 0 0 ...                                          % NUM_ELECTRODES-1 (i.e. before the RF) is the center electrode
                    0 0 0 0 0 ...                                          % electrodeMapping determines the pairing
                    0 0 0 0 0 ...                                          % manualElectrodes determines the electrodes which are under manual voltage control. It has NUM_ELECTRODES elements
                    0 0 0 0 0 ...                                          % (i.e. they are not connected to an arbitrary voltage, not to multipole knobs, )
                    0 0]';                                                 % all entries != 0 are under manual control, and entries = 0 are not under manual control  
if sum((electrodeMapping(:,1)>0))~=NUM_ELECTRODES
  fprintf('project_parameters warning: electrodeMapping variable does not agree with number of electrodes.\n');
end
if sum((electrodeMapping(:,2)>0))~=NUM_USED_ELECTRODES
  fprintf('project_parameters warning: electrodeMapping variable does not agree with number of used electrodes.\n');
end
if sum(electrodeMapping(:,2).*manualElectrodes)>0
  fprintf('project_parameters warning: some electrodes are both under multipole and manual control.\n')
end
usedMultipoles = ...                                                       %%%% used in trap_knobs
               [1 1 1 1 1 1 1 1];                                          % Select which multipoles you want to control
                                                                           % 0 means do not want to control something, 1 means control it. 
                                                                           % The sequence of multipoles is 
                                                                           %	(1)x,(2)y,(3)z,						(dipoles)
                                                                           %	(4)x^2-y^2,(5)2z^2-x^2-y^2,(6)xy,(7)yz,(8)zx		(quadrupoles)
                                                                           % Enter only 0 and 1 here, otherwise trap_knobs will make mistakes
                                                                           
axesPermutation = [2 3 1];                                                 %%%% used in import_data 
                                                                           % Axes permutation for importing the data from BEM solver. This reflects how the 
                                                                           % autocad drawing axes were defined
                                                                           % If drawing uses x - axial, y - radial, z - height, use perm = [2 3 1] (Eurotrap)
                                                                           % If drawing uses y - axial, x - radial, z - height, use perm = [1 3 2] (Sqip D trap, GG trap)
                                                                           % If drawing uses Nikos convention, use perm = [1 2 3]
                                                                           % This code assumes: 
                                                                           % x -> radial, parallel to trap
                                                                           % y -> radial, height above trap
                                                                           % z -> axial, parallel to trap

nStart = 1;                                                                %%%% used in import_data
nMatTot = 6;                                                               %%%% used in import_data, get_trapping_field
                                                                           % The trap simulation data is stored in text files describing overlapping 
                                                                           % data grids. The text files have names dataNames1.txt, etc
                                                                           % import_data.m can import multiple .txt files (see for loop). The number of
                                                                           % files to be imported can be adjusted via nMatTot.
                                                                           % nStart: index of the trap simulation file on which you want to start importing data
                                                                           % nMatTot: number of simulation files
scale = 1/1000;                                                            %%%% used by import_data, get_trapping_field, post_process_trap, setdc, spherharmxp
                                                                           % UNITS: if drawing (autocad) dimensions are in microns and output in mm then
                                                                           % scale = 1/1000 if drawing is in mm, scale = 1;
zMin = -630*scale;                                                         %%%% used in import_data, get_trapping_field
zMax = -510*scale;                                                         %%%% used in import_data, get_trapping_field
zStep = 20*scale;                                                          %%%% used in import_data, get_trapping_field
                                                                           % zMin,zMax: Range covered by all simulations 
                                                                           % zStep: Range covered by each .MAT file
trappingPosition = -615*scale;                                             %%%% used in get_trapping_field, expand_field, post_process_trap
                                                                           % trapping in mm along the trap axis.
regenerate_data = true;                                                    %%%% used in expand_field
                                                                           % Set this to true ifyou want expand_field to calculate the DC electrode potentials 
                                                                           % from the harmonic expansion coefficients                                                                        
driveAmplitude = 100;                                                      %%%% used in post_process_trap
driveFrequency = 40e6;                                                     %%%% used in post_process_trap
                                                                           
data = {};                                                                 % data stores all information about the trap
data.systemInformation.platform = platform;                                % This is used everywhere else. The rest is just redefinitions
data.systemInformation.timestarted = timestarted;
data.systemInformation.simulationDirectory = simulationDirectory;
data.systemInformation.projectName = projectName;
data.systemInformation.dataNames = dataNames;
data.systemInformation.matDataPath = matDataPath;
data.systemInformation.dateStarted = datestr(now);
data.systemInformation.axesPermutation = axesPermutation;

data.trapConfiguration.charge = qe;
data.trapConfiguration.mass = 40*mp;
data.trapConfiguration.NUM_AXIS = NUM_AXIS;
data.trapConfiguration.NUM_ELECTRODES = NUM_ELECTRODES;                    % Number of electrodes in the trap, includes all DC, center electrode, RF. Do not include ground
data.trapConfiguration.NUM_USED_ELECTRODES = NUM_USED_ELECTRODES;          % Number of electrodes to which DC voltage is applied. If two electrodes are shorted together, they count as one. 
                                                                           % If an electrode is kept at DC ground, do not count it here.
data.trapConfiguration.electrodeMapping = electrodeMapping;
data.trapConfiguration.manualElectrodes = manualElectrodes;
data.trapConfiguration.usedMultipoles = usedMultipoles;

data.trapConfiguration.nStart = nStart;                                    % nStart: index of the trap simulation file on which you want to start importing data
data.trapConfiguration.nMatTot = nMatTot;                                  % nMatTot: number of simulation files
data.trapConfiguration.zMin = zMin;                                        % zMin,zMax: Range covered by all simulations
data.trapConfiguration.zMax = zMax;                                        % 
data.trapConfiguration.zStep = zStep;                                      % zStep: Range covered by each .MAT file

data.trapConfiguration.trappingPosition = trappingPosition;                % position where you want to trap
data.trapConfiguration.scale = scale;                                      % scaling factor for the simulation file. 1 means mm, 1e-3 means micron, etc. 
data.trapConfiguration.r0 = 1;                                             % length scale parameter used for the spherical harmonic expansions r0 =1 means millimeters
data.trapConfiguration.expansionOrder = 9;                                 % this is the order to wich a spherical harmonic expansion is carried out
data.trapConfiguration.regeneratedPotentials = regenerate_data;            % this is true if expand_field is recalculating the potential data from the spherical harmonic expansion coefficients
data.trapConfiguration.Xcorrection = 0;                                    % Add a 'Xcorrection' offset to the RF saddle point, in case the simulation data was incorrect.  Obsolete. See expand_field
data.trapConfiguration.Ycorrection = 0;                                    % Add a 'Ycorrection' offset to the RF saddle point, in case the simulation data was incorrect.
data.trapConfiguration.regularize = true;

%%%%%% USED BY ANALYZETRAP %%%%%%

data.trapInstance.driveAmplitude = driveAmplitude;                         % Applied RF amplitude for post_process_trap analysis
data.trapInstance.driveFrequency = driveFrequency;                         % RF frequency for post_process_trap analysis
data.trapInstance.Efield = [0 0 0];                                        % Electric field that I want the trap to generate at the ion E = [Ex Ey Ez]
%params.scale = 1;

