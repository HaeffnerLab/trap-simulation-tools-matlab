function trap = project_parameters
% The electrode counting scheme is: first count electrodes on the left of the RF, bottom to top, 
% then count electrodes on the right of the RF, bottom to top, then Center electrode, then RF

% SI values for physical constants etc 
qe=1.60217646e-19; 
mp=1.67262158e-27;

platform = 'matlab';
projectName = 'eurotrap-pt';                                               %%%% used in  import_data
                                                                           % Usually this is the name of the trap. 
                                                                           % A directory '/simulationDirectory/projectName_time_post-processed' will be created. 
                                                                           % This is where the post processed trap files are stored.
%simulationDirectory = ...
%     '/media/Data/quantumlab/trap-simulations/G_trap_translume_v2/'; 	   %%%% used in  import_data
simulationDirectory = '/home/nikos/qlab/trap-simulations/eurotrap/';       %%%% used in  import_data. 
                                                                           % This is where the simulation files are stored.
DataNames= 'eurotrap-pt';                                                  %%%% used in  import_data	
                                                                           % Used to open BEM txt files, named: DataNames1.txt, DataNames2.txt,...
                                                                           % Also used as the name of the .mat structure where post-processed trap are saved. 
timestarted = datestr(now,'mm-dd-yyyy-HHhMM');                             %%%% used in  import_data
%file_name = sprintf('%s%c%s',projectName,'_',timestarted);                %%%% used in  import_data
                                                                           % This is the _time_ in '/simulationDirectory/projectName_time_post-processed
matDataPath = ...
  [simulationDirectory,projectName,'_',timestarted,'_','post-processed/']; %%%% used in  import_data
                                                                           % This is where the post processed trap files are stored.
NUM_AXIS  = 21;                                                            %%%% used in  import_data 
                                                                           % Number of trap points per axis used in  import_data
NUM_ELECTRODES  = 22;                                                      %%%% used in  import_data, get_trapping_field, expand_field, dc_potential, post_process_trap
                                                                           % Number of non-ground electrodes, including all DCs, center, RF

electrodeMapping = ...                                                     %%%% used in  import_data, expand_field
            [1 1; 2 2; 3 3; 4 4; 5 5; ...                                  % electrodeMapping determines which voltage goes to which electrodes. It has (NUM_ELECTRODES,2) elements
		    6 6; 7 7; 8 8; 9 9; 10 10; ...                                  % Here define the electrode combinations (for example electrode pairing)
		    11 11; 12 12; 13 13; 14 14; 15 15; ...                              % the convention is physical electrode -> functional electrode
		    16 16; 17 17; 18 18; 19 19; 20 20; ...                            % The last entry (NUM_ELECTRODES) is the RF electrode
		    21 21; 22 22];                                                 % Entry NUM_ELECTRODES-1 (i.e. before the RF) is the center electrode
                                                                           % Warning: this will only handle pairing up of adjacent electrodes
                                                                           %
manualElectrodes = ...                                                     % used in  import_data, expand_field, dc_potential
                   [0 0 0 0 0 ...                                          % manualElectrodes determines the electrodes which are under manual voltage control. It has NUM_ELECTRODES elements
                    0 0 0 0 0 ...                                          % The electrodes which are under manual control are connected to an arbitrary voltage source, and are not controlled 
                    0 0 0 0 0 ...                                          % via multipole knobs. 
                    0 0 0 0 0 ...                                          % All entries != 0 are under manual control, and entries = 0 are under multipole control 
                    0 0]';                                                 % If an electrode is grounded, set it under manual control, and do not set a manual voltage for it in set_voltages 
                                                                           % (voltage defaults to zero)
                                                                           %
manualVoltages = ...                                                       % manualVoltages detrmines the voltages on the electrodes which are under manual control.
                 [0 0 0 0 0 ...                                            % It is a vector with NUM_ELECTRODES elements. See examples below. 
                  0 0 0 0 0 ...                                            % [Note: If you decide to merge it with the manualElectrodes field, you will need to make sure you 
                  0 0 0 0 0 ...                                            % are not confusing the code into thinking that grounded electrodes are under multipole control] 
                  0 0 0 0 0 ...                                            %
                  0 0]';                                                  %
                                                                           %
                                                                           % Examples: 
                                                                           % If all electrodes are separated and are under mulitpole control, then enter 
                                                                           % electrodeMapping = [1 1; 2 2; 3 3; 4 4]   |
                                                                           % manualElectrodes = [0 0 0 0]              |------------> EL_DC1, EL_DC2, EL_DC3, EL_DC4 in trap.Simulation
                                                                           % manualVoltages = [0 0 0 0]                |
                                                                           % If electrodes 1 and 2 are combined into one electrode and are under mulitpole control, then enter 
                                                                           % electrodeMapping = [1 1; 2 1; 3 2; 4 3]   |
                                                                           % manualElectrodes = [0 0 0 0]              |------------> EL_DC1, EL_DC2, EL_DC3, EL_DC4 in trap.Simulation
                                                                           % manualVoltages = [0 0 0 0]                |
                                                                           % If electrodes 1 and 2 are combined into one electrode and are under manual control, then enter 
                                                                           % electrodeMapping = [1 1; 2 1; 3 2; 4 3]   |
                                                                           % manualElectrodes = [1 1 0 0]              |------------> EL_DC1, EL_DC2, EL_DC3, EL_DC4 in trap.Simulation
                                                                           % manualVoltages = [0.3 0.3 0 0 ]           |                **It is your responsibility to set the manual voltages
                                                                           %                                           |                to all paired up electrodes**
                                                                           % If electrodes 1 and 4 are not in use (grounded), electrode 2 is under manual control, then enter 
                                                                           % electrodeMapping = [1 1; 2 2; 3 3; 4 4]   |
                                                                           % manualElectrodes = [1 1 0 1]              |------------> EL_DC1, EL_DC2, EL_DC3, EL_DC4 in trap.Simulation                                                                                                                                                  
                                                                           % manualVoltages = [0 0.5 0 0 ]             |
                                                                           % If electrodes are paired up into twos and are under mulitpole control, then enter 
                                                                           % electrodeMapping = [1 1; 2 1; 3 2; 4 2]   |
                                                                           % manualElectrodes = [0 0 0 0]              |------------> EL_DC1, EL_DC2, EL_DC3, EL_DC4 in trap.Simulation
                                                                           % manualVoltages = [0 0 0 0]                |
                                                                           
if (sum((electrodeMapping(:,1)>0))~=NUM_ELECTRODES)||(sum((electrodeMapping(:,2)>0))~=NUM_ELECTRODES)
  fprintf('project_parameters warning: electrodeMapping variable does not agree with number of electrodes.\n');
end

usedMultipoles = ...                                                       %%%% used in trap_knobs
               [1 1 1 1 1 1 1 1];                                          % Select which multipoles you want to control
                                                                           % 0 means do not want to control something, 1 means control it. 
                                                                           % The sequence of multipoles is 
                                                                           %	(1)x,(2)y,(3)z,						(dipoles)
                                                                           %	(4)x^2-y^2,(5)2z^2-x^2-y^2,(6)xy,(7)yz,(8)zx		(quadrupoles)
                                                                           % Enter only 0 and 1 here, otherwise trap_knobs will make mistakes
                                                                           
axesPermutation = [2 3 1];                                                 %%%% used in  import_data 
                                                                           % Axes permutation for importing the trap from BEM solver. This reflects how the 
                                                                           % autocad drawing axes were defined
                                                                           % If drawing uses x - axial, y - radial, z - height, use perm = [2 3 1] (Eurotrap)
                                                                           % If drawing uses y - axial, x - radial, z - height, use perm = [1 3 2] (Sqip D trap, GG trap)
                                                                           % If drawing uses Nikos convention, use perm = [1 2 3]
                                                                           % This code assumes: 
                                                                           % x -> radial, parallel to trap
                                                                           % y -> radial, height above trap
                                                                           % z -> axial, parallel to trap

nStart = 1;                                                                %%%% used in  import_data
nMatTot = 6;                                                               %%%% used in  import_data, get_trapping_field
                                                                           % The trap simulation trap is stored in text files describing overlapping 
                                                                           % trap grids. The text files have names DataNames1.txt, etc
                                                                           %  import_data.m can import multiple .txt files (see for loop). The number of
                                                                           % files to be imported can be adjusted via nMatTot.
                                                                           % nStart: index of the trap simulation file on which you want to start importing trap
                                                                           % nMatTot: number of simulation files
scale = 1/1000;                                                            %%%% used by  import_data, get_trapping_field, post_process_trap, multipole_set_dc, spher_harm_exp
                                                                           % UNITS: if drawing (autocad) dimensions are in microns and output in mm then
                                                                           % scale = 1/1000 if drawing is in mm, scale = 1;
zMin = -630*scale;                                                         %%%% used in  import_data, get_trapping_field
zMax = -510*scale;                                                         %%%% used in  import_data, get_trapping_field
zStep = 20*scale;                                                          %%%% used in  import_data, get_trapping_field
                                                                           % zMin,zMax: Range covered by all simulations 
                                                                           % zStep: Range covered by each .MAT file
trappingPosition = -615*scale;                                             %%%% used in get_trapping_field, expand_field, post_process_trap
                                                                           % trapping in mm along the trap axis.
regeneratePotentials = true;                                                    %%%% used in expand_field
                                                                           % Set this to true ifyou want expand_field to calculate the DC electrode potentials 
                                                                           % from the harmonic expansion coefficients                                                                        
driveAmplitude = 100;                                                      %%%% used in post_process_trap
driveFrequency = 40e6;                                                     %%%% used in post_process_trap
                                                                           % These two parameters (RF amplitude and drive frequency) may be overwritten in set_voltages
                                                                           % Always check!
                                                                                                                                                   
trap = {};                                                                 % trap stores all information about the trap
trap.systemInformation.platform = platform;                                % This is used everywhere else. The rest is just redefinitions
trap.systemInformation.timestarted = timestarted;
trap.systemInformation.simulationDirectory = simulationDirectory;
trap.systemInformation.projectName = projectName;
trap.systemInformation.DataNames = DataNames;
trap.systemInformation.matDataPath = matDataPath;
trap.systemInformation.dateStarted = datestr(now);
trap.systemInformation.axesPermutation = axesPermutation;

trap.Configuration.charge = qe;
trap.Configuration.mass = 40*mp;
trap.Configuration.NUM_AXIS = NUM_AXIS;
trap.Configuration.NUM_ELECTRODES = NUM_ELECTRODES;                    % Number of electrodes in the trap, includes all DC, center electrode, RF. Do not include ground
trap.Configuration.electrodeMapping = electrodeMapping;
trap.Configuration.manualElectrodes = manualElectrodes;
trap.Configuration.manualVoltages = manualVoltages;
trap.Configuration.usedMultipoles = usedMultipoles;

trap.Configuration.nStart = nStart;                                    % nStart: index of the trap simulation file on which you want to start importing trap
trap.Configuration.nMatTot = nMatTot;                                  % nMatTot: number of simulation files
trap.Configuration.zMin = zMin;                                        % zMin,zMax: Range covered by all simulations
trap.Configuration.zMax = zMax;                                        % 
trap.Configuration.zStep = zStep;                                      % zStep: Range covered by each .MAT file

trap.Configuration.trappingPosition = trappingPosition;                % position where you want to trap
trap.Configuration.scale = scale;                                      % scaling factor for the simulation file. 1 means mm, 1e-3 means micron, etc. 
trap.Configuration.r0 = 1;                                             % length scale parameter used for the spherical harmonic expansions r0 =1 means millimeters
trap.Configuration.expansionOrder = 9;                                 % this is the order to wich a spherical harmonic expansion is carried out
trap.Configuration.regeneratePotentials = regeneratePotentials;            % this is true if expand_field is recalculating the potential trap from the spherical harmonic expansion coefficients
trap.Configuration.Xcorrection = 0;                                    % Add a 'Xcorrection' offset to the RF saddle point, in case the simulation trap was incorrect.  Obsolete. See expand_field
trap.Configuration.Ycorrection = 0;                                    % Add a 'Ycorrection' offset to the RF saddle point, in case the simulation trap was incorrect.
trap.Configuration.regularize = true;

%%%%%% USED BY ANALYZETRAP %%%%%%

trap.Instance.driveAmplitude = driveAmplitude;                         % Applied RF amplitude for post_process_trap analysis
trap.Instance.driveFrequency = driveFrequency;                         % RF frequency for post_process_trap analysis
trap.Instance.Efield = [0 0 0];                                        % Electric field that I want the trap to generate at the ion E = [Ex Ey Ez]
%params.scale = 1;

