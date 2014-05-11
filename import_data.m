function trap =  import_data(trap_in)
% function  import_data(trap_in)
%
% Imports bemsolver .txt trap to matlab .mat trap, creates the directory 
% trapout.systemInformation.matDataPath, and saves the mat files in there.
% Takes .txt file as input and converts it into a trap structure .dat. Saves
% mat-structure under the same filename as the .txt. The potentials for the 
% trap electrodes, the grid vectors, and the grid parameters are stored
% as fields in the structure "Simulation"  
%
%  import_data.m can import multiple .txt files (see for loop). The number of
% files to be imported can be adjusted via nMatTot.
%
% nStart: index of the trap simulation file on which you want to start 
% importing trap
% nMatTot: number of simulation files
% 
% All the conventions concerning which electrodes are being used and which 
% ones are bound together were defined in project_parameters but they are 
% not implemented here. Here all electrodes are converted from txt to mat
%
% Rules:
% *All the electrodes are assumed to be DC electrodes to begin with.
% *The sequence for counting DC electrodes runs through left side of the RF 
% (bottom to top), right side of the RF (bottom to top), center electrodes 
% inside of the RF (left center, then right center), and finalLy RF. 
% *To specify that an electrode is grounded, go to project_parameters and
% set the corresponding parameter in electrodeMapping to 0 
% 
% Base code by Mike, modified by Gebhard Oct 2010. The conventions where 
% redefined/cleaning up/combined with later developments by Nikos Jun 2013
%

% unwrap variables
platform = trap_in.systemInformation.platform;
timestarted = trap_in.systemInformation.timestarted;
simulationDirectory = trap_in.systemInformation.simulationDirectory;
projectName = trap_in.systemInformation.projectName;
DataNames = trap_in.systemInformation.DataNames;
matDataPath = trap_in.systemInformation.matDataPath;
perm = trap_in.systemInformation.axesPermutation;
scale = trap_in.Configuration.scale;
nStart = trap_in.Configuration.nStart;
nMatTot = trap_in.Configuration.nMatTot;
NUM_AXIS = trap_in.Configuration.NUM_AXIS;
NUM_ELECTRODES = trap_in.Configuration.NUM_ELECTRODES;

trap = trap_in;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%check if older files exist in aux.txt nMatTot and perm must agree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if the mattrappath directory exists, and if it does not, create it
% if it does exist, just load trap from it
simfiles = dir(trap_in.systemInformation.simulationDirectory);
newest = 1e6; % assume some ridiculously old simulation directory
newestDirIndex = 0;
for i = 1:numel(simfiles)
    L = length(simfiles(i).name);
    strending = simfiles(i).name(max(L-14,1):L);
    if strcmp(strending,'_post-processed')
        if newest>days_distance(simfiles(i).date,'dd-mmm-yyyy HH:MM:SS')
            newest = days_distance(simfiles(i).date,'dd-mmm-yyyy HH:MM:SS');
            newestDirIndex = i;
        end
    end 
end

print_underlined_message('start',' import_data');
if newestDirIndex    
    trap.systemInformation.matDataPath = [simulationDirectory,simfiles(newestDirIndex).name, '/'];%%%out
    fprintf('Recent simulation directory already exists in simulation directory:\n');
    fprintf(sprintf('%s\nSkipping txt->mat conversion.\n',matDataPath));
    print_underlined_message(' stop',' import_data');
    return
end

if strcmp(platform,'octave') 
  mkdir (simulationDirectory,[projectName,'_',timestarted,'_','post-processed/']);
else
  mkdir (sprintf('%s%s',simulationDirectory,[projectName,'_',timestarted,'_','post-processed/']));
end

for iterationNumber=nStart:nMatTot

    %%%%%
    %			PART I: read txt file
    %%%%%
    fprintf(['Importing ',sprintf('%s%i',DataNames,iterationNumber),'.txt \n']);
    
    DataFromTxt = load([sprintf('%s',simulationDirectory),sprintf('%s',DataNames), sprintf('%d',iterationNumber), '.txt']);  % this instance of load returns a matrix with the numbers in the specified text file

    %find the x y and z grids
    y = 0;
    x = 0;
    z = DataFromTxt(1:NUM_AXIS, 3)*scale;
    for i=0:(NUM_AXIS-1)
        y(i+1) = DataFromTxt(NUM_AXIS*i+1, 2)*scale;
        x(i+1) = DataFromTxt(NUM_AXIS^2*i +1, 1)*scale;
    end
    x = x';
    y = y';
    coord = [x y z];
    x = coord(:,perm(1));
    y = coord(:,perm(2));
    z = coord(:,perm(3));
    
    %find  min, max and spacing of the axes
    xmin = min(x);
    xmax = max(x);
    ymin = min(y);
    ymax = max(y);
    zmin = min(z);   
    zmax = max(z);
    deltax = (xmax - xmin)/(NUM_AXIS-1);
    deltay = (ymax - ymin)/(NUM_AXIS-1);
    deltaz = (zmax - zmin)/(NUM_AXIS-1);

    %loads all the voltages and E vector into struct using dynamic naming
    for el=0:(NUM_ELECTRODES-1)
         for i=0:(NUM_AXIS-1)
            for j = 0:(NUM_AXIS-1)
                lb = NUM_AXIS^3*el + NUM_AXIS^2*i + NUM_AXIS *j + 1; %lower bound
                ub = NUM_AXIS^3*el + NUM_AXIS^2*i + NUM_AXIS *j +NUM_AXIS; %upper bound .
                struct.(['EL_phi' num2str(el)])(i+1,j+1,1:NUM_AXIS)=DataFromTxt(lb:ub, 4);
                
                if (size(DataFromTxt,2)>4)
                    % i.e. Ex,Ey,Ez are calculated in bemsolver (old
                    % version), fast
                    struct.(['EL_Ex' num2str(el)])(i+1,j+1,1:NUM_AXIS)=DataFromTxt(lb:ub, 5);
                    struct.(['EL_Ey' num2str(el)])(i+1,j+1,1:NUM_AXIS)=DataFromTxt(lb:ub, 6);
                    struct.(['EL_Ez' num2str(el)])(i+1,j+1,1:NUM_AXIS)=DataFromTxt(lb:ub, 7);
                else
                    % i.e. Ex, Ey, Ez are NOT calculated in bemsolver (slow
                    % bemsolver, more exact). Erf will be calculated by the
                    % numerical gradient in post_process_trap.m
                    struct.(['EL_Ex' num2str(el)])(i+1,j+1,1:NUM_AXIS)=0;
                    struct.(['EL_Ey' num2str(el)])(i+1,j+1,1:NUM_AXIS)=0;
                    struct.(['EL_Ez' num2str(el)])(i+1,j+1,1:NUM_AXIS)=0;
                end
                
            end
        end
        struct.(['EL_phi' num2str(el)]) = permute(struct.(['EL_phi' num2str(el)]), perm);
        struct.(['EL_Ex' num2str(el)]) = permute(struct.(['EL_Ex' num2str(el)]), perm);
        struct.(['EL_Ey' num2str(el)]) = permute(struct.(['EL_Ey' num2str(el)]), perm);
        struct.(['EL_Ez' num2str(el)]) = permute(struct.(['EL_Ez' num2str(el)]), perm);

    end

    %%%%%
    %			PART II: organize the electrodes in Simulation 
    %%%%%
    clear DataFromTxt
    if size(x,2)>size(x,1)
        Simulation.X = x;
        Simulation.Y = y;
        Simulation.Z = z;
    else
        Simulation.X = x';
        Simulation.Y = y';
        Simulation.Z = z';
    end
       
    Simulation.grid = [xmin ymin zmin deltax deltay deltaz];
    Simulation.EL_RF = struct.EL_phi0;	% Need to edit this if you are using out of phase! layer1=DC1 etc if layer 2 is DC1 if +1
    % initialize NUM_USED_DC electrodes
    for iii=1:NUM_ELECTRODES
        Simulation.(['EL_DC' num2str(iii)]) = zeros(NUM_AXIS,NUM_AXIS,NUM_AXIS);
    end
    % add each electrode.
    for iii=1:(NUM_ELECTRODES)
        Simulation.(['EL_DC' num2str(iii)]) = ...                           
            Simulation.(['EL_DC' num2str(iii)]) + ...                       
            struct.(['EL_phi' num2str(rem(iii,NUM_ELECTRODES))]);           
    end
     fprintf('Saving trap to matDataPath:\n%s\n',[sprintf('%s',matDataPath), sprintf('%s', DataNames), '_',sprintf('%d',iterationNumber), '.mat']);
    if strcmp(platform,'octave')
        save ('-v7',[sprintf('%s',matDataPath), sprintf('%s', DataNames), '_', sprintf('%d',iterationNumber), '.mat'],'Simulation');
    else
        save([sprintf('%s',trap.systemInformation.matDataPath), sprintf('%s', DataNames), '_', sprintf('%d',iterationNumber), '.mat'],'Simulation');
    end
end
% plot the RF potential
Ef = Simulation.EL_RF;
for a = 1:NUM_AXIS
    for b = 1:NUM_AXIS
        E(a,b) = Ef(a,b,NUM_AXIS);
    end
end
figure; 
imagesc(x,y,E'); 
title('import_ trap: Plotting the RF potential');
pause;
print_underlined_message(' stop',' import_data');


