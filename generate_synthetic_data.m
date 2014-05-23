% Make synthetic data to test that the code works properly.
% Assume a 'trap' with NUM_ELECTRODES electrodes. Each electrode generates 
% a simple cobmination of multipoles. The potentials are saved to a text
% file, similar to what you would get out of a BEM solver simulation.

fprintf('generating synthetic data...');
NUM_ELECTRODES = 14;
simulationDirectory = '/home/nikos/qlab/trap-simulations/synthetic_trap/';      
DataNames= 'synth_trap';                                                   %%%% used in  import_data	

% electrode i generates the multipoles in multipoleCoeff(:,i)
% Be careful with the units here. The units are V/mm for electrid field and
% V/mm^2 for the quadrupole coefficients 
%                 E0...                     E13    remember that BEMsolver
%                                                  starts counting at RF                                                  
multipoleCoeff = [0 0 0 0 0 0 0 0 0 0 0 0 0 0; ... Potential offset
                  0 1 0 0 0 0 0 0 0 0 0 0 0 0; ... Ex
                  0 0 1 0 0 0 0 0 0 0 0 0 0 0; ... Ey
                  0 0 0 1 0 0 0 0 0 0 0 0 0 0; ... Ez
                  1 0 0 0 1 1 1 1 0 0 0 0 0 0; ... U1 = (x^2-y^2)/2
                  0 0 0 0 1 0 0 0 1 1 1 0 0 0; ... U2 = (2z^2-x^2-y^2)/2
                  0 0 0 0 0 1 0 0 1 0 0 1 1 0; ... U3 = xy
                  0 0 0 0 0 0 1 0 0 1 0 1 0 1; ... U4 = yz
                  0 0 0 0 0 0 0 1 0 0 1 0 1 1 ];% U5 = xz

scale = 1/1000;                                 % this is for conversion from mm to micron
NUM_AXIS = 21;                                  % points per axis
grid = [-10 -10 -10 1 1 1];                  % Xmin Ymin Zmin dX dY dZ
% in what follows, generate the synthetic potentials and write to a text
% file. Use for loops to avoid ambiguities with the vectorized
% implementations in Matlab. The point of this excercise is not speed, but
% consistency, e.g. are you using a right handed or left handed coordinate
% system?
synthFileName = [simulationDirectory, DataNames, sprintf('%d',1), '.txt'];
f_to = fopen(synthFileName,'w');
for el = 1:NUM_ELECTRODES
    for x = grid(1):grid(4):grid(1)+(NUM_AXIS-1)*grid(4)
        for z = grid(3):grid(6):grid(3)+(NUM_AXIS-1)*grid(6)
            for y = grid(2):grid(5):grid(2)+(NUM_AXIS-1)*grid(5)
                v = potential_from_multipoles(x,y,z,multipoleCoeff(:,el),scale);
                s = [sprintf('%3.1f ', x), sprintf('%3.1f ', z), sprintf('%3.1f ', y), sprintf('%10.8f\n', v)];
                fprintf(f_to,s);
            end
        end
    end
end
fclose(f_to);
fprintf('done!\n');