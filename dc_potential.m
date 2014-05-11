function Vout = dc_potential(trap,V,Ex,Ey,Ez,x,y,z)
%EC01function Vout = dc_potential(trap,VMULT,VMAN,Ex,Ey,Ez,x,y,z)
%EC01%function Vout = dc_potential(trap,VMULT,VMAN,Ex,Ey,Ez,x,y,z)
% calculate the dc potential if I give you the applied voltages and the stray field
% trap: 
%   structure containing the voltages for all electrodes
% VMULT: 
%   electrode voltages determined by the multipole control algorithm.
% VMAN: 
%   electrode voltages determined by manual user control 
%	e.g.VMAN  = [0 0 -2. 0 ...] applies -2 V to electrode 3
% IMAN: 
%   array marking by an entry of 1 the electrodes which are under manual control, 
%	e.g. IMAN = [0 0 1 0 ...] means that electrode 3 is under manual control
% 	*BOTH* above conditions are necessary to manually apply the -2 V to electrode 3
% EX,EY,EZ: 
%   stray electric field
% NUM_ELECTRODES: 
%   number of DC electrodes in use (defined in project_parameters.m)
% x,y,z: 
%   matrices describing the grid in x,y,z
% 
% Nikos, cleaned up June 2013

    NUM_ELECTRODES = trap.Configuration.NUM_ELECTRODES;
    Vout = zeros(size(trap.Simulation.EL_DC1));  
    for ii=1:NUM_ELECTRODES
        Vout = Vout + V(ii)*trap.Simulation.(['EL_DC' num2str(ii)]);
    end
    
    Vout = Vout-Ex*x-Ey*y-Ez*z;

end
