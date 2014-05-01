function Vout = dc_potential(data,VMULT,VMAN,Ex,Ey,Ez,x,y,z)
%function Vout = dc_potential(data,VMULT,VMAN,Ex,Ey,Ez,x,y,z)
% calculate the dc potential if I give you the applied voltages and the stray field
% data: 
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

    IMAN = data.trapConfiguration.manualElectrodes;
    NUM_ELECTRODES = data.trapConfiguration.NUM_ELECTRODES;
    Vout = zeros(size(data.Simulation.EL_DC1));  
    for ii = 1:NUM_ELECTRODES
        if IMAN(ii),
            Vout = Vout + VMAN(ii)*data.Simulation.(['mEL_DC' num2str(ii)]);
        end
    end  
    for ii=1:NUM_ELECTRODES
        Vout = Vout + VMULT(ii)*data.Simulation.(['EL_DC' num2str(ii)]);
    end
    
    Vout = Vout-Ex*x-Ey*y-Ez*z;

end
