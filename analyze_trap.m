% analyze_trap
% Generic trap analysis template. Coordianates more specific tools.
% Define the project parameters.
% Open up the data of a simulated trap file. Get the voltages uses for 
% a given set of multipoles and call post_process_trap to do the rest.
% data is a structure with trap simulation information (potentials etc.) 
% params is a structure with all the trap parameters: voltages, RF setup, 
% ion position, stray pseudofield, frequencies, axes tilt, depth, quadrupole
% coefficients, alpha and q parameters

get_trapping_field_plot = true;
trap_knobs_plot = true;
data = project_parameters;
data = import_data(data);
data = get_trapping_field(data,get_trapping_field_plot);
data = expand_field(data);
data = trap_knobs(data,trap_knobs_plot);

set_voltages;
dcplot = '1d plots';
rfplot = '1d plots';
pseudpotplot = '1d plots';
trappotplot = '1d plots';
data = post_process_trap(data,'analyzeTrap',dcplot,rfplot,pseudpotplot,trappotplot);

mesg = sprintf('The secular frequencies are: (%G, %G, %G) Hz.\n', data.trapInstance.frequency(1), data.trapInstance.frequency(2), data.trapInstance.frequency(3));
disp(mesg);
mesg = sprintf('The trapdepth is: %G eV \n', data.trapInstance.trapDepth);
disp(mesg);
mesg = sprintf('The ion sits at (%G,%G,%G) micron.\n', 1e3*data.trapInstance.ionPosition(1), 1e3*data.trapInstance.ionPosition(2), 1e3*data.trapInstance.ionPosition(3));
disp(mesg);

CC=data.trapConfiguration.multipoleControl;
dlmwrite(strcat(data.systemInformation.matDataPath,'Multipole_Control_File.txt'), CC, 'delimiter', ' ')



