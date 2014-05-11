% analyze_trap
% Generic trap analysis template. Coordianates more specific tools.
% Define the project parameters.
% Open up the data of a simulated trap file. Get the voltages uses for 
% a given set of multipoles and call post_process_trap to do the rest.
% trap is a structure with trap simulation information (potentials etc.) 
% voltages, RF setup,  ion position, stray pseudofield, frequencies, 
% axes tilt, depth, quadrupole coefficients, alpha and q parameters

get_trapping_field_plot = false;
trap_knobs_plot = false;
trap = project_parameters; 
trap = import_data(trap);
trap = get_trapping_field(trap,get_trapping_field_plot);
trap = expand_field(trap);
trap = trap_knobs(trap,trap_knobs_plot);

trap = set_voltages(trap);
dcplot = '1d plots';
rfplot = '1d plots';
pseudpotplot = '1d plots';
trappotplot = '1d plots';
trap = post_process_trap(trap,'analyzeTrap',dcplot,rfplot,pseudpotplot,trappotplot);

mesg = sprintf('The secular frequencies are: (%G, %G, %G) Hz.\n', trap.Instance.frequency(1), trap.Instance.frequency(2), trap.Instance.frequency(3));
disp(mesg);
mesg = sprintf('The trap depth is: %G eV \n', trap.Instance.trapDepth);
disp(mesg);
mesg = sprintf('The ion sits at (%G,%G,%G) micron.\n', 1e3*trap.Instance.ionPosition(1), 1e3*trap.Instance.ionPosition(2), 1e3*trap.Instance.ionPosition(3));
disp(mesg);

CC=trap.Configuration.multipoleControl;
dlmwrite(strcat(trap.systemInformation.matDataPath,'Multipole_Control_File.txt'), CC, 'delimiter', ' ')



