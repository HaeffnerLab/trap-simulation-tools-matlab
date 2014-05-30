function trap_out = set_voltages(trap)
% function trap_out = set_voltages(trap)
% create a vector of voltages to use for the post_process_trap analysis
% this script has to be run after trap_knobs (because set_dc uses the trap 
% knobs control field in trap.Configuration.multipoleControl) but 
% before post_process_trap 

trap_out = trap;
multipoleControls = true;                                                  % See set_dc for more explanations
regularizeDC = false;                                                      % 
az = 4.5e-3;                                                               % Mathieu alpha_z parameter (valid only if multipoleControls == false) 
ax = -0.002;                                                               % Mathieu alpha_x parameter (valid only if multipoleControls == false)
phi = .0;                                                                  % Angle of rotation of DC multipole wrt RF multipole (valid only if multipoleControls == false)
Ex = .0;                                                                   % Stray electric field at the ion E = [Ex Ey Ez] (valid if multipoleControls == true)
Ey = .0;
Ez = .0;
U1 = .0;                                                                   % U1-U5 are the DC Quadrupoles that I want the trap to generate at the ion                                                                            % (valid if multipoleControls == true)
U2 = 10.0; 
U3 = .0;  
U4 = .0; 
U5 = .0;

el = multipole_set_dc(trap,[-Ex,-Ey,-Ez]',[U1,U2,U3,U4,U5]',[ax,az,phi],multipoleControls,regularizeDC);

for ii = 1:trap.Configuration.NUM_ELECTRODES
    if trap.Configuration.manualElectrodes(ii)
        el(ii) = trap.Configuration.manualVoltages(ii);
    end
end

trap_out.Instance.driveAmplitude = 1000;                         % Applied RF amplitude for post_process_trap analysis
trap_out.Instance.driveFrequency = 20e6;                         % RF frequency for post_process_trap analysis
trap_out.Instance.dcVoltages = el;
trap_out.Instance.E_in = [Ex Ey Ez];
trap_out.Instance.U_DC_in = [U1 U2 U3 U4 U5]; 
