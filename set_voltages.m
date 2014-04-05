%%%%%%%%%%%%%%%%%%%%%%%%% The following are voltage controls in case I decide to run post_process_trap to analyze specific inputs 
multipoleControls = true;                                                  % See set_dc for more explanations
regularizedc = false;                                                      % 
az = 4.5e-3;                                                               % Mathieu alpha_z parameter (valid only if multipoleControls == false) 
ax = -0.002;                                                               % Mathieu alpha_x parameter (valid only if multipoleControls == false)
phi = .0;                                                                  % Angle of rotation of DC multipole wrt RF multipole (valid only if multipoleControls == false)
Ex = .0;                                                                   % Stray electric field at the ion E = [Ex Ey Ez] (valid if multipoleControls == true)
Ey = .0;
Ez = .0;
U1 = .0;                                                                   % U1-U5 are the DC Quadrupoles that I want the trap to generate at the ion                                                                            % (valid if multipoleControls == true)
U2 = 12.0; 
U3 = .0;  
U4 = .0; 
U5 = .0;

manualvoltages = false;                                                    % Set this to true if you want to control the electrode voltages copmpletely manually
if manualvoltages,
    el = 	[0.0000;	-0.2724;	1.1564;		2.7547;		-11.5354;...
		 2.7670;	1.0742;		-0.37280;	-0.3178;	0.0000;...
		 0.0000;	0.0000;		-0.49017;	-0.9589;	-4.9847;...
		-3.7588;	-4.7823;	-0.94308;	-0.5042;	-0.3218;...
		 0.0000;	0.0000;		-0.57479];
else
    el = set_dc(data,[-Ex,-Ey,-Ez]',[U1,U2,U3,U4,U5]',[ax,az,phi],multipoleControls,regularizedc);
end
data.trapInstance.dcVoltages = el;
