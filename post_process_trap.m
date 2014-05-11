function  out = post_process_trap(trap,operation,dcplot,rfplot,pseudpotplot,trappotplot,varargin) 
% post processing tool
% Start with input values for all the dc voltages and RF parameters. These
% must not include compensation parameters (added seperately).
% Find stray field that would be compensated in the given configuration.
% The axes order in the potential trap is:
% I=radial horizontal (X); J=radial vertical (Y); K=axial (Z)
% ie I->X, J->Y, K->Z
%
% trap: 
%   a cpo-simulation trap-structure with the electrode potentials in the
%   trapping region
% operation: 
%   determines the task performed:
%       'findEfield' determines the stray electric field for given dc
%       voltages
%       'findCompensation' determines compensation parameters for given 
%       starting voltages and given stray field
%       'analyzeTrap' do not optimize, just analyze the trap, assuming
%       everything is ok
% dcplot: 
%   What kind of plots to make for DC potential. Takes value in 
%   ('no plots', '2d plots','1d plots','2d and 1d plots') 
% rfplot:
%   What kind of plots to make for RF potentialat maximum voltage. Takes 
%   value in ('no plots', '2d plots','1d plots','2d and 1d plots') 
% pseudpotplot:
%   What kind of plots to make for pseudopotential. Takes value in 
%   ('no plots', '2d plots','1d plots','2d and 1d plots') 
% trappotplot: 
%   What kind of plots to make for total trapping potential. Takes value in 
%   ('no plots', '2d plots','1d plots','2d and 1d plots') 
% varargin: 
%   accepts two optional arguments passed to plot_potential, 
%   varargin is used to save pdf figures of plots
%   outPath: path to save figures to. This is usually set to matDataPath
%   pathSuffix: suffix to add to outPath
%
% Nikos, January 2009
% Cleaned up, March 2014

print_underlined_message('start','post_process_trap');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization
out = trap;
qe=trap.Configuration.charge;                                              % elementary charge in SI
mass = trap.Configuration.mass;                                            % mass of the ion
qualityCheck = true;                                                       % perform quality checks for the multipole expansions 
Zval = trap.Configuration.trappingPosition;
dcVoltages = trap.Instance.dcVoltages;
scale = trap.Configuration.scale;
grid = trap.Simulation.grid;
X = normalize(trap.Simulation.X); 
Y = normalize(trap.Simulation.Y); 
Z = normalize(trap.Simulation.Z);
[y x z] = meshgrid(Y,X,Z);
RFampl = trap.Instance.driveAmplitude;                                     % RF parameters
Freq = trap.Instance.driveFrequency;                    
Omega = 2*pi*Freq;   
r0 = trap.Configuration.r0;                                                % lengthscale of multipole expansion in millimeters
V0 = mass*(2*pi*Omega)^2*(r0*1e-3)^2/qe;                                   % god given voltage in SI 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end Intitalization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check RF potential
[Irf Jrf Krf] = find_saddle(trap.Simulation.EL_RF,X,Y,Z,2,Zval,'RF potential in post_process_trap', true);
fprintf('RF saddle indices: %d,%d,%d\n',Irf,Jrf,Krf);
pause;
warn('RF',trap.Simulation.EL_RF,Irf,Jrf,Krf);
Vrf = RFampl*trap.Simulation.EL_RF;
plot_potential(Vrf,Irf,Jrf,Krf,grid,rfplot,'RF potential','V_{rf} (Volt)',varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end check RF potential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check DC potential 
if ~isempty(trap.Instance.Efield),                                         % check if the initial guess for E is ok
    EE = trap.Instance.Efield;
    Udc = dc_potential(trap,dcVoltages,EE(1),EE(2),EE(3),x,y,z);
    % DC parameters                                                            
    [Idum Jdum Kdum] =  find_saddle(Udc,X,Y,Z,3,Zval);
    %plot_potential(Udc,Idum,Jdum,Kdum,grid,dcplot,'DC potential (stray field included)','U_{dc} (Volt)');
    if (warn('DC',Udc,Idum,Jdum,Kdum))&&(~strcmp(operation,'analyzeTrap')), 
        trap.Instance.Efield = []; 
        isempty(trap.Instance.Efield)
    end
end
% restore Ex,Ey,Ez to zero for d_e to run properly
Ex = 0; Ey = 0; Ez = 0;
Udc = dc_potential(trap,dcVoltages,Ex,Ey,Ez,x,y,z);

if strcmp(operation,'findEfield'),                                         % this option means find stray field                                                      
    fprintf('This option is no longer supported. Ask Nikos for an older version');
elseif strcmp(operation,'findCompensation'),                               % this option means find compensation voltages
    fprintf('This option is no longer supported. Ask Nikos for an older version');
elseif strcmp(operation,'analyzeTrap'),                                    % this option means do not optimize anything, and just analyze the trap
    fprintf('Running post_process_trap in plain analysis mode (no optimizations).\n');
    E = trap.Instance.Efield;
    dist = d_e(E);
    fprintf('Stray field is ( %G, %G, %G) V/m.\n',1e3*E(1),1e3*E(2),1e3*E(3));
    fprintf('With this field the compensation is optimized to %G micron.\n\n',1e3*dist);
else
    fprintf('\nInvalid ''operation'' input option. Quiting.\n');
    print_underlined_message(' stop','post_process_trap');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end check DC potential
% By this point you should have electrode voltages and dc potential values
% that you can directly plug into the auxiliary functions in order to
% analyze the trap instance. Auxiliary functions are exact_saddle, pfit, and
% exact_saddle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analyze trap instance 
Udc = dc_potential(trap,dcVoltages,E(1),E(2),E(3),x,y,z);
[XRF YRF ZRF] = exact_saddle(trap.Simulation.EL_RF,X,Y,Z,2,Zval);          % find secular frequencies etc.
[XDC YDC ZDC] = exact_saddle(Udc,X,Y,Z,3,Zval);
fprintf('RF saddle: (%f %f %f)\nDC saddle (%f %f %f).\n',XRF,YRF,ZRF,XDC,YDC,ZDC);

plot_potential(Udc,Irf,Jrf,Krf,grid,dcplot,'Compensated DC potential','U_{dc} (V)',varargin);
[IDC JDC KDC] = find_saddle(Udc,X,Y,Z,3,Zval);

[fx fy fz theta Depth rx ry rz xe ye ze U] = pfit(E(1),E(2),E(3));

Qrf = spher_harm_exp(Vrf,XRF,YRF,ZRF,7,X,Y,Z);                             % find RF multipole coefficients

if sqrt((XRF-XDC)^2+(YRF-YDC)^2+(ZRF-ZDC)^2)>0.008, 
    Qdc = spher_harm_exp(Udc,XRF,YRF,ZRF,7,X,Y,Z); 
else
    Qdc = spher_harm_exp(Udc,XDC,YDC,ZDC,7,X,Y,Z);
end    

Arf = 2*sqrt( (3*Qrf(8))^2+(3*Qrf(9))^2 );
thetaRF = 45*(sign(Qrf(9)))-90*atan((3*Qrf(8))/(3*Qrf(9)))/pi;
Adc = 2*sqrt( (3*Qdc(8))^2+(3*Qdc(9))^2 );
thetaDC = 45*(sign(Qdc(9)))-90*atan((3*Qdc(8))/(3*Qdc(9)))/pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end analyze trap instance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Return values
out.Instance.trapPotential = U;
out.Instance.E_out = E;
out.Instance.misCompensation = dist;
out.Instance.ionPosition = [XRF YRF ZDC];
out.Instance.ionPosIndex = [Irf Jrf Krf];
out.Instance.frequency = [fx fy fz];
out.Instance.theta = theta;
out.Instance.trapDepth = Depth/qe;
out.Instance.escapePosition = [xe ye ze];
out.Instance.U_RF_out = 2*[Qrf(8)*3 Qrf(5)/2 Qrf(9)*6 -Qrf(7)*3 -Qrf(6)*3];
out.Instance.U_DC_out = 2*[Qdc(8)*3 Qdc(5)/2 Qdc(9)*6 -Qdc(7)*3 -Qdc(6)*3];
out.Instance.Arf = Arf;
out.Instance.thetaRF = thetaRF;
out.Instance.Adc = Adc;
out.Instance.thetaDC = thetaDC;
T = [2 -2 0 0 0;...
    -2 -2 0 0 0;...
     0  4 0 0 0; ...
     0  0 1 0 0; ...
     0  0 0 1 0; ...
     0  0 0 0 1];
out.Instance.q = (1/V0)*T*out.Instance.U_RF_out';
out.Instance.alpha = (2/V0)*T*out.Instance.U_DC_out';
out.Instance.compensationError = [X(IDC)-XDC Y(JDC)-YDC Z(KDC)-ZDC]; 
if strcmp(operation,'findCompensation'),
    fprintf('Obsolete option!\n');
end
if qualityCheck
    out.Configuration.qualityRF = spher_harm_qlt(Vrf,Qrf,XRF,YRF,ZRF,7,X,Y,Z,'Spherical harmonic expansion: RF potential error');
    out.Configuration.qualityDC = spher_harm_qlt(Udc,Qdc,XDC,YDC,ZDC,7,X,Y,Z,'Spherical harmonic expansion: DC potential error');
end
print_underlined_message(' stop','post_process_trap');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end return values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Auxiliary functions 
%
    function f = d_e(Ei)
    % find the miscompensation distance, d_e, for the rf and dc potential 
    % given in the main program, in the presence of stray field Ei
        dm = Ei;
        E1 = dm(1); E2 = dm(2); E3 = dm(3);
        Vl = Udc-E1*x-E2*y-E3*z;
        [xrf yrf zrf] = exact_saddle(trap.Simulation.EL_RF,X,Y,Z,2,Zval);
        [xdc ydc zdc] = exact_saddle(Vl,X,Y,Z,3,Zval); 
        f = sqrt((xrf-xdc)^2+(yrf-ydc)^2+(zrf-zdc)^2);
    end
%
    function [fx fy fz theta Depth Xdc Ydc Zdc Xe Ye Ze U] = pfit(e1,e2,e3)
    % find the secular frequencies, tilt angle, and position of the dc 
    % saddle point for given combined input parameters. The stray field E 
    % has been defined in the body of the main function
        
        % find dc potential
        % keep until nikos fixes it %Vl = VDC1(w,n,cnt,ex,ey,ez);
        % remove % Vl = CalcVDC(trap,scale*dcVoltages,e1,e2,e3,NUM_DC,NUM_Center,x,y,z,truncVoltages,RF_offset);
        Vl = dc_potential(trap,dcVoltages,e1,e2,e3,x,y,z);
        [Xdc Ydc Zdc] = exact_saddle(Vl,X,Y,Z,3,Zval);
        
        % find pseudopotential
        [Ex,Ey,Ez] = gradient(Vrf,1e-3*grid(4),1e-3*grid(5),1e-3*grid(6));       
        Esq = Ex.^2 + Ey.^2 + Ez.^2;
        PseudoPhi = qe^2*Esq/(4*mass*Omega^2);    
        plot_potential(PseudoPhi/qe,Irf,Jrf,Krf,grid,pseudpotplot,'Pseudopotential','U_{ps} (eV)',varargin);
        
        % find total trap potential
        U = PseudoPhi+qe*Vl;     
        plot_potential(U/qe,Irf,Jrf,Krf,grid,trappotplot,'TrapPotential','U_{sec} (eV)',varargin);
        
        % find trap frequencies and tilt in radial directions
        Uxy = U(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf); 
        MU = max(max(Uxy));
        Uxy = Uxy/MU;
        dL = (x(Irf+3,Jrf,Krf)-x(Irf,Jrf,Krf));
        xr = (x(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-x(Irf,Jrf,Krf))/dL ;
        yr = (y(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-y(Irf,Jrf,Krf))/dL;
        [C1 C2 theta] = p2d(Uxy,xr,yr);                                    % look at TestMyFunctions for testiin        
        fx = (1e3/dL)*sqrt(2*C1*MU/(mass))/(2*pi); 
        fy = (1e3/dL)*sqrt(2*C2*MU/(mass))/(2*pi);
        
        % find trap frequency in axial direction
        Uz = projection(U,Irf,Jrf,3);
        l1 = max([Krf-6 1]);
        l2 = min([Krf+6 max(size(Z))]);
        if size(Z(l1:l2)-Z(Krf),1)~=size(Uz(l1:l2),1),
            Uz = Uz';
        end
        p = polyfit( (Z(l1:l2)-Z(Krf))/dL,Uz(l1:l2),6 );
        if 0,                                                              % set this to 1 to debug axial fit
            ((Z(l1:l2)-Z(Krf))/dL)'
            Uz(l1:l2)'
            ft = polyval(p,(Z-Z(Krf))/dL,6); 
            plot(Z,Uz); hold on; plot(Z(l1:l2),ft(l1:l2),'r'); 
            title('Potential in axial direction');
            xlabel('axial direction (mm)'); ylabel('trap potential (J)');hold off; pause(0.1);
        end
        fz= (1e3/dL)*sqrt(2*p(5)/(mass))/(2*pi);
        [Depth Xe Ye Ze] = trap_depth(U,X,Y,Z,Irf,Jrf,Krf);                
    end
%
    function out = normalize(in)
    % keep only the first 4 significant digits of the increment in vector
    % "in"
        dr = (in(size(in,1))-in(1))/(size(in,1)-1);
        p = 0; cnt = 0;
        while (cnt == 0)
            dr = 10*dr;
            cnt = fix(dr);
            p = p+1;
        end
        out = roundn(in,-p-4);       
    end
%   
    function f= warn(str,Vi,I,J,K)
    % print a warning if the saddle point is out of bounds   
        f = false;
        if (I == 1)||(I == size(Vi,1)),
                fprintf('%s saddle out of bounds (I=%i).\n',str,I);
                f = true;
        end
        if (J == 1)||(J == size(Vi,2)),
                fprintf('%s saddle out of bounds (J=%i).\n',str,J);
                f = true;
        end
        if (K == 1)||(K == size(Vi,3)),
                fprintf('%s saddle out of bounds (K=%i).\n',str,K);
                f = true;
        end  
    end
end        
