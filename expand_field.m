function datout = expand_field(trap)
% function datout = expand_field(trap)
% Return a field datout.Configuration.multipoleCoefficients, which contains the multipole coefficients for all electrodes.
% Also (if regenerate = true) regenerate the DC potential trap for all electrodes using multipole
% expansion to order expansionOrder. 
% The electrodes are ordered as E(1), ...,  E(NUM_ELECTRODES-1)=E(Center), E(NUM_ELECTRODES)=E(RF)
% i.e. the NUM_ELECTRODES-1 is the center electrode bias, and the NUM_ELECTRODES is the RF electrode bias
%   
%                           ( multipoles    electrodes ->       )
%                           (     |                             )
% multipoleCoefficients =   (     V                             )
%                           (                                   )
%
% Multipole coefficients only up to order 8 are kept, but the
% coefficients are calculated up to order expansionOrder.
%
% trap: the simulation trap structure. 
% Xcorrection, Ycorrection: ADVANCED SETTING, leave at 0,0 most of the time
%                                      optional correction offsets from the RF saddle point, in case that was wrong by some known offset
% position: the axial position where the ion sits.
% expansionOrder: order of the multipole expansion
% NUM_ELECTRODES: number of DC electrodes and number of center electrodes

% Nikos June 2009
% Cleaned up 26-05-2013, 10-23-2013

% the correction Xcorrection,Ycorrection are parameters allowing one to offset the RF
% saddle point, for example to fix a wrong RF simulation
print_underlined_message('start','expand_field');
Xcorrection = trap.Configuration.Xcorrection;
Ycorrection = trap.Configuration.Ycorrection;
position = trap.Configuration.trappingPosition;
expansionOrder = trap.Configuration.expansionOrder;
NUM_ELECTRODES = trap.Configuration.NUM_ELECTRODES;
if isfield(trap.Configuration,'regeneratePotentialsCompleted')
    if regeneratePotentialsCompleted
        regenerate = false;
    elseif trap.Configuration.regeneratePotentials % in case the process was interrupted previously
        regenerate = true; 
    else
        regenerate = false;% do nothing
    end
else
    regenerate = trap.Configuration.regeneratePotentials;
end
regenerateCounter = 0;

fprintf('expand_field: Correction of XRF: %f mm.\n',Xcorrection);
fprintf('expand_field: Correction of YRF: %f mm.\n',Ycorrection);

datout = trap;
X = normalize(trap.Simulation.X); 
Y = normalize(trap.Simulation.Y); 
Z = normalize(trap.Simulation.Z);
[y x z] = meshgrid(Y,X,Z);

ord = zeros(1,NUM_ELECTRODES);
ord(:)=expansionOrder;

fprintf('________pre-smoothing the RF potential________\n');
%xpand the rf about the grid center, regenerate trap from the expansion
Irf = floor(numel(X)/2); Jrf = floor(numel(Y)/2); Krf = floor(numel(Z)/2);
Xrf = X(Irf); Yrf=Y(Jrf); Zrf=Z(Krf);
Qrf = spher_harm_exp(trap.Simulation.EL_RF,Xrf,Yrf,Zrf,expansionOrder,X,Y,Z);  
tempEL_RF = spher_harm_cmp(Qrf,Xrf,Yrf,Zrf,expansionOrder,X,Y,Z);
% expand the rf about its saddle point at the trapping position, save the quadrupole 
% components 
fprintf('________looking for RF saddle point________\n');
[Xrf Yrf Zrf] = exact_saddle(trap.Simulation.EL_RF,X,Y,Z,2,position);
Qrf = spher_harm_exp(tempEL_RF,Xrf+Xcorrection,Yrf+Xcorrection,Zrf,expansionOrder,X,Y,Z);  
datout.Configuration.Qrf = 2*[Qrf(8)*3 Qrf(5)/2 Qrf(9)*6 -Qrf(7)*3 -Qrf(6)*3];
datout.Configuration.thetarf = 45*(sign(Qrf(9)))-90*atan((3*Qrf(8))/(3*Qrf(9)))/pi;
if regenerate
    fprintf('________regenerating the RF potential________\n');
    datout.Simulation.EL_RF = spher_harm_cmp(Qrf,Xrf,Yrf,Zrf,expansionOrder,X,Y,Z);
    regenerateCounter = regenerateCounter+1;
end

E = [0 0 0];
c = [ 1  0  0  0  0  0  0  0  0; ...
      0  0  1  0  0  0  0  0  0; ...
      0  0  0  1  0  0  0  0  0; ...
      0 -1  0  0  0  0  0  0  0; ...
      0  0  0  0  0  0  0  6  0; ...
      0  0  0  0  1  0  0  0  0; ...
      0  0  0  0  0  0  0  0 12; ...
      0  0  0  0  0  0 -6  0  0; ...
      0  0  0  0  0 -6  0  0  0];
 M = zeros((expansionOrder+1)^2,NUM_ELECTRODES);

% expand all dc electrode potentials around RF null, save the coefficients
% and regenerate the trap if asked
for el = 1:NUM_ELECTRODES % Expand all the electrodes and  regenerate the potentials from the multipole coefficients
    multipoleDCVoltages = zeros(1,NUM_ELECTRODES);
    multipoleDCVoltages(el) = 1;
    Vdc = dc_potential(trap,multipoleDCVoltages,E(1),E(2),E(3),x,y,z);
    %plot_potential(Vdc,Irf,Jrf,Krf,trap.grid,'1d plots',sprintf('EL. %i DC Potential',el),'V (Volt)');
    Q = spher_harm_exp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
    M(:,el) = Q(1:(expansionOrder+1)^2);
    if regenerate,
        fprintf('________regenerating EL_DC%i potential________\n',el);
        trap.Simulation.(['EL_DC' num2str(el)]) = spher_harm_cmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);
        regenerateCounter = regenerateCounter+1;
    end
end

fprintf(sprintf('expand_field: Size of the raw multipole coefficient matrix is (%i,%i).\n',size(M,1),size(M,2)));
datout.Configuration.multipoleCoefficients = vertcat(c*M(1:9,:),M(10:(expansionOrder+1)^2,:)); 
if regenerate
    if regenerateCounter == NUM_ELECTRODES+1 % the RF is regenerated twice...
        datout.Configuration.regeneratePotentialsCompleted = true;
    else
        fprintf('Something went wrong with regenerating the potentials from the harmonic expansions.\n');
        datout.Configuration.regeneratePotentialsCompleted = false;
    end
else
   fprintf('Skipped regenerating the potentials from the harmonic expansions.\n');
   fprintf('Run expand_field after switching regeneratePotentials to true.\n');
   datout.Configuration.regeneratePotentialsCompleted = false; 
end
print_underlined_message(' stop','expand_field');

%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
 
end
