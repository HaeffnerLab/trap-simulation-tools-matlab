function datout = expand_field(data)
% function datout = expand_field(data)
% Return a field datout.trapConfiguration.multipoleCoefficients, which contains the multipole coefficients for all electrodes.
% Also (if regenerate = true) regenerate the DC potential data for all electrodes using multipole
% expansion to order expansionOrder. 
% The electrodes are ordered as E(1), ..., E(NUM_ELECTRODES)=E(RF)
% i.e. the NUM_ELECTRODES-1 is the center electrode bias, and the NUM_ELECTRODES is the RF electrode bias
% (if center and RF are used)  
%       ( multipoles    electrodes ->       )
%       (     |                             )
% M =   (     V                             )
%       (                                   )
% Multipole coefficients only up to order 8 are kept, but the
% coefficients are calculated up to order expansionOrder.
%
% data: the simulation data structure. 
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
regenerate = data.trapConfiguration.regeneratedPotentials;
Xcorrection = data.trapConfiguration.Xcorrection;
Ycorrection = data.trapConfiguration.Ycorrection;
position = data.trapConfiguration.trappingPosition;
expansionOrder = data.trapConfiguration.expansionOrder;
NUM_ELECTRODES = data.trapConfiguration.NUM_ELECTRODES;

fprintf('expand_field: Correction of XRF: %f mm.\n',Xcorrection);
fprintf('expand_field: Correction of YRF: %f mm.\n',Ycorrection);

datout = data;
X = normalize(data.trapConfiguration.X); 
Y = normalize(data.trapConfiguration.Y); 
Z = normalize(data.trapConfiguration.Z);
[y x z] = meshgrid(Y,X,Z);

ord = zeros(1,NUM_ELECTRODES);
ord(:)=expansionOrder;

% expand the rf about the grid center, regenerate data from the expansion
Irf = floor(numel(X)/2); Jrf = floor(numel(Y)/2); Krf = floor(numel(Z)/2);
Xrf = X(Irf); Yrf=Y(Jrf); Zrf=Z(Krf);
Qrf = spherharmxp(data.trapConfiguration.EL_RF,Xrf,Yrf,Zrf,expansionOrder,X,Y,Z);  
datout.trapConfiguration.EL_RF = spherharmcmp(Qrf,Xrf,Yrf,Zrf,expansionOrder,X,Y,Z);
%eval(sprintf('datout.trapConfiguration.%s = spherharmcmp(Qrf,Xrf,Yrf,Zrf,expansionOrder,X,Y,Z);','EL_RF'));

% expand the rf about its saddle point at the trapping position, save the quadrupole 
% components 
[Xrf Yrf Zrf] = exactsaddle(data.trapConfiguration.EL_RF,X,Y,Z,2,position);
Qrf = spherharmxp(data.trapConfiguration.EL_RF,Xrf+Xcorrection,Yrf+Xcorrection,Zrf,expansionOrder,X,Y,Z);  
datout.trapConfiguration.Qrf = 2*[Qrf(8)*3 Qrf(5)/2 Qrf(9)*6 -Qrf(7)*3 -Qrf(6)*3];
datout.trapConfiguration.thetarf = 45*(sign(Qrf(9)))-90*atan((3*Qrf(8))/(3*Qrf(9)))/pi;

E = [0 0 0];
%L = expansionOrder;
%ord = [L L L L L L L L L L L L L L L L L L L L L];
%ord =  [5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];
c = [ 1  0  0  0  0  0  0  0  0; ...
      0  0  1  0  0  0  0  0  0; ...
      0  0  0  1  0  0  0  0  0; ...
      0 -1  0  0  0  0  0  0  0; ...
      0  0  0  0  0  0  0  6  0; ...
      0  0  0  0  1  0  0  0  0; ...
      0  0  0  0  0  0  0  0 12; ...
      0  0  0  0  0  0 -6  0  0; ...
      0  0  0  0  0 -6  0  0  0];
 M1 = zeros((expansionOrder+1)^2,NUM_ELECTRODES);

% expand all dc electrode potentials around RF null, save the coefficients
% and regenerate the data if asked
for el = 1:(NUM_ELECTRODES)
  if data.trapConfiguration.electrodeMapping(el,2)
    multipoleDCVoltages = zeros(1,NUM_ELECTRODES);
    multipoleDCVoltages(el) = 1;
    manualDCVoltages = zeros(1,NUM_ELECTRODES);
    Vdc = dcpotential(data,multipoleDCVoltages,manualDCVoltages,E(1),E(2),E(3),x,y,z);
    %plotpot(Vdc,Irf,Jrf,Krf,data.grid,'1d plots',sprintf('EL. %i DC Potential',el),'V (Volt)');
    Q = spherharmxp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
    M1(:,el) = Q(1:(expansionOrder+1)^2);
    %max(Q(2:(expansionOrder+1)^2))
    if regenerate,
        if isfield(data.trapConfiguration,['EL_DC' num2str(el)]),
            data.trapConfiguration.(['EL_DC' num2str(el)]) = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);
            % old command 12-10-2013: eval(sprintf('datout.%s = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);',str));
        end
    end
  elseif data.trapConfiguration.manualElectrodes(el)
    if regenerate,
        data_loc = data;
        data_loc.trapConfiguration.manualElectrodes(el) = 1;
        multipoleDCVoltages = zeros(1,NUM_ELECTRODES);
        manualDCVoltages = zeros(1,NUM_ELECTRODES);
        manualDCVoltages(el)  = 1;
        Vdc = dcpotential(data_loc,multipoleDCVoltages,manualDCVoltages,E(1),E(2),E(3),x,y,z);
        %plotpot(Vdc,Irf,Jrf,Krf,data.grid,'1d plots',sprintf('El. %i DC Potential',el),'V (Volt)');
        Q = spherharmxp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);
        data.trapConfiguration.(['mEL_DC' num2str(el)]) = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);
        % old command 12-10-2013: eval(sprintf('datout.%s = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);',str));  
    end
  end
end

fprintf(sprintf('expand_field: Size of the raw multipole coefficient matrix is (%i,%i).\n',size(M1,1),size(M1,2)));
datout.trapConfiguration.multipoleCoefficients = vertcat(c*M1(1:9,:),M1(10:(expansionOrder+1)^2,:)); 
fprintf('expand_field: Size of the resized multipole coefficient matrix is (%i,%i).\n',size(vertcat(c*M1(1:9,:),M1(10:(expansionOrder+1)^2,:)),1),size(vertcat(c*M1(1:9,:),M1(10:(expansionOrder+1)^2,:)),2));
fprintf('expand_field: ended successfully.\n');
print_underlined_message('stop_','expand_field');

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
