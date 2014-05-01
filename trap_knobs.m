function datout = trap_knobs(data,plotOption)
% function datout = trap_knobs(data,plotOption)
% update data.trapConfiguration with the matrix which controls the independent multipoles, and the kernel matrix
% Start from the matrix multipoleCoefficients, return a field multipoleControl with the linear combimations of trap 
% electrode voltages that give 1 V/mm, or 1 V/mm^2 of the multipole number i.
% Also return matrix multipoleKernel which is the kernel matrix of electrode linear combinations which do nothing to 
% the multipoles
% The order of multipole coefficients is:
% 1/r0 ^ [ x y z ] and 
% 1/r0^2* [ (x^2-y^2)/2 (2z^2-x^2-y^2)/2 xy/2 yz/2 xz/2 ], where r0 is 1 mm
% (unless rescaling is applied)
%

% data:     input data struct
% position: position of the ion inside the trap
% plotOption:  show plots, 1 = yes
% regularize:      sets whether regularization is done or not. 1 = yes
%
%
%	this is just a reminder from expand_field
% c = [ 0  1  0  0  0  0  0  0; ...
%       0  0  1  0  0  0  0  0; ...
%      -1  0  0  0  0  0  0  0; ...
%       0  0  0  0  0  0  6  0; ...
%       0  0  0  1  0  0  0  0; ...
%       0  0  0  0  0  0  0 12; ...
%       0  0  0  0  0 -6  0  0; ...
%       0  0  0  0 -6  0  0  0];

print_underlined_message('start','trap_knobs');
datout = data;
multipoleCoefficients = data.trapConfiguration.multipoleCoefficients;
regularize = data.trapConfiguration.regularize;

Mt = vertcat(multipoleCoefficients(2:9,:));%,eye(21));
imagesc(Mt); 
title('trap_ knobs: Multipole coefficients (Ex-U5) for trap electrodes');
pause;
numMultipoles = sum(data.trapConfiguration.usedMultipoles);
C = zeros(numMultipoles,size(Mt,2));
Mtt = zeros(numMultipoles,size(Mt,2));
for mm = 1:8
  if data.trapConfiguration.usedMultipoles(mm),
    Mtt(mm,:) = Mt(mm,:);
  end
end
for ii=1:numMultipoles
        Mf = zeros(numMultipoles,1);
        Mf(ii) = 1;
        P = Mtt\Mf;
        Mout = Mtt*P; err = Mf-Mout;
        if plotOption
            plot(err); title('Error of the fit elements');
            pause;
            plotN(P);
            pause;
        end
        %imagesc21(P,'');
        C(ii,:) = P';
end

K = null(Mtt);

if regularize
    for ii = 1:numMultipoles
        Cv = C(ii,:)';
        lambda = K\Cv;
        C(ii,:) = C(ii,:)-(K*lambda)';
    end
end

datout.trapConfiguration.multipoleKernel = K;
datout.trapConfiguration.multipoleControl = C';
print_underlined_message('stop_','trap_knobs');




