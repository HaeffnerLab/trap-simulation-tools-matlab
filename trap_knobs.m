function datout = trap_knobs(trap,plotOption)
% function datout = trap_knobs(trap,plotOption)
% update trap.Configuration with the matrix which controls the independent multipoles, and the kernel matrix
% Start from the matrix multipoleCoefficients, return a field multipoleControl with the linear combimations of trap 
% electrode voltages that give 1 V/mm, or 1 V/mm^2 of the multipole number i.
% Also return matrix multipoleKernel which is the kernel matrix of electrode linear combinations which do nothing to 
% the multipoles
% The order of multipole coefficients is:
% 1/r0 ^ [ x y z ] and 
% 1/r0^2* [ (x^2-y^2)/2 (2z^2-x^2-y^2)/2 xy/2 yz/2 xz/2 ], where r0 is 1 mm
% (unless rescaling is applied)
%
%                           ( multipoles    electrodes ->       )
%                           (     |                             )
% multipoleCoefficients =   (     V                             )
%                           (                                   )
%
%                           ( electrodes    multipoles ->       )
%                           (     |                             )
% multipoleControl =        (     V                             )
%                           (                                   )
%
% trap:         input trap structure
% position:     position of the ion inside the trap
% plotOption:   show plots, 1 = yes
% regularize:   sets whether regularization is done or not. 1 = yes
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
debug = true;
datout = trap;
NUM_ELECTRODES = trap.Configuration.NUM_ELECTRODES;
electrodeMapping = trap.Configuration.electrodeMapping;
manualElectrodes = trap.Configuration.manualElectrodes;
expansionOrder =  trap.Configuration.expansionOrder;
numTotMultipoles = length(trap.Configuration.usedMultipoles);
usedMultipoles = trap.Configuration.usedMultipoles;
numUsedMultipoles = sum(trap.Configuration.usedMultipoles);
multipoleCoefficients = trap.Configuration.multipoleCoefficients;
regularize = trap.Configuration.regularize;

MR = compact_matrix(multipoleCoefficients, NUM_ELECTRODES, (expansionOrder+1)^2, electrodeMapping, manualElectrodes);

datout.Configuration.multipoleCoefficientsReduced = vertcat(MR(1:9,:),MR(10:(expansionOrder+1)^2,:));
fprintf('trap_knobs: with electrode constrains, the coefficient matrix size is (%i,%i).\n',...
            size(vertcat(MR(1:9,:),MR(10:(expansionOrder+1)^2,:)),1),...
            size(vertcat(MR(1:9,:),MR(10:(expansionOrder+1)^2,:)),2));

allM = vertcat(MR(2:1+numTotMultipoles,:));%,eye(21));
imagesc(allM); 
colormap(cool);
title('trap_ knobs: Multipole coefficients (Ex...U5) for each trap electrode. cyan:(-) / red:(+)');
pause;
C = zeros(numUsedMultipoles,size(allM,2));
usedM = zeros(numUsedMultipoles,size(allM,2));
usmm = 1;
for mm = 1:numTotMultipoles % keep only the multipoles you specified in usedMultipoles
  if trap.Configuration.usedMultipoles(mm),
    usedM(usmm,:) = allM(mm,:);
    usmm = usmm+1;
  end
end
for ii=1:numUsedMultipoles
        multipoleBasisVector = zeros(numUsedMultipoles,1);
        multipoleBasisVector(ii) = 1;
        singleMultipoleControl = usedM\multipoleBasisVector;
        multipoleOutputVector = usedM*singleMultipoleControl; 
        err = multipoleBasisVector-multipoleOutputVector;
        if plotOption
            plot(err); 
            title(sprintf('trap_ knobs: difference between set multipoles and achieved multipoles.\nIntended output at multipole %i',ii));
            pause;
            plotN(expand_matrix_el(singleMultipoleContro, NUM_ELECTRODES, electrodeMapping, manualElectrodes));
            pause;
        end
        C(ii,:) = singleMultipoleControl';
end

K = null(usedM);

if regularize
    for ii = 1:numUsedMultipoles
        Cv = C(ii,:)';
        lmbda = K\Cv;
        C(ii,:) = C(ii,:)-(K*lmbda)';
    end
end
if debug
    datout.Configuration.multipoleControlReduced = C';
end
C = expand_matrix_mult(C,numTotMultipoles,usedMultipoles);
C = expand_matrix_el(C, NUM_ELECTRODES, electrodeMapping, manualElectrodes); 

datout.Configuration.multipoleKernel = K;
datout.Configuration.multipoleControl = C';
print_underlined_message(' stop','trap_knobs');

%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function RM = compact_matrix(MM, NUM_ELECTRODES,numMultipoles, electrodeMap, manualEl) 
        % multipole compaction operation: combine paired electrodes and remove
        % manually controlled electrodes form multipoleCoefficients matrix
        % to test this function:
        %clear all;
        %numMultipoles = 5;
        %NUM_ELECTRODES = 9;
        %electrodeMap = [1 1; 2 1; 3 2; 4 2; 5 3; 6 4; 7 5; 8 6; 9 7];
        %manualEl = [1 1 0 0 0 0 0 0 1];
        %MM = [1 0 1 0 1 0 1 0 1;...
        %     0 0 1 0 0 1 0 0 1;...
        %     1 0 0 1 0 0 1 0 0;...
        %     0 1 0 0 1 0 0 1 0;...
        %     1 0 0 1 0 0 0 1 0];
        % %then copy the following code to command line 
        MR1 = zeros(numMultipoles,electrodeMap(NUM_ELECTRODES,2));
        for ell = 1:NUM_ELECTRODES % combine paired electrodes and set manually controlled electrodes to zero
            MR1(:,electrodeMap(ell,2)) = MR1(:,electrodeMap(ell,2)) + (manualEl(ell)==0)*MM(:,ell);
            fprintf('trap_knobs/compact_matrix/combine electrodes: %i->%i\n',ell,electrodeMap(ell,2));
        end
        RM = zeros(numMultipoles,electrodeMap(NUM_ELECTRODES,2)-sum((manualEl==1)));
        RMel = 1;
        for ell = 1:size(MR1,2) % remove the manually controlled electrodes
            if (max(abs(MR1(:,ell)))>0)
                RM(:,RMel) = MR1(:,ell);
                fprintf('trap_knobs/compact_matrix/keep multipole controlled electrodes: %i\n',ell);
                RMel = RMel+1;
            end
        end
    end

    
    function em = expand_matrix_mult(RM,numTotMultipoles,usedMultipoles)
        % expand the multipole control martix to cover all the multipoles 
        % add zero rows to the positions where the multipoles are not
        % constrained by the code
        % RM is the reduced matrix
        % EM (output) is the expanded matrix
        em = zeros(numTotMultipoles,size(RM,2));
        currmult = 1;
        for cc = 1:numTotMultipoles
            if usedMultipoles(cc)
                em(cc,:) = RM(currmult,:);
                currmult = currmult+1;
            end
        end
    end

    function EM = expand_matrix_el(RM, NUM_ELECTRODES, electrodeMap, manualEl) 
        % expand a multipole control matrix from the functional electrode 
        % basis to the physical electrode basis. First step is to put back 
        % the grounded electrodes as 0s. Second step is to split up the 
        % paired up electrodes into their constituents.
        % RM is the reduced matrix
        % EM is the expanded matrix
        EM = zeros(numTotMultipoles,NUM_ELECTRODES);
        Nfunctional = 1;
        if manualEl(1) == 0
            EM(:,1) = RM(:,1);
            if electrodeMap(1,2) < electrodeMap(2,2)
                Nfunctional = Nfunctional+1;
            end
        end 
        for ee = 2:NUM_ELECTRODES
            if manualEl(ee) == 0
                EM(:,ee) = RM(:,Nfunctional);
                if ee < NUM_ELECTRODES
                    % Nfunctional increases only when the electrode is in
                    % multipole control and the map changes
                    if electrodeMap(ee,2) < electrodeMap(ee+1,2) 
                        Nfunctional = Nfunctional+1;
                    end
                end
            end
        end
    end

end




