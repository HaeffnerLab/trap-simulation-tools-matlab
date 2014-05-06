% example of pairing up and inversion
%       ( multipoles    electrodes ->       )
%       (     |                             )
% M =   (     V                             )
%       (                                   )
%
%            ( electrodes    multipoles ->       )
%            (     |                             )
% C =        (     V                             )
%            (                                   )
clear all;
regularize = false;

% five electrodes, three multipoles
M1 = [1 0 -1 1 1 0;...
      0 1 -1 0 0 0;...
      0 0  1 1 0 0];
I = eye(3);
C1 = M1\I;
if regularize
    K = null(M1);
    for ii = 1:3
        Cv = C1(ii,:)';
        lmbda = K\Cv;
        C1(ii,:) = C1(ii,:)-(K*lmbda)';
    end
end


% five electrodes, two are paired up, three multipoles
M2 = [1 0 -1 2 0 0;...
      0 1 -1 0 0 0;...
      0 0  1 1 0 0];
C2 = M2\I;
if regularize
    K = null(M2);
    for ii = 1:3
        Cv = C2(ii,:)';
        lmbda = K\Cv;
        C2(ii,:) = C2(ii,:)-(K*lmbda)';
    end
end