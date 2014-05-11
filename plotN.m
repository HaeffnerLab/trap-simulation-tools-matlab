function plotN(V,titleString)
% mesh the values of the dc voltages corresponding to the N electrodes
% of a planar trap, in a geometrically "correct" way
% Assume that N-2 are DC electrodes
% V is a vector of N elements
% titleString is the title to place on the plot
% Nikos, July 2009

N = max(size(V));
Ndc = N-2;
A = zeros(10*(N-2)/2,90);
% DC electrode segments
for i=1:Ndc/2
    A(10*(i-1)+1:10*i,21:30) = V(i);
    A(10*(i-1)+1:10*i,61:70) = V(Ndc/2+i);
end
% RF
A(:,31:40) = V(N)*ones(10*(N-2)/2,10);
A(:,51:60) = V(N)*ones(10*(N-2)/2,10);
% Center
A(:,41:50) = V(N-1)*ones(10*(N-2)/2,10);
mesh(A);
axis([0 size(A,1) 0 size(A,2) min(0,1.2*min(V)) max(0,1.2*max(V))]);
title(titleString);
