function f=spherharmq(V,C,Xc,Yc,Zc,Order,Xe,Ye,Ze,tit)
% function f=spherharmq(V,C,Xc,Yc,Zc,Order,Xe,Ye,Ze);
% This function determines the "quality" of the expansion of potential V in 
% spherical harmonics (V=C00*Y00+C10*Y10+C11c*Y11c+C11s*Y11s+... )
% here the Ynm are chosen to be real, and subscript c corresponds to 
% cos(m*phi) dependence, while s is sin(m*phi). The expansion is carried up
% to multipoles of order Order.
%
% The indices in V are V(i,j,k)<-> V(x,y,z). 
%
% V 
%   is the expanded potential
% C 
%   is the coefficient vector
% Xc,Yc,Zc 
%   are the coordinates of the center of the multipoles. 
% Order 
%   is the order of the expansion
% Xe,Ye,Ze 
%   are the vectors that define the grid in three directions
% tit 
% is a string describing the input potential for plot purposes
% ('RF', 'DC', etc.). if title=='noplots' no plots are made
%
% The coefficients are in the order:[C00 C10 C11c C11s ]'
% These correspond to the multipoles in cartesian coordinares: 
% [c z -x -y (z^2-x^2/2-y^2/2) -3zx -3yz 3x^2-3y^2 6xy]
%  1 2  3  4       5             6    7     8       9
% Nikos January 2009


s=size(V); nx=s(1); ny=s(2); nz=s(3);
Vfit = spherharmcmp(C,Xc,Yc,Zc,Order,Xe,Ye,Ze);
dV = V-Vfit;
e = reshape(dV,1,nx*ny*nz)/( max(max(max(V))) - min(min(min(V))) );
f = [max(abs(e)) mean(abs(e)) median(abs(e))];
if strcmp(tit,'noplots'), return; end
plot((abs(e))); title(tit); pause(0.2);