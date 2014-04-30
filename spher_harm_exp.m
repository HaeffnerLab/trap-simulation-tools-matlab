function f=spher_harm_exp(V,Xc,Yc,Zc,Order,X,Y,Z)
% function f=spher_harm_exp(V,Xc,Yc,Zc,Order,X,Y,Z);
% V is the matrix containing an electri poiential (which must satisfy Laplace's equation)
% Xc,Yc,Zc are the coordinates of the center of the multipoles. 
% Order is the order of the expansion
% X,Y,Z are the vectors that define the grid in three directions
%
% This function expands the potential V in sspherical harmonics, i.e.:
% V=C00*Y00+C10*Y10+C11c*Y11c+C11s*Y11s+...
% here the Ynm are chosen to be real, and subscript c corresponds to 
% cos(m*phi) dependence, while s is sin(m*phi). The expansion is carried up
% to multipoles of order Order.
%
% The indices in V are V(i,j,k)<-> V(x,y,z). 
%
% The function returns the coefficients in the order:[C00 C10 C11c C11s ]'
% These correspond to the multipoles in cartesian coordinares: 
% [c z -x -y (z^2-x^2/2-y^2/2) -3zx -3yz 3x^2-3y^2 6xy ... 15(3xy^2-x^3) 15(y^3-3yx^2)         ...]
%  1 2  3  4       5             6    7     8       9  ...       15           16          17   ...
%
% Q(10)  0.5[2z^3-3(x^2+y^2)z]
% Q(11) -1.5[4xz^2-x(x^2+y^2)]
% Q(12) -1.5[4yz^2-(x^2+y^2)y]
% Q(13)  15[z(x^2-y^2)]
% Q(14)  30xyz
% Q(15)  15(3xy^2-x^3)
% Q(16)  15(y^3-3x^2y)
%
% Nikos January 2009


s=size(V); nx=s(1); ny=s(2); nz=s(3); %size(V)
[y x z] = meshgrid(Y-Yc,X-Xc,Z-Zc); %size(x)

x=reshape(x,1,nx*ny*nz); y=reshape(y,1,nx*ny*nz); z=reshape(z,1,nx*ny*nz);
r=sqrt(x.^2+y.^2+z.^2); rt=sqrt(x.^2+y.^2); theta=atan2(rt,z); phi=atan2(y,x);
scale = max(max(max(r)));
r = r/scale;

% make the spherical harmonic matrix in sequence of [Y00 Y10 Y11c Y11s Y20 Y21c Y21s...]
% In other words: Calculate the basis vectors of the sph. harm. expansion:
N=nx*ny*nz; Q=(1:N)'; Q(:)=1; 
for n=1:Order
	p=legendre(n,cos(theta));
	c=r.^n.*p(1,:); c=c'; Q=horzcat(Q,c);
    for m=2:n+1
		c=r.^n.*p(m,:).*cos((m-1)*phi); c=c'; Q=horzcat(Q,c);
        s=r.^n.*p(m,:).*sin((m-1)*phi); s=s'; Q=horzcat(Q,s);
    end;
end;

% the dc potential:
W=reshape(V,1,nx*ny*nz)';

% numerically invert, here the actual expansion takes place and we obtain
% the expansion coefficients M_{ji}.
f=Q\W;

% rescale to original units
i = 1;
for n = 1:Order
    for m = 1:2*n+1
        i = i+1;
        f(i) = f(i)/scale^(n);
    end
end

