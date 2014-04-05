function V = spherharmcmp(C,Xc,Yc,Zc,Order,Xe,Ye,Ze)
% function V=spherharmcmp(C,Xc,Yc,Zc,Order,Xe,Ye,Ze);
% This function computes the potential V from the spherical harmonic 
% coefficients, i.e.: V=C1*Y00+C2*Y10+C3*Y11c+C4*Y11s+...
% here the Ynm are chosen to be real, and subscript c corresponds to 
% cos(m*phi) dependence, while s is sin(m*phi). The expansion is carried up
% to multipoles of order Order. If the size of the coefficient vector C is
% not Order^2, a warning message is displayed
%
% The indices in V are V(i,j,k)<-> V(x,y,z). 
%
% C = [C1 C2 ...]'  the vector of coefficients.
% Xc,Yc,Zc:         the coordinates of the center of the multipoles. 
% Order:            the order of the expansion.
% Xe,Ye,Ze:         the vectors that define the grid in three directions.
%
% The input coefficients are given in the order:[C00 C10 C11c C11s ]'
% These correspond to the multipoles in cartesian coordinares: 
% [c z -x -y (z^2-x^2/2-y^2/2) -3zx -3yz 3x^2-3y^2 6xy]
%  1 2  3  4       5             6    7     8       9
% Nikos June 2009

V = [];
%size(C)
%return;
if (size(C,1)~=(Order+1)^2),
    while 1
        st = input('spherharrmcmp.m warning:\nSize of coefficient vector not equal to Order^2. Proceed? (y/n)\n','s');
        if strcmp(st,'n'), return;
        elseif strcmp(st,'y'), break;
        end
    end
end
[y x z] = meshgrid(Ye-Yc,Xe-Xc,Ze-Zc);
s=size(x); nx=s(1); ny=s(2); nz=s(3);

x=reshape(x,1,nx*ny*nz); y=reshape(y,1,nx*ny*nz); z=reshape(z,1,nx*ny*nz);
r=sqrt(x.^2+y.^2+z.^2); rt=sqrt(x.^2+y.^2); theta=atan2(rt,z); phi=atan2(y,x);

% make the spherical harmonic matrix in sequence of [Y00 Y10 Y11c Y11s
% Y20 Y21c Y21s...]
N=nx*ny*nz; Q=(1:N)'; Q(:)=1; 
for n=1:Order
	p=legendre(n,cos(theta));
	c=r.^n.*p(1,:); c=c'; Q=horzcat(Q,c);
    for m=2:n+1
		c=r.^n.*p(m,:).*cos((m-1)*phi); c=c'; Q=horzcat(Q,c);
        s=r.^n.*p(m,:).*sin((m-1)*phi); s=s'; Q=horzcat(Q,s);
    end;
end;
%W=reshape(V,1,nx*ny*nz)';
size(Q);
% now multiply and reshape to the original
W = C'*Q';
V = reshape(W,nx,ny,nz);
