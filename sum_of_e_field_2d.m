function f = sum_of_e_field_2d(r,z0,A,X,Y,Z)
%function f = sum_of_e_field_2d(r,A,X,Y,Z)
% r: x,y center position for the spherical harmonic expansion
% z0: z center position for the spherical harmonic expansion
% A: potential matrix
% X,Y,Z: vectors specifying the grid
% find the weight of high order multipole terms compared to the weight 
% of second order multipole terms in matrix V, when the center of the
% multipoles is at x0,y0,z0
% Used by exact_saddle for 2-d saddle search
	x0 = r(1); y0 = r(2);
	c = spher_harm_exp(A,x0,y0,z0,4,X,Y,Z);
	s = c.^2; 
	f = (sum(s(2:4)))/sum(s(5:9));
end
