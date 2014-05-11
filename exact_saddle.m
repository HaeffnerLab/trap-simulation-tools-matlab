function [Xs Ys Zs] = exact_saddle(V,X,Y,Z,dim,Z0)
% [Xs Ys Zs] = exact_saddle(V,X,Y,Z,dim,Z0)
% This version finds the approximate saddle point using the
% pseudopotential, does a multipole expansion around it, and finds the 
% exact saddle point by maximizing the quadrupole terms
% Data are stored in matrix V. The axes are (by convention): 
% x = radial horizontal -> I
% y = radial vertical   -> J
% z = axial             -> K
% X,Y,Z are vectors defining the grid in the three directions.
% dim is the dimensionality 2, or 3
% Z0 is the coordinate where a saddle point will be sought if dim==2
% Nikos Daniilidis 9-1-09
% Had issues with function nesting and variable scope definitions in Octave
% Revisited for Octave Compatibilty 5-25-13
% Needs Octave >3.6, packages general, miscellaneous, struct, optim

Xs = 0; Ys = 0; Zs = 0; 

if dim==3,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize 
    grid = [X(1) Y(1) Z(1) X(2)-X(1) Y(2)-Y(1) Z(2)-Z(1)];
    [I J K] = find_saddle(V,X,Y,Z,3,Z0,'exact_saddle',true);
    if (I<3)||(I>size(V,1)-2), 
        fprintf('exact_saddle.m: Saddle point out of bounds in radial direction.\n');
        return; 
    end
    if (J<3)||(J>size(V,2)-2), 
        fprintf('exact_saddle.m: Saddle point out of bounds in vertical direction.\n');
        return; 
    end
    if (K<3)||(K>size(V,3)-2), 
        fprintf('exact_saddle.m: Saddle point out of bounds in axial direction.\n');
        return; 
    end
  
    A = V(I-2:I+2,J-2:J+2,K-2:K+2); 
    xo = X(I); yo = Y(J); zo = Z(K);
    Xn = double(X(I-2:I+2)); Yn = double(Y(J-2:J+2)); Zn = double(Z(K-2:K+2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Minimize
    
    r = fminunc(@(r) sum_of_e_field (r,V,X,Y,Z),[xo, yo, zo]);

    Xs = r(1); 
    Ys = r(2); 
    Zs = r(3); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dim==2,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize
    if numel(size(V))==3,
        K = find(Z<Z0,1,'last');                                           % Extrapolate to value Z0
        if K==size(Z,2), return; end;
        v1 = V(:,:,K);
        v2 = V(:,:,K+1);
        V2 = v1+(v2-v1)*(Z0-Z(K))/(Z(K+1)-Z(K));
    end
    [I J] = find_saddle(V2,X,Y,Z,2,Z0,'exact_saddle',true);
    if (I<3)||(I>size(V,1)-2), 
        fprintf('exact_saddle.m: Saddle point out of bounds in radial direction.\n');
        return; 
    end
    if (J<3)||(J>size(V,2)-2), 
        fprintf('exact_saddle.m: Saddle point out of bounds in vertical direction.\n');
        return; 
    end
    for k=1:9,
        A(:,:,k) = V2(I-4:I+4,J-4:J+4);
    end
    xo = X(I); yo = Y(J); z0 = Z0;
    Xn = X(I-4:I+4); Yn = Y(J-4:J+4); Zn = Z(K-4:K+4);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Minimize
    h = @D2;
    r = fminunc(@(r) sum_of_e_field_2d (r,Z0,V,X,Y,Z),[xo, yo]);

    Xs = r(1); 
    Ys = r(2); 
    Zs = Z0; 
end


%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sum_of_e_field
% sum_of_e_field_2d
    
end


