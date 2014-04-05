function Aout = projection(Ain,index1,index2,dim)
% Aout = projection(Ain,index1,index2,dim)
% Project along dimension dim, on the line passing through 
% index1, index2.
% Ain is a i,j,k matrix (3 dimensions)
% dim = 1,2,3 specifies the dimention of the projection
% index1 and index2 are always part of a cyclic permutation 
% of i,j,k 

if (max(dim == 1:3) == 0), 
    Aout = NaN; 
    return;
end
reorder = [1 2 3; 2 3 1; 3 1 2];
a = permute(Ain,reorder(dim,:));
Aout=a(:,index1,index2);