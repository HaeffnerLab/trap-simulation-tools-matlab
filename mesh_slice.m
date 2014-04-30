function mesh_slice(A,n,grid)
% mesh_slice(A,n,grid)
% plot successive slices of matrix A in the direction n=1(I),2(J),3(K)

order = [2 3 1; 3 1 2; 1 2 3];
q = permute(A,order(n,:));

i = grid(1+mod(n,3)) : grid(4+mod(n,3)) : ...
    ( grid(1+mod(n,3))+(size(A,1+mod(n,3))-1)*grid(4+mod(n,3)) );
j = grid(1+mod(n+1,3)) : grid(4+mod(n+1,3)) : ...
    ( grid(1+mod(n+1,3))+(size(A,1+mod(n+1,3))-1)*grid(4+mod(n+1,3)) );

for ii=1:size(q,3)
    slice = q(:,:,ii);
    mesh(j,i,slice); title(ii);
    if n==1, xlabel('horizontal axial (mm)'); ylabel('height (mm)'); end;
    if n==2, ylabel('horizontal axial (mm)'); xlabel('horizontal radial (mm)'); end;
    if n==3, ylabel('horizontal radial (mm)'); xlabel('height (mm)'); end;
    pause;
end;