function [D x y z] = trapdepth(V,X,Y,Z,Im,Jm,Km)
% [D x y z] = trapdepth(V,X,Y,Z,Im,Jm,Km)
% find the trap depth for trap potential V
% D is the trap depth
% x,y,z are the coordinates of the escape position
% X,Y,Z are vectors defining the grid in X,Y, and Z directions
% Im, Jm, Km are the indices of the trap potential minimum (ion position)

if (numel(size(V)) ~= 3), 
    fprintf('Problem with trapdepth dimensionalities.\n'); 
    return; 
end;

m = min(min(min(V(Im-1:Im+1,Jm-1:Jm+1,Km-1:Km+1))));
%V = max(V,2*m-V);

N1 = size(V,1); N2 = size(V,2); N3 = size(V,3);
[Ygr Xgr Zgr] = meshgrid(Y,X,Z);

f = V/max(max(max(V)));
[qx qy qz] = gradient(f,X(2)-X(1),Y(2)-Y(1),Z(2)-Z(1));
q = sqrt(qx.^2+qy.^2+qz.^2)/max(max(max(sqrt(qx.^2+qy.^2+qz.^2))));

f = reshape(f,[1 N1*N2*N3]);
q = reshape(q,[1 N1*N2*N3]); 
Xgr = reshape(Xgr,[1 N1*N2*N3]); 
Ygr = reshape(Ygr,[1 N1*N2*N3]); 
Zgr = reshape(Zgr,[1 N1*N2*N3]); 
v = reshape(V,[1 N1*N2*N3]); 

[srt ix] = sort(f);
%plot(srt); title('Finding trap depth. Sorted trap potential data shown.');pause(1);
qs = q(ix);
Xs = Xgr(ix);
Ys = Ygr(ix);
Zs = Zgr(ix);
vs = v(ix);
%plot(qs); title('Finding trap depth. Modulus of gradient of (sorted) trap potential shown.');
[der I] = min(qs(27:size(qs,2)));  % the nearest neighbors of the minimum cannot be in the escape range
if I==27, 
    fprintf('trapdepth.m:\nEscape point too close to trap minimum.\nImprove grid resolution, or extend grid.\n');
end
if der>0.2,
    fprintf('trapdepth.m:\nEscape point parameter too high.\nImprove grid resolution, or extend grid.\n');
end
D = vs(I+26)-m; x =Xs(I+26); y= Ys(I+26); z = Zs(I+26);

