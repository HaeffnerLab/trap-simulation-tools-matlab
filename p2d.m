function [Af Bf theta]=p2d(V,x,y)
% fit a 2d polynomial to the trap in V. 
% x and y are the 2d coordinate matrices on which the function is defined
% Af Bf and theta are the curvaturess in the Xr axis, Yr axes, where theta 
% is the angle between X and Xr. Not sure if the angle is correct when x
% and y are not centered around zero

N=size(V,2);
con=x; con(:,:)=1;
xx=x.*x; yy=y.*y; xy=x.*y; xxx = xx.*x; yyy = yy.*y; xxy = xx.*y; xyy = x.*yy;
xxxx = xx.*xx; yyyy = yy.*yy; xxxy = xxx.*y; xxyy = xx.*yy; xyyy = x.*yyy;

V2=reshape(V,1,N^2)';
xxxx2 =reshape(xxxx,1,N^2)'; yyyy2 =reshape(yyyy,1,N^2)'; xxxy2 =reshape(xxxy,1,N^2)';
xxyy2 =reshape(xxyy,1,N^2)'; xyyy2 =reshape(xyyy,1,N^2)'; xxx2 =reshape(xxx,1,N^2)';
yyy2 =reshape(yyy,1,N^2)'; xxy2 =reshape(xxy,1,N^2)'; xyy2 =reshape(xyy,1,N^2)';
xx2=reshape(xx,1,N^2)'; yy2=reshape(yy,1,N^2)'; xy2=reshape(xy,1,N^2)';
x2=reshape(x,1,N^2)'; y2=reshape(y,1,N^2)'; con2=reshape(con,1,N^2)';
Q=horzcat(xxxx2,yyyy2,xxxy2,xxyy2,xyyy2,xxx2,yyy2,xxy2,xyy2,xx2,yy2,xy2,x2,y2,con2);
c=Q\V2;

%mesh(x,y,V); pause(0.5);
Vfit = reshape(Q*c,N,N);
%mesh(x,y,(V-Vfit)/max(max(V)));
%title('Normalized error of the 2d fit');
%pause(0.5);

theta=-0.5*atan(c(12)/(c(11)-c(10)));
Af=0.5*(c(10)*(1+1/cos(2*theta))+c(11)*(1-1/cos(2*theta)));
Bf=0.5*(c(10)*(1-1/cos(2*theta))+c(11)*(1+1/cos(2*theta)));
theta=180*theta/pi;

