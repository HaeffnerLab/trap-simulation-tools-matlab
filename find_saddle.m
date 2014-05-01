function [Is Js Ks]=find_saddle(A,X,Y,Z,dim,Z0,varargin)
% [Is Js Ks]=find_saddle(A,grid,dim,Z0,varargin)
% returns the indices of the local extremum or saddle point of the scalar A. 
% A is  stored as  a MxNxP matrix if dim==3 and as a MxN matrix if dim==2.
% Z0 is the Z axis value on which the saddle point is evaluated, if dim==2.
% varargin are 2 optional arguments (both must be passed): 
% * a string indicating where find_saddle was called from
% * a boolean indicating whether it should be run in debug mode
% For dim ==2 the values of A are linearly extrapolated from [Z0] and [Z0]+1 
% to those corresponding to Z0 and Ks is such that Z(Ks)<Z0, Z(Ks+1)>=Z0.

debug = false;
threshold = 1e-3;
if length(varargin)>1,
    called_from_string = varargin{1};
    debug = varargin{2};
else
    called_from_string = '';
end

if (dim==3)
	if (numel(size(A)) ~= 3), 
		fprintf('%s: Problem with find_saddle.m dimensionalities.\n',called_from_string); 
		return; 
	end;
	
	f = A/max(max(max(A)));
	[qx qy qz] = gradient(f,X(2)-X(1),Y(2)-Y(1),Z(2)-Z(1));
	q = sqrt(qx.^2+qy.^2+qz.^2);
	%[m,Ks]=min(min(min(q)));
	%q=permute(q,[1 3 2]);
	%[m,Js]=min(min(min(q)));
	%q=permute(q,[3 2 1]);
	%[m,Is]=min(min(min(q)));
	m = q(1,1,1);
	Is = 1; Js = 1; Ks = 1;
	for i=1:size(q,1)
		for j=1:size(q,2)
		for k=1:size(q,3)
			if (q(i,j,k)<m)
			m = q(i,j,k);
			Is = i;
			Js = j;
			Ks = k;
			end
		end
		end
	end
	
	if debug,
		f = reshape(q,1,size(q,1)*size(q,2)*size(q,3));
		[f ind] = sort(f); 
		v = reshape(A,1,size(A,1)*size(A,2)*size(A,3)); v = v(ind);
		plot(f/max(f),'b'); hold on; plot(v/max(v),'r'); 
        legend('red: normalized potential','blue: normalized gradient of the potential');
        titl = sprintf('find_saddle called from %s.\nDisable by setting debug argument to false',called_from_string); 
        title(titl);
        hold off; pause; 
	end
	
	if (Is==1)||(Is==size(A,1)), 
		fprintf('%s/find_saddle: Saddle out of bounds in  x (i) direction.\n',called_from_string); 
	end
	if (Js==1)||(Js==size(A,2)), 
		fprintf('%s/find_saddle: Saddle out of bounds in  y (j) direction.\n',called_from_string); 
	end
	if (Ks==1)||(Ks==size(A,3)), 
		fprintf('%s/find_saddle: Saddle out of bounds in  z (k) direction.\n',called_from_string); 
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (dim==2)
	% extrapolate to the values of A at Z0
	A2 = A;
	if numel(size(A))==3,
		Ks = find(Z<=Z0,1,'last'); % this  makes Z(Ks)<Z0, Z(Ks+1)>=Z0
		if Ks==size(Z,2), 
            fprintf('Ks = last element!\n'); 
            return; 
        end;
		a1 = A(:,:,Ks);
		a2 = A(:,:,Ks+1);
		A2 = a1+(a2-a1)*(Z0-Z(Ks))/(Z(Ks+1)-Z(Ks));
	end
	if (size(size(A2),2) ~= 2), 
		fprintf('%s: Problem with find_saddle.m dimensionalities.\n',called_from_string); 
		return; 
	end;
	f = A2/max(max(A2));
	[qx qy] = gradient(f,X(2)-X(1),Y(2)-Y(1));
	q = sqrt(qx.^2+qy.^2);
    if 0,% add a plotting flag here?
        imagesc(q);
        title('find_saddle plot: (grad(U_{ps}))^2')
        pause(0.1);
    end
	m = min(min(q));
	if m>threshold
		if debug,
		Is = NaN; Js = NaN;
		fprintf('%s: Warning, there might be no saddle point\n',called_from_string);
		fprintf('Consider running find_saddle with debug = true and increasing the threshold parameter.\n');
		fprintf(sprintf('threshold = %3.1E\n',threshold));
		%return;
		end
	end
	mr = q(1,1);
	for i=1:size(q,1)
		for j=1:size(q,2)
		if (q(i,j)<mr)
			mr = q(i,j);
			Is = i;
			Js = j;
			%fprintf('val,Is,Js = %f,%i,%i\n',mr,i,j);
		end
		end
	end
	if (Is==1)||(Is==size(A2,1)), 
		fprintf('%s/find_saddle: Saddle out of bounds in  x (i) direction.\n',called_from_string); 
	end
	if (Js==1)||(Js==size(A2,2)), 
		fprintf('%s/find_saddle: Saddle out of bounds in  y (j) direction.\n',called_from_string); 
	end
% 	mf = min(min(f));
% 	if f(Is,Js)>mf,
% 		fprintf('find_saddle warning: I seem to be stuck in the escape point. Reverting to minimum position.\n');
% 		[mna jm] = min(f);
% 		[mna Js] = min(mna);
% 		Is = jm(Js);
% 	end % commented by Nikos on oct.11.2010: why had i been asking for the
% 	absolute minimum?
end

end

