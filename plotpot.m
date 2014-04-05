function plotpot(A,I,J,K,grid,func,titleString,yLabel,varargin) % varargin = {outPath,pathSuffix}
% plotpot(A,I,J,K,grid,f,titleString,yLabel,varargin)
% Make 2d mesh plots and 1d plots of A around the point I,J,K
% remember that I->X, J->Y, K->Z
%
% A: 
%   input 3 dimensional array
% I,J,K: 
%   index around which the plots are made, e.g. plot A(:,J,K) vs X
% grid: 
%   summary of the grid information 
%   grid = [Xmin Ymin Zmin Xstep Ystep Zstep]
% func : 
%   What kind of plot to produce. Takes one of the values in:
%   ('no plots', '2d plots','1d plots','2d and 1d plots') 
% titleString: 
%   the title of the produced plots
% yLabel:
%   the label on the y axis of the produced plot
% varargin: 
%   accepts two optional arguments used to save pdf figures of plots
%   outPath: path to save figures to. This is usually set to matDataPath
%   pathSuffix: suffix to add to outPath
%
% Nikos 2009
% Cleaned up by Nikos: March 2014

if strcmp(func,'no plots'), return; end;

if ( strcmp(func,'2d plots')||strcmp(func,'2d plots and 1d plots') ),
    meshslice(A,1,grid);
    meshslice(A,2,grid);
    meshslice(A,3,grid);
end

if ( strcmp(func,'1d plots')||strcmp(func,'2d plots and 1d plots') ),
    % plot along I
    PrI = A(:,J,K); 
    t = grid(1):grid(4):grid(1)+(size(A,1)-1)*grid(4);
    close all; 
    h = figure;
    axes1=axes('Parent',h,'FontSize',12);
    box(axes1,'on');
    hold(axes1,'all');
    plot(t,PrI); plot(t(I),PrI(I),'r*'); 
    xlabel('x (mm)','FontSize',12); ylabel(yLabel,'Interpreter','tex');
    title(titleString);
    xlim([min(t) max(t)])
    %print(h,'-depsc',[sprintf('%s%s/',outpath,newfilename) 'pot_radial.eps'])
    %saveas(h,[sprintf('%s%s/',outpath,newfilename) 'pot_radial.fig'])
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0.1 3.2 2.4])
    set(gcf, 'PaperSize', [3.1 2.4]);
    %print(h,'-dpdf',[outpath 'eps/' sprintf('%i',newfilename) sprintf('_%s_radial.pdf',titleString)])
    %saveas(h,[outpath 'fig/' sprintf('%i',newfilename) sprintf('_%s_radial.fig',titleString)])
    hold off; 
    pause;
    
    % plot along J
    PrJ = A(I,:,K); 
    t = grid(2):grid(5):grid(2)+(size(A,2)-1)*grid(5);
    close all;     h = figure;
    axes1=axes('Parent',h,'FontSize',12);  
    box(axes1,'on');
    hold(axes1,'all');
    plot(t,PrJ); plot(t(J),PrJ(J),'r*');
    xlabel('y (mm)','FontSize',12); ylabel(yLabel,'Interpreter','tex');
    title(titleString);
    xlim([min(t) max(t)])
    if length(varargin{1})>0,
        outPath = varargin{1};
        pathSuffix = varargin{2};
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0.1 3.2 2.4])
        set(gcf, 'PaperSize', [3.1 2.4]);
        print(h,'-dpdf',[outPath 'eps/' sprintf('%i',pathSuffix) sprintf('_%s_height.pdf',titleString)])
        saveas(h,[outPath 'fig/' sprintf('%i',pathSuffix) sprintf('_%s_height.fig',titleString)])
    end
    hold off; 
    pause;
    
    % plot along K
    aaa = permute(A,[3 1 2]);
    PrK = aaa(:,I,J);
    t = grid(3):grid(6):grid(3)+(size(A,3)-1)*grid(6);
    close all;     
    h = figure;
    axes1=axes('Parent',h,'FontSize',12);
    box(axes1,'on');
    hold(axes1,'all');
    plot(t,PrK); plot(t(K),PrK(K),'r*');
    xlabel('z (mm)','FontSize',12); ylabel(yLabel,'Interpreter','tex'); 
    title(titleString);
    xlim([min(t) max(t)])
    if  length(varargin{1})>0,
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0.1 3.2 2.4])
        set(gcf, 'PaperSize', [3.1 2.4]);
        print(h,'-dpdf',[outPath 'eps/' sprintf('%i',pathSuffix) sprintf('_%s_axial.pdf',titleString)])
        saveas(h,[outPath 'fig/' sprintf('%i',pathSuffix) sprintf('_%s_axial.fig',titleString)])
    end
    hold off; 
    pause;
end