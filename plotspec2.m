function [H,C] = plotspec2(SM,ptype,maxf,nf,islab,isclb)

%MODIFIED from DIWASP V1.4 function
%plotspec2: plots the spectral matrix in 3D or polar form
%
%plotspec2(SM,ptype,maxf,nfig,islab,isclb)
%
%Inputs:
%    SM     A spectral matrix structure
% ptype     plot type:
%     1     3D surface plot
%     2     polar type plot 
%     3     3D surface plot (compass bearing angles direction from)
%     4     polar type plot (compass bearing angles direction from)
%  maxf     Maximum frequency to be displayed in polar plot
%  nfig     1 to call a new figure, zero not to call 
% islab     Define between showing axes labels (1, default) or not (0)
% isclb     Define between showing colorbar (1, default) or not (0)
%
%Outputs:
%   H: Handles for the plot.
%   C: Handles for the colorbar, in case ptype is 3 or 4 (polar).
%
%The 3D surface plot type is a MATLAB surface plot with SM.freqs on the x axis, SM.dirs on the y axis and the spectral density, SM.S as the z value. 
%The polar type plot is a MATLAB polar plot with the direction showing values in SM.dirs, the radius showing values in SM.freqs 
%and contours representing the spectral density, SM.S. An example of the polar type plot is shown on the front cover of the manual.
%For plot types 1 and 2, the direction is the direction of propagation relative to the Cartesian axis. 
%For options 3 and 4 the direction is coming from as a true compass bearing (this has changed from previous versions). 
%Directions are corrected internally from the SM.xaxisdir and SM.dunit
%fields that define the orientation of the axes and directional units in the spectral matrix. 
%
%"help data_structures" for information on the DIWASP data structures

%Copyright (C) 2002 Coastal Oceanography Group, CWR, UWA, Perth

% Version 2 modified from original by Rafael to avoid plotting a new
% figure, as well as correcting gap on polar contourf plots, and outputting
% the handles for the contour object and colorbar.


if ispc
    fn='Arial';
elseif isunix
    fn = 'Liberation Sans';
end

if nargin<3 || nf==1
    figure
end
if ~exist('isclb','var')
    clb = 1;
end
if ~exist('islab','var')
    clb = 1;
end

SM = check_data(SM,2);
if isempty(SM) 
    return
end

[SM,sfac] = spectobasis(SM);    % Convert to basis matrix

  dirs = SM.dirs;
ffreqs = SM.freqs/(2*pi);
     S = 2*pi^2*real(SM.S)/180;
     
%==========================================================================
% RAFAEL UPDATE AUGUST 2012 - Close matrix to avoid gap on circunference
S = [S,S(:,1)];
dirs = [dirs;dirs(1)];
%==========================================================================

%Convert directrions to nautical
if ptype==3 || ptype==4
      if isfield(SM,'xaxisdir')
         xaxisdir = SM.xaxisdir;
      else
         xaxisdir = 90;
      end
      dirs = dirs+pi+pi*(90-xaxisdir)/180;
end

%Surface plots
if ptype==1 || ptype==3
    if ptype == 3;
        dirs = mod(dirs,2*pi);
    end
    [dirs,order] = sort(180*dirs/pi);
       [ddir,df] = meshgrid(dirs,ffreqs);
               S = S(:,order);
    H = surf(df,ddir,real(S));
    C = 0; % For not crashing when trying to output colorbar handles
    shading interp;
    xlabel('frequency [Hz]','FontName',fn);
    if ptype == 1
        ylabel('direction [degrees]','FontName',fn);
        axis([0 (max(ffreqs)) -180 180 0 (max(max(S)))]);
    else
        ylabel('direction [bearing]','FontName',fn);
        axis([0 (max(ffreqs)) 0 360 0 (max(max(S)))]);
    end
    zlabel('m^2s / deg','FontName',fn);

    %Polar plots
elseif ptype==2 || ptype==4
%     h = polar([0 2*pi], [0 0.8*max(ffreqs)]);
    h = polar([0 2*pi], [0 maxf]);
    delete(h);

    [df,ddir] = meshgrid(ffreqs,dirs);

    % Uses the existing polar figure function and replaces numbering of 
    % angles for compass directions. Will probably be changed in future versions.
    if ptype == 4
        set(0,'ShowHiddenHandles','on')
        chhs = get(gca,'Children');
        for i = 1:size(chhs,1);
            obj = chhs(i);
            if strcmp(get(obj,'Type'),'text')
                num = str2num(get(obj,'String'));
                if ~isempty(num)
                    if mod(num,30) == 0
                        num = 90-num;
                        num = (num<0)*360+num;
                        set(obj,'String',num2str(num));
                    end
                end
            end
        end
        set(0,'ShowHiddenHandles','off')
    end

    hold on;
    
    [px,py] = pol2cart(ddir,df);
    H = contour(px,py,real(S'),20);

    caxis([0.000 max(max(S))]);
    if isclb
        C = colorbar('vert');
    else
        C = 0;
    end
    
    if islab
        if ptype == 2
            ylabel('direction [degrees] / frequency [Hz]','FontName',fn);
        else
            ylabel('direction [bearing] / frequency [Hz]','FontName',fn);
        end
        xlabel('m^2s / deg','FontName',fn);
    end
    hold off
end

set(gca,'Color','none');


