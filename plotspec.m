function plotspec(SM,ptype);

%DIWASP V1.4 function
%plotspec: plots the spectral matrix in 3D or polar form
%
%plotspec(SM,ptype)
%
%Inputs:
% SM   		A spectral matrix structure
% ptype		plot type:
%   1	3D surface plot
%   2	polar type plot 
%   3	3D surface plot (compass bearing angles direction from)
%   4	polar type plot (compass bearing angles direction from)
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

fig=figure;

SM=check_data(SM,2);if isempty(SM) return;end;
[SM,sfac]=spectobasis(SM);%Convert to basis matrix
dirs=SM.dirs;ffreqs=SM.freqs/(2*pi);S=2*pi^2*real(SM.S)/180;

%Convert directrions to nautical
if (ptype==3|ptype==4)
      if isfield(SM,'xaxisdir')
         xaxisdir=SM.xaxisdir;
      else
         xaxisdir=90;
      end
      dirs=dirs+pi+pi*(90-xaxisdir)/180;
end

%Surface plots
if(ptype==1|ptype==3)
    if ptype==3;dirs=mod(dirs,2*pi);end
    [dirs,order]=sort(180*dirs/pi);
    [ddir,df]=meshgrid(dirs,ffreqs);
    S=S(:,order);
    surf(df,ddir,real(S));
   shading interp;
   xlabel('frequency [Hz]');
   if(ptype==1)
      ylabel('direction [degrees]');
      axis([0 (max(ffreqs)) -180 180 0 (max(max(S)))]);

   else
      ylabel('direction [bearing]');
      axis([0 (max(ffreqs)) 0 360 0 (max(max(S)))]);
   end
   zlabel('m^2s / deg');
   
%Polar plots
elseif(ptype==2|ptype==4)
   h = polar([0 2*pi], [0 0.8*max(ffreqs)]);
   delete(h);
   
   [df,ddir]=meshgrid(ffreqs,dirs);
   
%uses the existing polar figure function and replaces numbering of angles for compass directions. Will probably be changed in future versions.
   if(ptype==4)
   set(0,'ShowHiddenHandles','on')
   chhs=get(gca,'Children');
   for i=1:size(chhs,1);
      obj=chhs(i);
      if strcmp(get(obj,'Type'),'text')
         num=str2num(get(obj,'String'));
         if~(isempty(num))
         if mod(num,30)==0
            num=90-num;
            num=(num<0)*360+num;
            set(obj,'String',num2str(num));
         end
         end
      end
   end
   set(0,'ShowHiddenHandles','off')
   end
         
   hold on;
    
	[px,py]=pol2cart(ddir,df);
	contour(px,py,real(S'),20);

   caxis([0.000 max(max(S))]);
   colorbar('vert');
   if(ptype==2)
      ylabel('direction [degrees] / frequency [Hz]');
   else
      ylabel('direction [bearing] / frequency [Hz]');
   end
	xlabel('m^2s / deg');
   hold off;
end

set(gca,'Color','none');


