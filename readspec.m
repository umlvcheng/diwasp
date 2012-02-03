function [SM]=readspec(filename,funit,dunit)

%DIWASP V1.4 function
%readspec: reads in spectrum file in DIWASP format and displays surface plot
%
%[SM]=readspec(filename)
%
%Outputs:
%SM   		A spectral matrix structure containing the file data
%
%Inputs:
%filename	filename for the file in DIWASP format including file extension
%
%"help data_structures" for information on the DIWASP data structures

%Copyright (C) 2002 Coastal Oceanography Group, CWR, UWA, Perth

if nargin<2;funit='hz';end
if nargin<3;dunit='cart';end

datain=load(filename); 

SM.xaxisdir=datain(1);
nfreq=datain(2);
ndirs=datain(3);
SM.freqs=datain(4:nfreq+3);SM.funit=funit;
SM.dirs=datain(nfreq+4:nfreq+3+ndirs);SM.dunit=dunit;
headercheck=datain(nfreq+ndirs+4);
if headercheck~=999
   error('corrupt file header');
end
mat=datain(nfreq+ndirs+5:nfreq+ndirs+4+(nfreq*ndirs));

S=reshape(mat,ndirs,nfreq);
SM.S=S';

plotspec(SM,1);
