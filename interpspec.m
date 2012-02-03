function SMout=interpspec(SMin,SMout)

%DIWASP V1.4 function
%interpspec: interpolates between spectral matrix bases
%
%SMout=interpspec(SMin,SMout)
%
%Outputs:
%SMout		Output spectral matrix structure with interpolated power density
%
%Inputs:
%SMin   	A spectral matrix structure containing the original spectra
%SMout      A spectral matrix defining the new spectral matrix
%
%SMout only needs to have the frequency and directional axes defined -
%spectral density ignored
%
%"help data_structures" for information on the DIWASP data structures

Hs1=hsig(SMin)

[SMin,facin]=spectobasis(SMin);
[SMtmp,facout]=spectobasis(SMout);

s1=SMin.freqs(:)*sin(SMin.dirs(:)');c1=SMin.freqs(:)*cos(SMin.dirs(:)');
s2=SMtmp.freqs(:)*sin(SMtmp.dirs(:)');c2=SMtmp.freqs(:)*cos(SMtmp.dirs(:)');

Stmp=griddata(s1,c1,SMin.S,s2,c2);
Stmp(isnan(Stmp))=0;
SMout.S=Stmp/facout;

%check Hsig of mapped spectrum and check sufficiently close to original
Hs2=hsig(SMout);
if (Hs2-Hs1)/Hs1 >0.02
    warning('User defined grid may be too coarse; try increasing resolution of ''SM.freqs'' or ''SM.dirs''');
end








