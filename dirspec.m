function [SMout,EP]=dirspec(ID,SM,EP,varargin)

%DIWASP V1.4 function
%dirspec: main spectral estimation routine
%
%[SMout,EPout]=dirspec(ID,SM,EP,{options})
%
%Outputs:
%SMout	    A spectral matrix structure containing the results
%Epout		The estimation parameters structure with the values actually used for the computation including any default settings.
%
%Inputs:
%ID			An instrument data structure containing the measured data
%SM   		A spectral matrix structure; data in field SM.S is ignored.
%EP		    The estimation parameters structure. For default values enter EP as []
%[options]  options entered as cell array with parameter/value pairs: e.g.{'MESSAGE',1,'PLOTTYPE',2};
%                Available options with default values:
%                    'MESSAGE',1,    Level of screen display: 0,1,2 (increasing output)
%                    'PLOTTYPE',1,   Plot type: 0 none, 1 3d surface, 2 polar type plot, 3 3d surface(compass angles), 4 polar plot(compass angles)
%                    'FILEOUT',''  	 Filename for output file: empty string means no file output
%                
%Input structures ID and SM are required. Either [EP] or [options] can be included but must be in order if both are included.
%"help data_structures" for information on the DIWASP data structures

%All of the implemented calculation algorithms are as described by:
%Hashimoto,N. 1997 "Analysis of the directional wave spectrum from field data" 
%In: Advances in Coastal Engineering Vol.3. Ed:Liu,P.L-F. Pub:World Scientific,Singapore
%
%

%Copyright (C) 2002 Coastal Oceanography Group, CWR, UWA, Perth

Options=struct('MESSAGE',1,'PLOTTYPE',1,'FILEOUT','');

if nargin<3
   error('All inputs other than OPTIONS required');
elseif nargin>=4
   nopts=length(varargin{1});
else
   nopts=[];
end

ID=check_data(ID,1);if isempty(ID) return;end;
SM=check_data(SM,2);if isempty(SM) return;end;
EP=check_data(EP,3);if isempty(EP) return;end;

if ~isempty(nopts)
if(rem(nopts,2)~=0)
   	warning('Options must be in Name/Value pairs - setting to defaults');
else
   	for i=1:(nopts/2)
      	arg=varargin{1}{(2*i)};
      	field=varargin{1}{(2*i-1)};
      	Options=setfield(Options,field,arg);
    end;
end   
end

ptype=Options.PLOTTYPE;displ=Options.MESSAGE;

disp(' ');disp('calculating.....');disp(' ');disp('cross power spectra');

data=detrend(ID.data);
[ndat,szd]=size(ID.data);

%get resolution of FFT - if not specified, calculate a sensible value
if isempty(EP.nfft)
    nfft=2^(8+round(log2(ID.fs)));
    EP.nfft=nfft;
else
    nfft=EP.nfft;
end
if nfft>ndat;error(['Data length of ' num2str(ndat) ' too small']);end
    
%calculate the cross-power spectra
for m=1:szd
for n=1:szd
   [xpstmp,Ftmp]=diwasp_csd(data(:,m),data(:,n),nfft,ID.fs);
   xps(m,n,:)=xpstmp(2:(nfft/2)+1);
 end
end

F=Ftmp(2:(nfft/2)+1);nf=nfft/2;

%calculate wavenumbers
disp('wavenumbers')
wns=wavenumber(2*pi*F,ID.depth*ones(size(F)));
pidirs=[-pi:(2*pi/EP.dres):pi-(2*pi/EP.dres)];

%calculate transfer parameters
disp('transfer parameters');
disp(' ');
for m=1:szd
	trm(m,:,:)=feval(ID.datatypes{m},2*pi*F,pidirs,wns,ID.layout(3,m),ID.depth);
 for n=1:szd
    kx(m,n,:,:)=wns*((ID.layout(1,n)-ID.layout(1,m))*cos(pidirs)+(ID.layout(2,n)-ID.layout(2,m))*sin(pidirs));
end
end

for m=1:szd
   tfn(:,:)=trm(m,:,:);
   Sxps(1:nf)=xps(m,m,:);
   Ss(m,:)=Sxps./(max(tfn').*conj(max(tfn')));
end

ffs=(F>=min(SM.freqs) & F<=max(SM.freqs));
SM1.freqs=F(ffs);
SM1.funit='Hz';
SM1.dirs=pidirs;
SM1.dunit='rad';

% call appropriate estimation function
disp(['directional spectra using' blanks(1) EP.method ' method']);disp(' ');
SM1.S=feval(EP.method,xps(:,:,ffs),trm(:,ffs,:),kx(:,:,ffs,:),Ss(:,ffs),pidirs,EP.iter,displ);
SM1.S(isnan(SM1.S) | SM1.S<0)=0.;

%Interpolate onto user specified matrix
disp(' ');disp('interpolating onto specified matrix...');disp(' ');
SMout=interpspec(SM1,SM);

%smooth spectrum
if(strcmp(EP.smooth,'ON'))
    disp(' ');disp('smoothing spectrum...');disp(' ');
    SMout=smoothspec(SMout,[1 0.5 0.25;1 0.5 0.25]);
end

infospec(SMout);

%write out spectrum matrix in DIWASP format
filename=Options.FILEOUT;
if(size(filename,2)>0)
   disp('writing out spectrum matrix to file');
   writespec(SMout,filename);
end

%plot spectrum
if(ptype>0)
    disp('finished...plotting spectrum');
    plotspec(SMout,ptype);
    T=['Directional spectrum estimated using ' blanks(1) EP.method ' method'];title(T);
end