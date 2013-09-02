function DDS=check_data(DDS,type)
% internal DIWASP1.1 function
% checks data structures
%
%   DDS=check_data(DDS,type)
%       DDS:  the data structure
%       type: 1, Instrument data structure;
%             2, Spectral matrix structure;
%             3, Estimation parameters structure.

% Updated on 29/04/2013 by r.guedes to fix some errors that were not being
% detected (e.g., ID.datatype not being a cell array).


%--------------------------------------------------------------------------
% Defaults
%--------------------------------------------------------------------------
SM.xaxisdir=90;
EP.dres=180; EP.nfft=[]; EP.method='EMEP'; EP.iter=100;
error = '';

switch type
    %----------------------------------------------------------------------
    % Instrument data structure
    %----------------------------------------------------------------------
    case 1
        if ~isstruct(DDS)
            disp('DIWASP data_check: Instrument data type is not a structure');
        end
        nc=1;
        if isfield(DDS,'layout')
            [nr,nc]=size(DDS.layout);
            if nr<3
                if nr==2
                    DDS.layout(3,:)=0;
                else
                    error='layout';
                end
            end
        end
        
        if ~isfield(DDS,'datatypes') || ~iscell(DDS.datatypes) || min(size(DDS.datatypes))~=1 || max(size(DDS.datatypes))~=nc
            error='datatypes';
        else
            DDS.datatypes=DDS.datatypes(:)'; % make sure to be arow cell
        end
        
        if ~isfield(DDS,'depth') || length(DDS.depth)~=1
            error='depth';
        end
        
        if ~isfield(DDS,'fs') || length(DDS.fs)~=1
            error='fs';
        end
        
        if isfield(DDS,'data')
            if size(DDS.data,2)~=nc
                error='data';
            end
        else
            DDS.data=zeros(1,nc);
        end
        
        if ~isempty(error)
            disp(' ')
            disp(['Instrument data structure error: field [' error '] not specified correctly']);
            DDS=[];
            return;
        end
        
    %----------------------------------------------------------------------
    % Spectral matrix
    %----------------------------------------------------------------------
    case 2
        if ~isstruct(DDS)
            disp('DIWASP data_check: Spectral matrix data type is not a structure');
        end
        
        if isfield(DDS,'freqs') && min(size(DDS.freqs))==1
            nf=length(DDS.freqs);
        else
            error='freqs';
        end
        
        if isfield(DDS,'dirs') && min(size(DDS.dirs))==1
            nd=length(DDS.dirs);
        else
            error='dirs';
        end
        
        if isfield(DDS,'S')
            if size(DDS.S,1)~=nf || size(DDS.S,2)~=nd %|| isempty(DDS.S)
                error='S';
            end
        else
            DDS.S=[];
        end
        
        if isfield(DDS,'xaxisdir')
            if length(DDS.xaxisdir)~=1
                error='xaxisdir';
            end
        else
            DDS.xaxisdir=SM.xaxisdir;
        end
        
        if ~isfield(DDS,'dunit')
            DDS.dunit='cart';
        end
        
        if ~isfield(DDS,'funit')
            DDS.funit='hz';
        end
        
        if ~isempty(error)
            disp(' ')
            disp(['Spectral matrix structure error: field [' error '] not specified correctly']);
            DDS=[];
            return;
        end
        
    %----------------------------------------------------------------------
    % Estimation parameters
    %----------------------------------------------------------------------
    case 3
        if ~isstruct(DDS)
            disp('DIWASP data_check: Estimation parameter data type is not a structure');
        end
        
        if isfield(DDS,'dres')
            if length(DDS.dres)~=1
                error='dres';
            elseif DDS.dres<10
                DDS.dres=10;
                warning('dres is too small and has been set to 10')
            end
        else
            DDS.dres=EP.dres;
        end
        
        if isfield(DDS,'nfft')
            if length(DDS.nfft)~=1
                error='nfft';
            elseif DDS.nfft<64
                DDS.nfft=64;
                warning('nfft is too small and has been set to 64')
            end
        else
            DDS.nfft=EP.nfft;
        end
        
        if isfield(DDS,'iter')
            if length(DDS.iter)~=1
                error='iter';
            end
        else
            DDS.iter=EP.iter;
        end
        
        if isfield(DDS,'smooth')
            if ~strcmp(DDS.smooth,'OFF')
                DDS.smooth='ON';
            end
        else
            DDS.smooth='ON';
        end
        
        if isfield(DDS,'method')
            if ~(any(strcmp(DDS.method,{'DFTM','EMLM','IMLM','EMEP','BDM'})))
                error='method';
            end
        else
            DDS.method=EP.method;
        end
        
        if ~isempty(error)
            disp(' ')
            disp(['Estimation parameters structure error: field [' error '] not specified correctly']);
            DDS=[];
            return;
        end
        
    otherwise
        disp(' ')
        warning('DIWASP data_check: Data type unknown');
        DDS=[];
end

