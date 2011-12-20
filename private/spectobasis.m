function [SM,Sfac]=spectobasis(SM)
%Converts any spectral matrix to rad/s and cartesian radians

Sfac=1.0;
if isfield(SM,'funit') && strcmpi(SM.funit,'hz')
    SM.freqs=2*pi*SM.freqs;
    Sfac=Sfac/(2*pi);
end

r2d=pi/180;
if isfield(SM,'dunit')
    if strcmpi(SM.dunit(1:3),'car')
        SM.dirs=r2d*SM.dirs;
        Sfac=Sfac/r2d;
    elseif strcmpi(SM.dunit(1:3),'nau')
        if isfield(SM,'xaxisdir');
            SM.dirs=SM.dirs+(90-SM.xaxisdir);
        end
        SM.dirs=r2d*(270-SM.dirs);
        Sfac=Sfac/r2d;
    end
end
SM.S=SM.S*Sfac;
