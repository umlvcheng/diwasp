function trm=accx(ffreqs,dirs,wns,z,depth)

Kz=cosh(z*wns)./sinh(depth*wns);
%include a maximum cuttoff for the velocity response function
Kz(find(Kz<0.1))=0.1;
trm=-i*(ffreqs.*ffreqs.*Kz)*sin(dirs);