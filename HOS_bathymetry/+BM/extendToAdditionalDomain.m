function [fAdd,zAdd]=extendToAdditionalDomain(f,z,hAddRatio)

%Extended z-vector
h=abs(min(z));
dz=z(end)-z(end-1);
Nz=length(z)-1;
zAdd=[z,dz*(1:(hAddRatio*Nz))];

%Non-dimensional vertical coordinate 
hAdd=dz*(hAddRatio*Nz);
zm=(hAdd-h)/(2*h);
zHat=zAdd(Nz+2:hAddRatio*Nz)/(h*zm)-1;

%Find polynomial coefficients for C2 match, cf Bonnefoy et al. (2006), App. Oc. Res. 28.
f0=f(:,end);
df0=(f(:,end)-f(:,end-1))/dz;
C1=(17*f0+3*df0*(h*zm))/18;
C3=-(79*f0+15*df0*(h*zm))/18;
C5=(22*f0+6*df0*(h*zm))/9;

%Extend function to additional domain over surface elevation
fAdd=zeros(size(f,1),length(zAdd));
fAdd(:,1:Nz+1)=f;
fAdd(:,hAddRatio*Nz+1:end)=-fliplr(f);
fAdd(:,Nz+2:hAddRatio*Nz)=C1*zHat+C3*zHat.^3+C5*zHat.^5;

end