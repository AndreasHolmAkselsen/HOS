function [X,z]=initializeFlapMotion_SingleFlap(theta,h,hF,Nz)

%Initialize z vector
dz=h/Nz;
z=dz*(0:Nz)-h;

%Compute flap position matrix
X=zeros(length(theta),Nz);
X(:,z>-hF)=tan(theta/180*pi)*(z(z>-hF)+hF);

end
