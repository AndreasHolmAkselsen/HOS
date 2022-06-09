function flapMotion=computeFlapMotionDerivatives(X,dt,dz)

flapMotion=zeros(size(X,1),size(X,2),5);
flapMotion(:,:,1)=X;

% [~,flapMotion(:,:,2)]=gradient(X,dt);
% flapMotion(:,:,4)=gradient(X,dz);
[flapMotion(:,:,4),flapMotion(:,:,2)] = gradient(X,dz,dt);
[~,flapMotion(:,:,3)]=gradient(flapMotion(:,:,2),dt);
[~,flapMotion(:,:,5)]=gradient(flapMotion(:,:,4),dt);

end