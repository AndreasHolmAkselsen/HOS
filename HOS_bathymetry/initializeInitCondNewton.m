function zz = initializeInitCondNewton(fz,dfz,z_target,xi0,nIt)

% eta = 0*xi0;




% xi = xi0;
% zz = xi0+1i*.1;
zz = z_target/max(abs(dfz(xi0)));
xx = real(zz);yy=imag(zz);
for i = 1:nIt      
    f = fz(zz);
    df = dfz(zz);
    
%     zz = zz - (f-z_target).*df./abs(df).^2;
    
    xx = xx - real(f-z_target).*real(df)./abs(df).^2;
    yy = yy - imag(f-z_target).*imag(df)./abs(df).^2;
    zz = xx+1i*yy;
end

figure('color','w');hold on
plot(z_target,'k.');plot(zz,'.')

end
