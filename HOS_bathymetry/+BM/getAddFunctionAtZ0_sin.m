function f=getAddFunctionAtZ0_sin(A,k,h,chiLargeX,chiLargeK,indK)

indKLargeX=indK&(k<=k(size(chiLargeX,1)));
indKChiLargeX=indK(1:size(chiLargeX,1));
indKLargeK=indK&(k>k(size(chiLargeX,1)));
indKChiLargeK=indK(size(chiLargeX,1)+1:end);
fLargeX=sum(A(indKLargeX).'.*(chiLargeX(indKChiLargeX,:).*sin(k(indKLargeX).'*h)),1)/(length(A)-1);
fLargeK=[sum(A(indKLargeK).'.*(chiLargeK(indKChiLargeK,:).*sin(k(indKLargeK).'*h)),1)/(length(A)-1),zeros(1,size(chiLargeX,2)-size(chiLargeK,2))];
f=fLargeX+fLargeK;

% f=sum(A(indKeep).'.*(chi(indKeep,:).*cos(k(indKeep).'*h)))/(length(A)-1);

end

