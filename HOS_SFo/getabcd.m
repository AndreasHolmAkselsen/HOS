function abcd=getabcd(nmax)

abcd=cell(nmax+1,nmax+1);
abcd{1,1}=[0 0 0 0 1];
for n=1:nmax
    for p=0:n       
        if p==0
            abcd{n+1,p+1}=zeros(2*size(abcd{n,p+1},1),5);
            for ic=1:size(abcd{n,p+1},1)
                abcd{n+1,p+1}(2*(ic-1)+1,:)=abcd{n,p+1}(ic,:);
                abcd{n+1,p+1}(2*(ic-1)+1,2)=abcd{n+1,p+1}(2*(ic-1)+1,2)+1;
                abcd{n+1,p+1}(2*(ic-1)+2,:)=abcd{n,p+1}(ic,:);
                abcd{n+1,p+1}(2*(ic-1)+2,4)=abcd{n+1,p+1}(2*(ic-1)+2,4)+1;
            end
        else
            abcd{n+1,p+1}=zeros(2*size(abcd{n,p},1),5);
            for ic=1:size(abcd{n,p},1)
                abcd{n+1,p+1}(2*(ic-1)+1,:)=abcd{n,p}(ic,:);
                abcd{n+1,p+1}(2*(ic-1)+1,1)=abcd{n+1,p+1}(2*(ic-1)+1,1)+1;
                abcd{n+1,p+1}(2*(ic-1)+2,:)=abcd{n,p}(ic,:);
                abcd{n+1,p+1}(2*(ic-1)+2,3)=abcd{n+1,p+1}(2*(ic-1)+2,3)+1;
            end
        end
        [temp,ia,ic]=unique(abcd{n+1,p+1}(:,1:4),'rows');
        for i0=1:length(ia)
            temp(i0,5)=sum(abcd{n+1,p+1}(ic==i0,5));
        end
        abcd{n+1,p+1}=temp;
    end
end

end