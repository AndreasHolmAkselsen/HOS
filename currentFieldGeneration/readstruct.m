function Res=readstruct(FileName)

fid=fopen(FileName,'r');
Content=textscan(fid,'%s','delimiter','\n');
for iL=1:length(Content{1})
    eval(['Res.' Content{1}{iL}]);
end
fclose(fid);

end