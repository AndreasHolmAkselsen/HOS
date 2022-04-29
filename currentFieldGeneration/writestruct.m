function writestruct(fid, varname)
% writestruct(fid, varname)
%
% Write the contents of a structure to a file which
% has been opened with file id fid.
% 

% Test if fid is valid:
try
    a=ftell(fid);
catch
    fprintf(1, '%d is not a valid file id');
    return;
end

% Get list of elements in object
names = evalin('base',['fieldnames(', varname, ')']);

% Check if structure is an array of structures (only 1-D array supported)
n = evalin('base',['length(',varname,')']);
for s=1:n
    for i=1:length(names)
        fld = names{i};
        if (n>1)
            subvarname = [varname, '(', num2str(s), ').', fld];
        else
            subvarname = [varname,'.', fld];
        end
        
        if (evalin('base',['isstruct(',subvarname,')']))
            writestruct(fid, subvarname);
        elseif (evalin('base',['iscell(',subvarname,')']))
            for iC=1:evalin('base',['length(',subvarname,')'])
                if (evalin('base',['isstruct(',subvarname,'{' num2str(iC) '})']))
                    writestruct(fid, [subvarname '{' num2str(iC) '}']);
                else
                    writearray(fid, [subvarname '{' num2str(iC) '}']);
                    
                end
            end
        else
            writearray(fid, subvarname);
        end
    end
end





function writearray(fid, varname)
% writearray(fid, varname)
%
%  Output contents of a non-structure variable as an assignment statement
%  that can be executed later to re-create the variable.
%
%  Large variables (arbitrarily, the cutoff is at 200 elements) are saved
%  to a MAT file rather than written out inline.
%
%   R. Poe, 10/25/01
var = evalin('base',varname);
% Trivial cases: empty and scalar
if (isempty(var))
    fprintf(fid, '%s = [];\n', varname);
    return;
end
dims = size(var);
nd = length(dims);

if (prod(dims)==1)
    fprintf(fid, '%s = ', varname);
    if (ischar(var))
        fprintf(fid, '''%s'';\n', var);
    else
        fprintf(fid, '%g;\n', var);
    end
    return
end

% Handle the case of a large array
if (prod(dims)>200)
    fprintf(fid,'%% Variable contains a large amount of data. Read in from file\n');
    % Generate some random integers for workspace and var names
    filenm = ['data',sprintf('%04d',fix(9999*rand(1)))];
    while (exist([filenm,'.mat'], 'file'))
        filenm = ['data',sprintf('%04d',fix(9999*rand(1)))];
    end
    evalin('base',[filenm, ' = ',varname,';']);
    evalin('base',['save ',filenm,' ',filenm]);
    fprintf(fid, 'load %s\n%s = %s;\n', filenm, varname, filenm);
    return;
end

ii = ones(nd-2,1);
carry = 0;
while (~carry)
    if (nd>2)
        varref = [varname,'(:,:',sprintf(',%d', ii),')'];
    else
        varref = varname;
    end
    m = dims(1); n=dims(2);
    subvar = evalin('base',varref);
    cflag = ischar(subvar);
    fprintf(fid, '%s = ', varref);
    if (~cflag || m>1) fprintf(fid, '['); end
    for i=1:m
        if (cflag)
            fprintf(fid, '''%s''', subvar(i,:));
        else
            fprintf(fid, '%g ',subvar(i,:));
        end
        if (i<m) fprintf(fid,';'); end
    end
    if (~cflag || m>1) fprintf(fid, ']'); end
    fprintf(fid, ';\n');
    carry = 1;
    for k=1:nd-2
        ii(k) = ii(k)+carry;
        carry = (ii(k)>dims(k+2));
        if (carry) ii(k)=1; end
    end
end
