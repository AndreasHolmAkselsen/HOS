function [h_fig, h_ax] = multi_axes(n,m,h_fig,margin,sep)
% multi_axes(n, m, h_fig, margin, sep), 
% margin = [xLeft, xRight, yBottom, yTop], 
% sep = [sep_x,sep_y]

if exist('h_fig','var') && ~isempty(h_fig)
    figure(h_fig);
else
    h_fig = figure;
    h_fig.Color = 'white';
    h_fig.Position = [0         724        1586        1111];
end

if ~exist('margin','var') || isempty(margin)
    margin = [.05,.04,.05,.05];
end
margin_xL = margin(1);
margin_xR = margin(2);
margin_yB = margin(3);
margin_yT = margin(4);


if exist('sep','var') && ~isempty(sep)
    sep_x = sep(1);
    sep_y = sep(2);
else  
    sep_x = .05;
    sep_y = .05;
end

width =  (1-margin_xL-margin_xR-(m-1)*sep_x)/m;
height = (1-margin_yB-margin_yT-(n-1)*sep_y)/n;

for i = n:-1:1
    for j = m:-1:1
        h_ax(n-i+1,j) = axes('Position', ...
            [margin_xL+(j-1)*(width +sep_x),...
             margin_yB+(i-1)*(height+sep_y),...
             width, height]);
         hold on
    end
end

% for i = 1:n,  for j = 1:m, axes(h_ax(i,j)); text(.5,.5,num2str([i,j])); end; end

end