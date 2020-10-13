function toPrint(handle, rmTitleFlag, fontSz, options)
%
% function toPrint(handle, rmTitleFlag, fontSz)
%
% Given a figure indicated by handle
% remove its title and change font size
%
% Defaults: handle=gcf; rmTitleFlag=1; fontSz=16
%

% 14/2/01, 11.11.2010 (opt), J. Gaspar

if nargin<4
    options= [];
end
if nargin<3 || isempty(fontSz)
   fontSz=16;
end
if nargin<2 || isempty(rmTitleFlag)
   rmTitleFlag=1;
end
if nargin<1 || isempty(handle)
   handle= gcf;
end

figure(handle);
if rmTitleFlag,
   title('');
else
    h= get(gca, 'Title'); set(h,'FontName', 'Times New Roman'); set(h,'FontSize',fontSz);
end

% enlarge the fonts:
%
set(gca,'FontName', 'Times New Roman'); set(gca,'FontSize',fontSz);
h= get(gca,'XLabel'); set(h,'FontName', 'Times New Roman'); set(h,'FontSize',fontSz);
h= get(gca,'YLabel'); set(h,'FontName', 'Times New Roman'); set(h,'FontSize',fontSz);
h= get(gca,'ZLabel'); set(h,'FontName', 'Times New Roman'); set(h,'FontSize',fontSz);
h= get(gca,'Legend'); set(h,'FontName', 'Times New Roman'); set(h,'FontSize',fontSz);

% Obtain objects in legends 
[~,icons] = legend();
for iIcon = 1:length(icons)
    if strcmpi(icons(iIcon).Type,'text')
        icons(iIcon).FontName = 'Times New Roman';
        icons(iIcon).FontSize = fontSz;
    elseif strcmpi(icons(iIcon).Type,'Line')
        icons(iIcon).LineWidth  = 1;
        icons(iIcon).MarkerSize = fontSz;
    end
end

if ~isfield(options, 'norepos')
    % Get current axes
    currentAxis   = gca;
    outerPosition = currentAxis.OuterPosition;
    tightInset    = currentAxis.TightInset; 
    leftPoint     = outerPosition(1) + tightInset(1);
    bottomPoint   = outerPosition(2) + tightInset(2);
    axisWidth     = outerPosition(3) - tightInset(1) - tightInset(3);
    axisHeight    = outerPosition(4) - tightInset(2) - tightInset(4);
    currentAxis.Position = [leftPoint bottomPoint axisWidth axisHeight];
end
