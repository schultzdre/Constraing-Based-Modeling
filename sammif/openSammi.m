% USAGE:
% openSAMMI(htmlName)
% 
% OPTIONAL INPUT:
%   htmlName: Name of the html file previously written using the sammi
%   function. Defaults to 'html_load'.
% 
% OUTPUT:
%   No MATLAB output, opens the visualization in a new browser tab.
function openSammi(htmlName)
    %Set default
    if nargin < 1
        htmlName = 'index_load.html';
    end
    %Get COBRA directory
    global CBTDIR;
    %Open file
    web([CBTDIR '/external/visualization/sammif/' htmlName],'-browser')
end