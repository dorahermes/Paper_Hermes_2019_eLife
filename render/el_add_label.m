function el_add_label(els, labels, varargin)
% EL_ADD_LABEL add 3-D labels of electrodes
% 
% Syntax:
%   el_add_label(els, labels, 'Name', value)
% 
% Input(s)
%   els             - [num] 3D position of the electrodes
%   labels          - [string] labels of electrodes
%   ax              - [obj] (opt) axes object (default = gca)
%   FontSize        - [num] (para) font size (default = 12)
%   Color           - [str|num] (para) color of text (default = 'k')
% 
%  Output(s)
% 
% See also .

% Copyright 2020 Richard J. Cui. Created: Sun 03/08/2020 10:49:56.385 PM
% $Revision: 0,1 $  $Date: Sun 03/08/2020 10:49:56.385 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(els, labels, varargin{:});
els = q.els;
labels = q.labels;
ax = q.ax;
if isempty(ax), ax = gca; end % if
font_size = q.FontSize;
color = q.Color;

% =========================================================================
% main
% =========================================================================
tf = ishold(ax);
if tf == false
    hold(ax, 'on')
end % if

num_els = size(els, 1);
for k = 1:num_els
    els_k = els(k, :);
    label_k = labels(k);
    text(els_k(1), els_k(2), els_k(3), label_k, 'FontSize', font_size,...
        'Color', color);
end % for

if tf == false
    hold(ax, 'off')
end % if

end % function

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults
default_ax = [];
default_fs = 12;
default_cl = 'k';

% parse rules
p = inputParser;
p.addRequired('els', @isnumeric);
p.addRequired('labels', @(x) ischar(x) || isstring(x));
p.addOptional('ax', default_ax, @(x) isempty(x) || isobject(x));
p.addParameter('FontSize', default_fs, @isnumeric);
p.addParameter('Color', default_cl, @(x) isnumeric(x) || isstr(x));

% parse and return results
p.parse(varargin{:});
q = p.Results;

end % funciton

% [EOF]