function [prf_size] = xy2prfsize(xys,v_area)

% function xy2prfsize calculates estimate the population receptive
% field (pRF) size from the pRF x and y positions in degrees. It uses the
% relation between pRF eccentricity and pRF size from Kay et al., 2013
% Journal of Neurophysiology, Figure 7A.
%
% input: 
%   xys = [x y optional]; x and y position of pRF in degrees of visual angle 
%   v_area = [1]; use 2 for V1, 3 for V2, 4 for V3 and 6 for V4 as in freesurfer label
%
% output:
%   prf_size: prf size in degrees of visual angle.
%
% D Hermes 2018


% slope from eccentricity to size
if v_area(1) == 2 %V1
    this_slope = 0.8/5;    
elseif v_area(1) == 3 %V2
    this_slope = 1/5;
elseif v_area(1) == 4 % V3
    this_slope = 1.6/5;
elseif v_area(1) == 6 %V4
    this_slope = 3/5;
end

% x, y to eccentricity:
xy_ecc = sqrt(xys(1).^2 + xys(2).^2); % diagonal = sqrt(x.^2 + y.^2)

% eccentricity to size:
prf_size = this_slope * xy_ecc;


