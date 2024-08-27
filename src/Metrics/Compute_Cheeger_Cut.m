function [hs] = Compute_Cheeger_Cut(Lap,indic,varargin)
% Compute the conductance (Cheeger Cut) between two partitions
% Input:
% Lap: Graph Laplacian
% indic: 1/0 binary indicator vector denoting partition membership
% variable:
%         min: consider the minimum in the denominator between the volume
%              of a partition and its complement
%       nomin: consider the volume of the partition indicated with 1's as the denominator 
% Output:
% hs: the conductance value
% 
% Adapted from https://github.com/estbautista/Label_Propagation_SSL/blob/master/SSL_methods/graph_tools/Cheeger_ratio.m

variable = 'min';
for i = 1 : length(varargin)
	variable = varargin{i};
end

% Degree matrix (diag Laplacian entries)
Deg = diag(diag(Lap));

% Volume of entire graph
volG = sum(diag(Deg));

% Volume of the partition under test
volS = indic'*Deg*indic;

% boundary (cut edges)
boundary  = indic'*Lap*indic;

% Cheeger ratio/Conductance
switch variable
    case 'nomin'
        hs = boundary/volS;
    case 'min'
        hs = boundary/min(volS, volG - volS);
    otherwise
end