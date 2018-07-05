% Function to read.FOLD files and convert into a format MERLIN2 can use.
% Created by Steve Grey 25/06/2018 - University of Bristol
% Email: steven.grey@bristol.ac.uk

function [Node,Panel] = FOLD2MERLIN(filename)

% Obtain the nodes and panel data
[~,~,vertices,~,faces,~,~] = FOLD_reader(filename);

% Convert from cell to array
Node = cell2mat(vertices.coords);

% Convert from structure to cell
Panel = faces.vertices;
