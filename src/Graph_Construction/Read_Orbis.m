close all; clear all;

%% Load Orbis data

addpath ../Input/ORBIS/
% Read nodes
Nodes_Tab = readtable('Input/ORBIS/gorbit-nodes_label.csv');
label     = table2array(Nodes_Tab(:,6));
Nodes_Tab = [Nodes_Tab(1:end,1),Nodes_Tab(1:end,3),Nodes_Tab(1:end,4),Nodes_Tab(1:end,5)];
Nodes          = table2array(Nodes_Tab);
Nodes(:,end+1) = (1:size(Nodes,1))';
num_nodes      = size(Nodes,1);
[~,~,ib] = unique(label, 'stable');
label    = ib;

% Read edges
EdgeList_Tab = readtable('Input/ORBIS/gorbit-edges.csv');
EdgeList_Tab = [EdgeList_Tab(1:end,1),EdgeList_Tab(1:end,2),...
         EdgeList_Tab(1:end,3),EdgeList_Tab(1:end,4),...
         EdgeList_Tab(1:end,5)];
EdgeList = table2array(EdgeList_Tab(:,1:2));
Km        = table2array(EdgeList_Tab(:,3));
Days      = table2array(EdgeList_Tab(:,4));
Cost      = table2array(EdgeList_Tab(:,5));
num_edges = size(EdgeList,1);

EdgeList_re = zeros(num_edges,2);
% Renumber the edgelist
for i = 1:num_edges
    index_1 = EdgeList(i,1);
    index_2 = EdgeList(i,2);
    index_1_re = find(Nodes(:,1)==index_1);
    index_2_re = find(Nodes(:,1)==index_2);
    EdgeList_re(i,:) = [index_1_re,index_2_re];
end

Connectivity_Km   = accumarray(EdgeList_re, Km, [num_nodes, num_nodes]);
Connectivity_Days = accumarray(EdgeList_re, Days, [num_nodes, num_nodes]);
Connectivity_Cost = accumarray(EdgeList_re, Cost, [num_nodes, num_nodes]);

%% Visualization

% Visualize the sparsity pattern of the network
figure;
spy(Connectivity_Km,'k.')

% Visualize the location of settlements
marker_size   = 20;
font_size     = 18;

% Coordinates of all settlements
coords = [Nodes(:,4),Nodes(:,3)];

figure;
gx = geoaxes;
geoscatter(Nodes(:,4),Nodes(:,3),marker_size,'k','filled');
hold on;
gx.Basemap = 'colorterrain';
gx.Scalebar.Visible = 'off';
gx.TickDir = 'none';
gx.LatitudeAxis.TickLabels = {};
gx.LongitudeAxis.TickLabels = {};
gx.LatitudeLabel.String = '';
gx.LongitudeLabel.String = '';
grid off;
% geolimits([17.811456,-14.150391],[57, -14.150391]);
% exportgraphics(gx,'Rome_Locations.pdf','ContentType','vector')
print('Rome_Locations','-dpng','-r200')


% hold on to plot the graph on top
coords_curr = [Nodes(:,4),Nodes(:,3)];

[i,j]       = find(Connectivity_Days);
[ignore, p] = sort(max(i,j));
i = i(p);
j = j(p);

% Create a long, NaN-seperated list of line segments,
% rather than individual segments.

X = [ coords_curr(i,1) coords_curr(j,1) NaN*ones(size(i))]';
Y = [ coords_curr(i,2) coords_curr(j,2) NaN*ones(size(i))]';
X = X(:);
Y = Y(:);
hold on;
geoplot(X, Y,'Linewidth',2, 'color','black');
print('Rome_Graph','-dpng','-r200')

%% Build adjacency matrices
C        = (Connectivity_Km + Connectivity_Km')/2;
[W_Km]   = Get_Adj_from_Cost(C);
W        = (W_Km + W_Km')/2;
save('Orbis_Km_sym_three.mat','W','C','label','coords');
C        = (Connectivity_Days + Connectivity_Days')/2;
[W_Days] = Get_Adj_from_Cost(C);
W        = (W_Days + W_Days')/2;
save('Orbis_Days_sym_three.mat','W','C','label','coords');
C        = (Connectivity_Cost + Connectivity_Cost')/2;
[W_Cost] = Get_Adj_from_Cost(C);
W        = (W_Cost + W_Cost')/2;
save('Orbis_Cost_sym_three.mat','W','C','label','coords');

% Starting nodes
% Rome 50327, Constantinopolis 50129, Londinium 50235
Start_idx_Rome  = find(Nodes(:,1)==50327);
Start_idx_Const = find(Nodes(:,1)==50129);
Start_idx_Lond  = find(Nodes(:,1)==50235);