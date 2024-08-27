% Worms Plot for the Erlangen Newsletter


% script to visualize the Artificial Gauss datasets for Section 3.
close all;


% load High_Moon_2.mat;
addpath wgPlot/wgPlot/

% Worms D3
% load worms_constD3_10NN.mat;
% load worms_constD3_all.mat;
% load worms_constD3_clusters.mat;

% 
load Gauss8_55_10NN.mat;
load Gauss8_sol.mat;

start_node = 439;
coords = data;
% coords   = worms_constD3_data;
clusters = best_cluster+1;


n        = size(clusters,1);

font_size = 74;

% PLOT PURE INPUT DATA
figure('Renderer', 'painters', 'Position', [0 0 10 7]*250);
% figure;
set(gcf,'color','w');

how_big = 30;
gray_is = [160 160 160]/255;
s = scatter(coords(:,1),coords(:,2),how_big,gray_is,'filled');
s.MarkerEdgeColor = 'k';
s.LineWidth       = 1; 
h = get(gca);
axis off;
% xlabel('Initial data','interpreter','latex','fontsize',font_size);

h.XAxis.Label.Visible='on';
tightfig;
saveas(gcf,'Gauus_Initial','pdf');



% PLOT P-SPEC CLUSTERS
% color_node_partitions = [0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250];
blue_is   = [30, 136, 229]/255;
green_is  = [0, 77, 64]/255;
yellow_is = [255, 193, 7]/255; 
red_is    = [220,20,60]/255;
purple_is = [172,79,198]/255; 
black_is  = [0,0,0]/255;
check = [0         0    1.0000];
color_node_partitions = [gray_is; green_is; check];


% SCATTTER THE POINTS

figure('Renderer', 'painters', 'Position', [0 0 10 7]*250);
% PLOT CUT EDGES
hold on;
color_edgecut = red_is;


[i,j] = find(W);
f = find(clusters(i) > clusters(j));
cut_edges = sparse(i(f),j(f),1,n,n);

[i,j] = find(cut_edges);
[ignore, p] = sort(max(i,j));
i = i(p);
j = j(p);


X = [ coords(i,1) coords(j,1) NaN*ones(size(i))]';
Y = [ coords(i,2) coords(j,2) NaN*ones(size(i))]';
X = X(:);
Y = Y(:);
plot (X, Y, '-', 'Linewidth',7, 'color',color_edgecut);

hold on;
% PLOT EDGES WITHIN THE CLUSTERS
parts = setfilter(clusters);
ncolor = 1;
color_edges_within = gray_is;
for partnumber = 2
    c = color_node_partitions(partnumber,:);

    part_curr = find(clusters == partnumber);
    
    % Input     
    W_curr      = W(part_curr,part_curr);
    coords_curr = coords(part_curr,:);
    
    
    [i,j]       = find(W_curr);
    [ignore, p] = sort(max(i,j));
    i = i(p);
    j = j(p);

    % Create a long, NaN-seperated list of line segments,
    % rather than individual segments.

    X = [ coords_curr(i,1) coords_curr(j,1) NaN*ones(size(i))]';
    Y = [ coords_curr(i,2) coords_curr(j,2) NaN*ones(size(i))]';
    X = X(:);
    Y = Y(:);
        
    h =  plot(X, Y,'-','Linewidth',5, 'color',color_edges_within);
       
end



hold on;
set(gcf,'color','w');
how_big = 500;
label = label +1;
map_color = colormap(color_node_partitions);
s = scatter(coords(:,1),coords(:,2),how_big,clusters,'filled');
s.MarkerEdgeColor = 'k';
s.LineWidth       = 1;     
hold on;
scatter(coords(start_node,1),coords(start_node,2),3000,yellow_is,'diamond','filled');

k = get(gca);
axis off;
% xlabel('$p$-spectral clusters','interpreter','latex','fontsize',font_size);
k.XAxis.Label.Visible='on';


tightfig;
saveas(gcf,'Gauss_Clusters_re','pdf');



