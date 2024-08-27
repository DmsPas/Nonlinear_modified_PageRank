function Visualize_Orbis_Results(W,labels_full,labels_dist_all,coords,best_cluster_all)
% Visualize the classification results on the Orbis dataset with NPR clustering
% Input
% W               : adjacency matrix of the graph
% labels_full     : labels of the nodes
% labels_dist_all : labels of the nodes based on the shortest path distance from the seed node
% coords          : coordinates of the nodes
% best_cluster_all: cluster indicator vector

% Define colors
blue_is   = [30, 136, 229]/255;
green_is  = [0, 77, 64]/255;
yellow_is = [255, 193, 7]/255; 
red_is    = [220,20,60]/255;
purple_is = [172,79,198]/255;

%% Visualize all locations (Roman settlements) in the Orbis dataset
marker_size   = 50;
figure;
gx = geoaxes;
geoscatter(coords(:,1),coords(:,2),marker_size,'k','filled');
hold on;
gx.Basemap = 'colorterrain';
gx.Scalebar.Visible = 'off';
gx.TickDir = 'none';
gx.LatitudeAxis.TickLabels = {};
gx.LongitudeAxis.TickLabels = {};
gx.LatitudeLabel.String = '';
gx.LongitudeLabel.String = '';
grid off;
set(get(gca, 'Title'), 'String', 'Locations of all settlements', 'Interpreter', 'Latex','Fontsize',24);
print('NPR_Rome_Locations','-dpng','-r200')


%% Starting node is 126-Constantinoupolis

cluster_curr       = best_cluster_all(:,1);
labels_curr        = labels_dist_all(:,1);
adder              = labels_curr + cluster_curr ;
TP                 = find(adder==2);
coords_TP          = coords(TP,:);
subtr              = labels_curr - cluster_curr;
FP                 = find(subtr == -1);
coords_FP          = coords(FP,:);
FN                 = find(subtr == 1);
coords_FN          = coords(FN,:);
start_node         = 126;
marker_size        = 150;

figure;
gx = geoaxes;
geoscatter(coords_TP(:,1),coords_TP(:,2),marker_size,'black','filled');
hold on;
geoscatter(coords_FP(:,1),coords_FP(:,2),marker_size,blue_is,'filled');
hold on;
geoscatter(coords_FN(:,1),coords_FN(:,2),marker_size,red_is,'filled');
hold on;
geoscatter(coords(start_node,1),coords(start_node,2),marker_size+250,yellow_is,'diamond','filled');
legend('True Positive','False Positive','False Negative','Seed Node', 'Interpreter', 'Latex','Fontsize',16, 'Location','southwest');
gx.Basemap = 'colorterrain';
gx.Scalebar.Visible = 'off';
gx.TickDir = 'none';
gx.LatitudeAxis.TickLabels = {};
gx.LongitudeAxis.TickLabels = {};
gx.LatitudeLabel.String = '';
gx.LongitudeLabel.String = '';
grid off;
set(get(gca, 'Title'), 'String', '\textbf{Seed  node: Constantinopolis}', 'Interpreter', 'Latex','Fontsize',24);
geolimits([34.533712,48.144098],[21.740038, 42.960080]);
print('NPR_Orbis_Const','-dpng','-r200')

%% Starting node is 232-Londinium

cluster_curr       = best_cluster_all(:,2);
labels_curr        = labels_dist_all(:,2);
adder              = labels_curr + cluster_curr ;
TP                 = find(adder==2);
coords_TP          = coords(TP,:);
subtr              = labels_curr - cluster_curr;
FP                 = find(subtr == -1);
coords_FP          = coords(FP,:);
FN                 = find(subtr == 1);
coords_FN          = coords(FN,:);
start_node         = 232;
marker_size        = 150;

figure;
gx = geoaxes;
geoscatter(coords_TP(:,1),coords_TP(:,2),marker_size,'black','filled');
hold on;
geoscatter(coords_FP(:,1),coords_FP(:,2),marker_size,blue_is,'filled');
hold on;
geoscatter(coords_FN(:,1),coords_FN(:,2),marker_size,red_is,'filled');
hold on;
geoscatter(coords(start_node,1),coords(start_node,2),marker_size+250,yellow_is,'diamond','filled');
legend('True Positive','False Positive','False Negative','Seed Node', 'Interpreter', 'Latex','Fontsize',16, 'Location','southwest');
gx.Basemap = 'colorterrain';
gx.Scalebar.Visible = 'off';
gx.TickDir = 'none';
gx.LatitudeAxis.TickLabels = {};
gx.LongitudeAxis.TickLabels = {};
gx.LatitudeLabel.String = '';
gx.LongitudeLabel.String = '';
grid off;
set(get(gca, 'Title'), 'String', '\textbf{Seed  node: Londinium}', 'Interpreter', 'Latex','Fontsize',24);
geolimits([42.244785,56.182254],[-7.908365, 5.958885]);
print('NPR_Orbis_Lond','-dpng','-r200')


%% Starting node is 324-Rome

cluster_curr       = best_cluster_all(:,3);
labels_curr        = labels_dist_all(:,3);
adder              = labels_curr + cluster_curr ;
TP                 = find(adder==2);
coords_TP          = coords(TP,:);
subtr              = labels_curr - cluster_curr;
FP                 = find(subtr == -1);
coords_FP          = coords(FP,:);
FN                 = find(subtr == 1);
coords_FN          = coords(FN,:);
start_node         = 324;
marker_size        = 150;

figure;
gx = geoaxes;
geoscatter(coords_TP(:,1),coords_TP(:,2),marker_size,'black','filled');
hold on;
geoscatter(coords_FP(:,1),coords_FP(:,2),marker_size,blue_is,'filled');
hold on;
geoscatter(coords_FN(:,1),coords_FN(:,2),marker_size,red_is,'filled');
hold on;
geoscatter(coords(start_node,1),coords(start_node,2),marker_size+250,yellow_is,'diamond','filled');
legend('True Positive','False Positive','False Negative','Seed Node', 'Interpreter', 'Latex','Fontsize',16, 'Location','southwest');
gx.Basemap = 'colorterrain';
gx.Scalebar.Visible = 'off';
gx.TickDir = 'none';
gx.LatitudeAxis.TickLabels = {};
gx.LongitudeAxis.TickLabels = {};
gx.LatitudeLabel.String = '';
gx.LongitudeLabel.String = '';
grid off;
set(get(gca, 'Title'), 'String', '\textbf{Seed  node: Rome}', 'Interpreter', 'Latex','Fontsize',24);
geolimits([29.420460,47.650588],[-3.375215,22.736990]);
print('NPR_Orbis_Roma','-dpng','-r200')

end %function