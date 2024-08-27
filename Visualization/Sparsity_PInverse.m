function Sparsity_PInverse(Inv,factor)
% Visualize the sparsity pattern of
% the pseudoinverse of the Incidence



figure;
% imagesc(log10(abs(Inv)+eps));
   
imagesc(Inv>mean(mean(abs(Inv)/factor)));
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
axis square
colormap(flipud(bone))
colorbar;
set(gca,'Fontsize',15);
set(gcf,'position',[0 0 1.2 1]*250)
tightfig;
% saveas(gcf,'Theta_inv_clust_proper','pdf');
%%%%%%%%%%%%%%%%

end