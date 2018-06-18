% Script that allow to analyze the features vectors and order the
% the targetted gene by degree of disruption

% copyright 2017 All rights reserved
%
% Pascaline Parisot (pascaline.parisot.ucl@gmail.com) 
% Christophe De Vleeschouwer (christophe.devleeschouwer@uclouvain.be)
% ISPGroup, Universite catholique de Louvain (Belgium)
% http://sites.uclouvain.be/ispgroup/
%
% Denis L.J. Lafontaine (denis.lafontaine@ulb.ac.be)
% RNA Molecule Biology, Universite Libre de Bruxelles (Belgium)
% http://www.LafontainLab.com
% http://www.RibosomalProteins.com
% http://www.RibosomeSynthesis.com

clear all
close all

addpath('../data')
addpath('../results')


markersize = 10; % default = 6
linewidth = 2; % default = 0.5
fontsize = 20; % default = 10

% Creation of the subdirectory to save the pictures
if ~exist(['../results/images'])
  mkdir('../results/', 'images')
end

REP_IMAGE = '../results/images/';

% Minimum number of nuclei per well to consider the well for the analysis
thres_nb_nuclei_min = 200; 

% load the information about the database 
% => DATA, gene_name, compute_all_features
dataconfig

% Initialization
DATA_gene_index = cat(1,DATA.gene_index);

nb_gene_name = size(gene_name,1);

% Number of nuclei by well
load nb_nuclei_by_well.mat

% All the distances for all features
if (compute_all_features==1)
  load all_the_distances_all.mat
else
  load all_the_distances.mat
end

% the_vec_dist
% row = well index
% column = feature

% The selected features are:
%  1) AA_lcc:   area of the largest connected component
%  2) SE_lcc:   elliptical regularity of the largest connected component
%  3) TH_lcc:   percent of pixels with an intensity smaller the X after 
%               one erosion on the largest connected component
%  4) TLM_lcc:  intensity of the smallest local minimum of 
%               the largest connected component
%  5) TV_lcc:   number of local minima of the largest connected component

%  6) AS:       ratio white/gray
%  7) SE_lcc:   aspect ratio of the largest connected component
%  8) TU_lcc:   greatest number of white areas
%  9) AN_lcc:   number of connected components
% 10) SC_lcc:   concavity of the largest connected component
% 11) TH2_lcc:  percent of pixels with an intensity smaller the X after
%               two erosions on the largest connected component
% 12) TH3_lcc:  percent of pixels with an intensity smaller the X 
%               on the three last erosions on the largest connected component
% 13) TP_lcc:   number of local maxima of the largest connected component

the_distances_vectors = [
  distance_area_lcc ...
  distance_elliptical_regularity_lcc ...
  distance_texture_psX_1eros_lcc ...
  distance_intensity_ppminloc_lcc ...
  distance_nb_minloc_lcc
];
if (compute_all_features==1) % compute all features
  the_distances_vectors = [
    the_distances_vectors ...
    distance_ratio_white_gray ...
    distance_aspect_ratio_lcc ...
    distance_nb_max_white_area ...
    distance_nb_connected_comp ...
    distance_concavity_lcc ...
    distance_texture_psX_2eros_lcc ...
    distance_texture_psX_3lastlevel_eros_lcc ...
    distance_nb_maxloc_lcc
  ];
end

nb_features = size(the_distances_vectors, 2);

the_vec_dist = the_distances_vectors;

color = {'b', 'g', 'r', 'c', 'm', 'y', 'k', [0 0.5 0.8]};
symbol = {'.', 'o', 'x', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p', 'h'};

color_groupe = {'*r', 'db', 'vm', 'hg'};
color_groupe_ko = {'*k', 'dk', 'vk', 'hk'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The features of interest
sum_dist = sum(abs(the_vec_dist),2);
[val_sort ind_sort] = sort(sum_dist);

% sort by well
h=figure;
hold on
j=1;
for i=ind_sort'
  if (isnan(sum_dist(i)) | nb_nuclei_by_well(i)<thres_nb_nuclei_min)
    plot(j, 0, '*k', 'MarkerSize', 10, 'LineWidth', 2)
  else
    plot(j, sum_dist(i), '*r', 'MarkerSize', 10, 'LineWidth', 2)
  end
  j=j+1;
end
hold off
h_ylabel = ylabel('\Delta_k');
title(['Well by well'], 'FontSize', fontsize, 'FontName', 'Times')
axis([0 nb_wells+1 0 1])
axis 'auto y'
set(gca,'XTick',[0:nb_wells+1])
set(gca, 'XTickLabel', [{' ' the_name_legend{ind_sort} ' '}]);
set(gca,'FontSize', fontsize, 'FontName', 'Times')
set(h_ylabel,'FontSize', fontsize, 'FontName', 'Times');

saveas(h,[REP_IMAGE '/fig_distance_sumabsdist_' num2str(nb_features) 'features_fisher.png'])
saveas(h,[REP_IMAGE '/fig_distance_sumabsdist_' num2str(nb_features) 'features_fisher.fig'])

% sort by target gene
mean_sum_dist_gene = zeros(nb_gene_name,1);
for i=1:nb_gene_name
  ind = find(DATA_gene_index==i);
  ind_ok = find(~isnan(sum_dist(ind)) & nb_nuclei_by_well(ind)>thres_nb_nuclei_min);
  ind_ok = ind(ind_ok);
  ind_ko = setdiff(ind,ind_ok);
  [val indx] = sort(sum_dist(ind_ok));
  ind_puits_ok{i} = ind_ok(indx);
  ind_puits_ko{i} = ind_ko;
  mean_sum_dist_gene(i) = mean(sum_dist(ind_ok));
end
[val_sort_mean ind_mean_sort] = sort(mean_sum_dist_gene);

h=figure;
hold on
j=1;
for i=ind_mean_sort'
  if (isnan(mean_sum_dist_gene(i)))
    if ~isempty(ind_puits_ko{i})
    plot(j, 0, '*k', 'MarkerSize', 10, 'LineWidth', 2)
    end
  else
    if ~isempty(ind_puits_ok{i})
      plot(j, sum_dist(ind_puits_ok{i}), '*r', 'MarkerSize', 10, 'LineWidth', 2)
    end
    if ~isempty(ind_puits_ko{i})
      plot(j, sum_dist(ind_puits_ko{i}), '*k', 'MarkerSize', 10, 'LineWidth', 2)
    end 
  end
  j=j+1;
end
hold off
%h_ylabel = ylabel('mean(\Delta_k) per target gene');
h_ylabel = ylabel('Index of nucleolar disruption (iNo) per target gene');
title(['Target gene'], 'FontSize', fontsize, 'FontName', 'Times')
axis([0 nb_gene_name+1 0 1])
axis 'auto y'
set(gca,'XTick',[0:nb_gene_name+1])
set(gca, 'XTickLabel', [{' ' gene_name{ind_mean_sort} ' '}]);
set(gca,'FontSize', fontsize, 'FontName', 'Times')
set(h_ylabel,'FontSize', fontsize, 'FontName', 'Times');

disp('gene name sorted by index of nucleolar disruption (iNo):')
for i=1:nb_gene_name
  disp(sprintf('\t%s \t %2.6f', gene_name{ind_mean_sort(i)}, val_sort_mean(i)))
end

saveas(h,[REP_IMAGE '/fig_distance_sumabsdist_' num2str(nb_features) 'features_fisher_sorted_by_gene.png'])
saveas(h,[REP_IMAGE '/fig_distance_sumabsdist_' num2str(nb_features) 'features_fisher_sorted_by_gene.fig'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA analysis
% Keep only the wells with enough nuclei
ind_OK = find(nb_nuclei_by_well>=thres_nb_nuclei_min & DATA_gene_index>=1);
nb_wells_OK = size(ind_OK,1);

if (nb_wells_OK<30)
  disp(' ');
  disp(['WARNING : Number of data/wells to compute the PCA is not sufficient to make it statistically relevant. Number of data/wells (here ' num2str(nb_wells_OK) ') should be greater than 30.']);
  disp(' ');
end

[COEFF, SCORE, LATENT] = pca(the_vec_dist(ind_OK,:));

h=figure;
hold on
ii_save = [];
for i=1:nb_gene_name
  is = ceil(i/8);
  ic = i-8*(is-1);
  ii = find(i==DATA_gene_index(ind_OK));
  ii_save = [ii_save ; ii];
  for j=ii'
  plot3(SCORE(j,1), SCORE(j,2), SCORE(j,3), 'Color', color{ic}, 'Marker', symbol{is}, 'MarkerSize', 10, 'LineWidth', 2)
  %text(SCORE(j,1), SCORE(j,2), SCORE(j,3), the_name_legend{j}, 'Color', color{ic})
  end
end
hold off
h_xlabel = xlabel('x PCA');
h_ylabel = ylabel('y PCA');
%h_zlabel = zlabel('z PCA')
if (nb_wells_OK<30)
  title(['WARNING : number of data to compute the PCA is not sufficient to make it statistically relevant.'],'Color','r')
end
set(gca,'FontSize', fontsize, 'FontName', 'Times')
set(h_xlabel,'FontSize', fontsize, 'FontName', 'Times');
set(h_ylabel,'FontSize', fontsize, 'FontName', 'Times');
%set(h_zlabel,'FontSize', fontsize, 'FontName', 'Times');

saveas(h,[REP_IMAGE '/fig_distance_pca_' num2str(nb_features) 'features_fisher.png'])
saveas(h,[REP_IMAGE '/fig_distance_pca_' num2str(nb_features) 'features_fisher.fig'])

% Figure with name of gene
h=figure;
hold on
ii_save = [];
for i=1:nb_gene_name
  is = ceil(i/8);
  ic = i-8*(is-1);
  ii = find(i==DATA_gene_index(ind_OK));
  the_name_legend_OK = the_name_legend(ind_OK);
  ii_save = [ii_save ; ii];
  for j=ii'
  plot3(SCORE(j,1), SCORE(j,2), SCORE(j,3), 'Color', color{ic}, 'Marker', symbol{is}, 'MarkerSize', 10, 'LineWidth', 2)
  text(SCORE(j,1), SCORE(j,2), SCORE(j,3), the_name_legend_OK{j}, 'Color', color{ic})
  end
end
hold off
h_xlabel = xlabel('x PCA');
h_ylabel = ylabel('y PCA');
%h_zlabel = zlabel('z PCA')
if (nb_wells_OK<30)
  title(['WARNING : number of data to compute the PCA is not sufficient to make it statistically relevant.'],'Color','r')
end
set(gca,'FontSize', fontsize, 'FontName', 'Times')
set(h_xlabel,'FontSize', fontsize, 'FontName', 'Times');
set(h_ylabel,'FontSize', fontsize, 'FontName', 'Times');
%set(h_zlabel,'FontSize', fontsize, 'FontName', 'Times');

saveas(h,[REP_IMAGE '/fig_distance_pca_' num2str(nb_features) 'features_fisher_with_legend.png'])
saveas(h,[REP_IMAGE '/fig_distance_pca_' num2str(nb_features) 'features_fisher_with_legend.fig'])
