% Script that computes the features and the distances
% to optimize the parameters

% IMPORTANT: It is assumed that the segmentation of the images has been done.
% (see script s_compute_segmentation.m)

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
addpath('tools')

COMPUTE_FEATURES = 1;
COMPUTE_DISTANCES = 1;

% load the information about the database
% => DATA, index_normal_wells, compute_all_features
dataconfig

% Index of wells for the optimization of the parameters
index_wells_for_opti = [3 5 8 13 16];
nb_wells_for_opti = size(index_wells_for_opti,2);

% Creation of the subdirectories to save the results
if ~exist(['../results/features_for_parameters_optimization'])
  mkdir('../results/', 'features_for_parameters_optimization')
end
for plateNameIndex=[index_normal_wells index_wells_for_opti]
  if ~exist(['../results/features_for_parameters_optimization/' DATA(plateNameIndex).plate])
    mkdir('../results/features_for_parameters_optimization/', DATA(plateNameIndex).plate)
  end
end

percentile_for_shape = 99.9;

step = 2;

connexity_threshold_range  =  60:step:200;
intensity_threshold_range  = 100:step:255;
if (compute_all_features==1)
  gray_threshold_range       =  60:step:250;
  white_threshold_range      =  60:step:255;
  white_area_threshold_range = 100:step:250;
end % compute all features

% For the saving of the features
the_step = step;
the_connexity_threshold_range  = connexity_threshold_range;
the_intensity_threshold_range  = intensity_threshold_range;
if (compute_all_features==1)
  the_gray_threshold_range       = gray_threshold_range;
  the_white_threshold_range      = white_threshold_range;
  the_white_area_threshold_range = white_area_threshold_range;
end % compute all features

if (COMPUTE_FEATURES)

  if (compute_all_features==0)

    for wellIndex=[index_normal_wells index_wells_for_opti]

      plateName = DATA(wellIndex).plate;
      wellName = DATA(wellIndex).well;

      for siteIndex = 1:nb_site_per_well

        siteName = ['s' num2str(siteIndex)];

        disp(['Compute features for opti: ' plateName ' ' wellName ' ' siteName]);

tic
        [ nb_nuclei_valid, area_lcc, ...
          elliptical_regularity_lcc, texture_psX_1eros_lcc, ...
          intensity_ppminloc_lcc, nb_minloc_lcc ...
        ] = f_get_nucleus_features_for_parameters_optimization(compute_all_features, ...
                                                               plateName, wellName, siteName, ...
							       percentile_for_shape, ...
							       connexity_threshold_range, ...
							       intensity_threshold_range);
toc

        save(['../results/features_for_parameters_optimization/' plateName '/res_'  plateName wellName siteName '_features_for_parameters_optimization.mat'], 'nb_nuclei_valid', 'area_lcc', 'elliptical_regularity_lcc', 'texture_psX_1eros_lcc', 'intensity_ppminloc_lcc', 'nb_minloc_lcc', 'the_connexity_threshold_range', 'the_intensity_threshold_range', 'the_step')
 
      end
    end

  else % compute all features

    for wellIndex=[index_normal_wells index_wells_for_opti]

      plateName = DATA(wellIndex).plate;
      wellName = DATA(wellIndex).well;

      for siteIndex = 1:nb_site_per_well

        siteName = ['s' num2str(siteIndex)];

        disp(['Compute features for opti: ' plateName ' ' wellName ' ' siteName]);
tic
        [ nb_nuclei_valid, area_lcc, ...
          elliptical_regularity_lcc, texture_psX_1eros_lcc, ...
          intensity_ppminloc_lcc, nb_minloc_lcc, ...
          ratio_white_gray, aspect_ratio_lcc, ...
          nb_max_white_area, nb_connected_comp, concavity_lcc, ...
          texture_psX_2eros_lcc, texture_psX_3lastlevel_eros_lcc, ...
          nb_maxloc_lcc, intensity_ppmaxloc_lcc ...
        ] = f_get_nucleus_features_for_parameters_optimization(compute_all_features, ...
                                                               plateName, wellName, siteName, ...
							       percentile_for_shape, ...
							       connexity_threshold_range, ...
							       intensity_threshold_range, ...
							       gray_threshold_range, ...
							       white_threshold_range, ...
							       white_area_threshold_range);
toc

    save(['../results/features_for_parameters_optimization/' plateName '/res_'  plateName wellName siteName '_features_for_parameters_optimization_all.mat'], 'nb_nuclei_valid', 'area_lcc', 'elliptical_regularity_lcc', 'texture_psX_1eros_lcc', 'intensity_ppminloc_lcc', 'nb_minloc_lcc', 'ratio_white_gray', 'aspect_ratio_lcc', 'nb_max_white_area', 'nb_connected_comp', 'concavity_lcc', 'texture_psX_2eros_lcc', 'texture_psX_3lastlevel_eros_lcc', 'nb_maxloc_lcc', 'intensity_ppmaxloc_lcc', 'the_connexity_threshold_range', 'the_intensity_threshold_range', 'the_gray_threshold_range', 'the_white_threshold_range', 'the_white_area_threshold_range', 'the_step')
 
      end
    end

  end % compute all features

end % COMPUTE_FEATURES

% Compute all distances
if (COMPUTE_DISTANCES)

  disp(['Compute distances for opti: ']);

  if (compute_all_features==0)
tic
    [ distance_area_lcc_for_opti, ...
      distance_elliptical_regularity_lcc_for_opti, ...
      distance_texture_psX_1eros_lcc_for_opti, ... 
      distance_intensity_ppminloc_lcc_for_opti, ...
      distance_nb_minloc_lcc_for_opti ...
    ] = f_get_distances_for_parameters_optimization(compute_all_features, ...
						    DATA, index_wells_for_opti, ...
						    index_normal_wells, ...
						    nb_site_per_well, ...
						    connexity_threshold_range, ...
						    intensity_threshold_range);
toc

    save(['../results/all_the_distances_for_parameters_optimization.mat'], 'index_normal_wells', 'index_wells_for_opti', 'distance_area_lcc_for_opti', 'distance_elliptical_regularity_lcc_for_opti', 'distance_texture_psX_1eros_lcc_for_opti', 'distance_intensity_ppminloc_lcc_for_opti', 'distance_nb_minloc_lcc_for_opti', 'connexity_threshold_range', 'intensity_threshold_range', 'step');

  else % compute all features

tic
    [ distance_area_lcc_for_opti, ...
      distance_elliptical_regularity_lcc_for_opti, ...
      distance_texture_psX_1eros_lcc_for_opti, ... 
      distance_intensity_ppminloc_lcc_for_opti, ...
      distance_nb_minloc_lcc_for_opti, ...
      distance_ratio_white_gray_for_opti, ...
      distance_aspect_ratio_lcc_for_opti, ...
      distance_nb_max_white_area_for_opti, ...
      distance_nb_connected_comp_for_opti, ...
      distance_concavity_lcc_for_opti, ...
      distance_texture_psX_2eros_lcc_for_opti, ...
      distance_texture_psX_3lastlevel_eros_lcc_for_opti, ...
      distance_nb_maxloc_lcc_for_opti, ...
      distance_intensity_ppmaxloc_lcc_for_opti ...
    ] = f_get_distances_for_parameters_optimization(compute_all_features, ...
						    DATA, index_wells_for_opti, ...
						    index_normal_wells, ...
						    nb_site_per_well, ...
						    connexity_threshold_range, ...
						    intensity_threshold_range, ...
						    gray_threshold_range, ...
						    white_threshold_range, ...
						    white_area_threshold_range);
toc

    save(['../results/all_the_distances_for_parameters_optimization_all.mat'], 'index_normal_wells', 'index_wells_for_opti', 'distance_area_lcc_for_opti', 'distance_elliptical_regularity_lcc_for_opti', 'distance_texture_psX_1eros_lcc_for_opti', 'distance_intensity_ppminloc_lcc_for_opti', 'distance_nb_minloc_lcc_for_opti', 'distance_ratio_white_gray_for_opti', 'distance_aspect_ratio_lcc_for_opti', 'distance_nb_max_white_area_for_opti', 'distance_nb_connected_comp_for_opti', 'distance_concavity_lcc_for_opti', 'distance_texture_psX_2eros_lcc_for_opti', 'distance_texture_psX_3lastlevel_eros_lcc_for_opti', 'distance_nb_maxloc_lcc_for_opti', 'distance_intensity_ppmaxloc_lcc_for_opti', 'connexity_threshold_range', 'intensity_threshold_range', 'gray_threshold_range', 'white_threshold_range', 'white_area_threshold_range', 'step');

  end % compute all features

end % COMPUTE_DISTANCES

%% Display the distance to manually select the thresholds
SAVE_REP  = '../results/figures_parameters_optimization/';
if ~exist(SAVE_REP)
  mkdir('../results/', 'figures_parameters_optimization')
end

% load all the distances
if (compute_all_features==0)
  load(['../results/all_the_distances_for_parameters_optimization.mat']);
else
  load(['../results/all_the_distances_for_parameters_optimization_all.mat']);
end

supAxes = [.08 .08 .87 .87]; 
the_titles{ 1} = 'AA_{lcc} area lcc';                          the_features_name{ 1} = 'area_lcc'; 
the_titles{ 2} = 'SR_{lcc} elliptical regularity lcc';         the_features_name{ 2} = 'elliptical_regularity_lcc';
the_titles{ 3} = 'TH_{lcc} texture psX 1eros lcc';             the_features_name{ 3} = 'texture_psX_1eros_lcc';
the_titles{ 4} = 'TLM_{lcc} intensity ppminloc lcc';           the_features_name{ 4} = 'intensity_ppminloc_lcc';
the_titles{ 5} = 'TV_{lcc} nb minloc lcc';                     the_features_name{ 5} = 'nb_minloc_lcc';
if (compute_all_features==1)
  the_titles{ 6} = 'AS ratio white gray';                        the_features_name{ 6} = 'ratio_white_gray';
  the_titles{ 7} = 'SE_{lcc} aspect ratio lcc';                  the_features_name{ 7} = 'aspect_ratio_lcc';
  the_titles{ 8} = 'TU_{lcc} nb max white area';                 the_features_name{ 8} = 'nb_max_white_area';
  the_titles{ 9} = 'AN_{lcc} nb connected component';            the_features_name{ 9} = 'nb_connected_component';
  the_titles{10} = 'SC_{lcc} concavity lcc';                     the_features_name{10} = 'concavity_lcc';
  the_titles{11} = 'TH2_{lcc} texture psX 2eros lcc';            the_features_name{11} = 'texture_psX_2eros_lcc';
  the_titles{12} = 'TH3_{lcc} texture psX 3lastlevel eros lcc';  the_features_name{12} = 'texture_psX_3lastlevel_eros_lcc';
  the_titles{13} = 'TP_{lcc} nb maxloc lcc';                     the_features_name{13} = 'nb_maxloc_lcc';
  the_titles{14} = 'TLMM_{lcc} intensity ppmaxloc lcc';          the_features_name{14} = 'intensity_ppmaxloc_lcc';
end % compute all features

nb_ind_title = size(the_titles,2);
ind_title = [1:nb_ind_title]

nb_connexity_threshold  = size(connexity_threshold_range, 2);
nb_intensity_threshold  = size(intensity_threshold_range, 2);
if (compute_all_features==1)
  nb_gray_threshold       = size(gray_threshold_range, 2);
  nb_white_threshold      = size(white_threshold_range, 2);
  nb_white_area_threshold = size(white_area_threshold_range, 2);
end

nbl = ceil(sqrt(nb_wells_for_opti));
nbc = nbl;

for ind=1:nb_wells_for_opti

  wellIndex = index_wells_for_opti(ind);
  plateName = DATA(wellIndex).plate;
  wellName  = DATA(wellIndex).well;
  
  figure(1)
  subplot(nbl,nbc,ind)
  plot(connexity_threshold_range, abs(distance_area_lcc_for_opti{ind}))
  xlabel('thres connexity area (\tau_a)')
  title([plateName wellName])
   
  figure(2)
  subplot(nbl,nbc,ind)
  plot(connexity_threshold_range, abs(distance_elliptical_regularity_lcc_for_opti{ind}))
  xlabel('thres connexity shape (\tau_s)')
  title([plateName wellName])

  figure(3)
  subplot(nbl,nbc,ind)
  mesh(intensity_threshold_range, connexity_threshold_range, abs(distance_texture_psX_1eros_lcc_for_opti{ind}))
  xlabel('thres intensity 1eros (\alpha)')
  ylabel('thres connexity texture psX (\tau_t)')
  title([plateName wellName])
   
  figure(4)
  subplot(nbl,nbc,ind)
  plot(connexity_threshold_range, abs(distance_intensity_ppminloc_lcc_for_opti{ind}))
  xlabel('thres connexity minmaxloc (\tau_t)')
  title([plateName wellName])
   
  figure(5)
  subplot(nbl,nbc,ind)
  plot(connexity_threshold_range, abs(distance_nb_minloc_lcc_for_opti{ind}))
  xlabel('thres connexity minmaxloc (\tau_t)')
  title([plateName wellName])

  if (compute_all_features==1)

    figure(6)
    subplot(nbl,nbc,ind)
    mesh(white_threshold_range, gray_threshold_range, abs(distance_ratio_white_gray_for_opti{ind}))
    xlabel('thres white (\tau_{i2})')
    ylabel('thres gray (\tau_{i1})')
    title([plateName wellName])
   
    figure(7)
    subplot(nbl,nbc,ind)
    plot(connexity_threshold_range, abs(distance_aspect_ratio_lcc_for_opti{ind}))
    xlabel('thres connexity shape (\tau_s)')
    title([plateName wellName])

    figure(8)
    subplot(nbl,nbc,ind)
    mesh(white_area_threshold_range, connexity_threshold_range, abs(distance_nb_max_white_area_for_opti{ind}))
    xlabel('thres white area (\beta)')
    ylabel('thres connexity white area (\tau_t)')
    title([plateName wellName])

    figure(9)
    subplot(nbl,nbc,ind)
    plot(connexity_threshold_range, abs(distance_nb_connected_comp_for_opti{ind}))
    xlabel('thres connexity nb connected comp (\tau_a)')
    title([plateName wellName])
   
    figure(10)
    subplot(nbl,nbc,ind)
    plot(connexity_threshold_range, abs(distance_concavity_lcc_for_opti{ind}))
    xlabel('thres connexity shape (\tau_s)')
    title([plateName wellName])

    figure(11)
    subplot(nbl,nbc,ind)
    mesh(intensity_threshold_range, connexity_threshold_range, abs(distance_texture_psX_2eros_lcc_for_opti{ind}))
    xlabel('thres intensity 2eros (\alpha 2)')
    ylabel('thres connexity texture psX (\tau_t)')
    title([plateName wellName])

    figure(12)
    subplot(nbl,nbc,ind)
    mesh(intensity_threshold_range, connexity_threshold_range, abs(distance_texture_psX_3lastlevel_eros_lcc_for_opti{ind}))
    xlabel('thres intensity 3lastlevel eros (\alpha 3)')
    ylabel('thres connexity texture psX (\tau_t)')
    title([plateName wellName])

    figure(13)
    subplot(nbl,nbc,ind)
    plot(connexity_threshold_range, abs(distance_nb_maxloc_lcc_for_opti{ind}))
    xlabel('thres connexity minmaxloc (\tau_t)')
    title([plateName wellName])
   
    figure(14)
    subplot(nbl,nbc,ind)
    plot(connexity_threshold_range, abs(distance_intensity_ppmaxloc_lcc_for_opti{ind}))
    xlabel('thres connexity minmaxloc (\tau_t)')
    title([plateName wellName])
  end % compute all features
end


for i=1:nb_ind_title
  figure(i)
  [ax,h3]=suplabel(the_titles{ind_title(i)}, 't', supAxes);

  saveas(i,[SAVE_REP '/opti_threshold_feature_' the_features_name{i} '.fig'])

end
