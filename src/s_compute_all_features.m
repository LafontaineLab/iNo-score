% Script that computes all the features for all the images
% and save the features for each well in a file

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

% load the information about the database
% => DATA, nb_wells, nb_site_per_well
dataconfig

% Creation of the subdirectories to save the results
if ~exist(['../results/features'])
  mkdir('../results/', 'features')
end
for plateNameIndex=1:nb_wells
  if ~exist(['../results/features/' DATA(plateNameIndex).plate])
    mkdir('../results/features/', DATA(plateNameIndex).plate)
  end
end

% Number of nuclei by well
nb_nuclei_by_well = zeros(nb_wells,1);

% Initialization of the thresholds
% The values of those thresholds have been optimized based on 
% the analysis of different values on those thresholds such that 
% to maximize the distance between the distribution of one feature 
% for "abnormal" wells compared to the distribution of the "normal" wells

thres_connexity_area = 100;             % AA_lcc (\tau_a) area
thres_connexity_shape = 110;            % SE_lcc (\tau_s) aspect ratio
                                        % SR_lcc (\tau_s) elliptical regularity
                                        % SC_lcc (\tau_s) concavity
thres_connexity_texture_psX = 100;      % TH_lcc (\tau_t) 1 eros
                                        % TH2_lcc (\tau_t) 2 eros
                                        % TH3_lcc (tau_t) 3 last level eros
thres_connexity_minmaxloc = 95;         % TLM_lcc (\tau_t) ppminloc
                                        % TV_lcc (\tau_t) nb_min_loc
                                        % TLMM_lcc (\tau_t) ppmaxloc
                                        % TP_lcc (\tau_t) nb_max_loc
thres_intensity_1eros = 200;            % TH_lcc (\alpha) 1 eros

if (compute_all_features==1)
  thres_connexity_nb_connected_comp = 95; % AN_lcc (\tau_a)
  thres_connexity_white_area = 100;       % TU_lcc (\tau_t)
  thres_gray = 100;                       % AS (\tau_i1)
  thres_white = 230;                      % AS (\tau_i2)
  thres_intensity_2eros = 210;            % TH2_lcc (\alpha_2)
  thres_intensity_3derniv_eros = 230;     % TH3_lcc (\alpha_3)
  thres_white_area = 200;                 % TU_lcc (\beta)
end

percentile_for_shape = 99.9;

if (compute_all_features==0)

  for wellIndex=1:nb_wells

    plateName = DATA(wellIndex).plate;
    wellName = DATA(wellIndex).well;

    for siteIndex = 1:nb_site_per_well

      siteName = ['s' num2str(siteIndex)];

      disp(['Compute features: ' plateName ' ' wellName ' ' siteName]);
tic
      [ nb_nuclei_valid, area_lcc, ...
	elliptical_regularity_lcc, texture_psX_1eros_lcc, ...
	intensity_ppminloc_lcc, nb_minloc_lcc ...
      ] = f_get_nucleus_features(compute_all_features, ...
                                 plateName, wellName, siteName, ...
				 percentile_for_shape, ...
				 thres_connexity_area, ...
				 thres_connexity_shape, ...
				 thres_connexity_texture_psX, ...
				 thres_connexity_minmaxloc, ...
				 thres_intensity_1eros);
toc
      nb_nuclei_by_well(wellIndex) = nb_nuclei_by_well(wellIndex) + nb_nuclei_valid;

      save(['../results/features/' plateName '/res_'  plateName wellName siteName '_features.mat'], 'nb_nuclei_valid', 'area_lcc', 'elliptical_regularity_lcc', 'texture_psX_1eros_lcc', 'intensity_ppminloc_lcc', 'nb_minloc_lcc')
 
    end % for site indexes
  end % for well indexes

else % compute all features

  for wellIndex=1:nb_wells

    plateName = DATA(wellIndex).plate;
    wellName = DATA(wellIndex).well;

    for siteIndex = 1:nb_site_per_well

      siteName = ['s' num2str(siteIndex)];

      disp(['Compute features: ' plateName ' ' wellName ' ' siteName]);
tic
      [ nb_nuclei_valid, area_lcc, ...
	elliptical_regularity_lcc, texture_psX_1eros_lcc, ...
	intensity_ppminloc_lcc, nb_minloc_lcc, ...
        ratio_white_gray, aspect_ratio_lcc, ...
        nb_max_white_area, nb_connected_comp, concavity_lcc, ...
        texture_psX_2eros_lcc, texture_psX_3lastlevel_eros_lcc, ...
	nb_maxloc_lcc, intensity_ppmaxloc_lcc ...
      ] = f_get_nucleus_features(compute_all_features, ...
                                 plateName, wellName, siteName, ...
				 percentile_for_shape, ...
				 thres_connexity_area, ...
				 thres_connexity_shape, ...
				 thres_connexity_texture_psX, ...
				 thres_connexity_minmaxloc, ...
				 thres_intensity_1eros, ...
				 thres_connexity_nb_connected_comp, ...
				 thres_connexity_white_area, ...
				 thres_gray, thres_white, ...
				 thres_intensity_2eros, ...
				 thres_intensity_3derniv_eros, ...
				 thres_white_area);
toc
      nb_nuclei_by_well(wellIndex) = nb_nuclei_by_well(wellIndex) + nb_nuclei_valid;

      save(['../results/features/' plateName '/res_'  plateName wellName siteName '_features_all.mat'], 'nb_nuclei_valid', 'area_lcc', 'elliptical_regularity_lcc', 'texture_psX_1eros_lcc', 'intensity_ppminloc_lcc', 'nb_minloc_lcc', 'ratio_white_gray', 'aspect_ratio_lcc', 'nb_max_white_area', 'nb_connected_comp', 'concavity_lcc', 'texture_psX_2eros_lcc', 'texture_psX_3lastlevel_eros_lcc', 'nb_maxloc_lcc', 'intensity_ppmaxloc_lcc')
 
    end % for site indexes
  end % for well indexes

end

save(['../results/nb_nuclei_by_well.mat'], 'nb_nuclei_by_well')

% Display number of nuclei by well
disp('Number of nuclei by well:')
for wellIndex=1:nb_wells
  disp(sprintf('%s%s:\t%6d', DATA(wellIndex).plate, DATA(wellIndex).well, nb_nuclei_by_well(wellIndex)))
end
