% Script that computes, for each well and for each feature, 
% the distance measures between the histogram of that feature 
% from normal nuclei well and the histogram of the same feature 
% from 'abnormal' nuclei well

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

addpath('../data')
addpath('../results')
addpath('tools')

clear all
close all

% load the information about the database
% => DATA, nb_wells, nb_site_per_well, index_normal_wells
%    compute_all_distances
dataconfig

tic
% Initialization
distance_area_lcc                        = zeros(nb_wells,1);
distance_elliptical_regularity_lcc       = zeros(nb_wells,1);
distance_texture_psX_1eros_lcc           = zeros(nb_wells,1);
distance_intensity_ppminloc_lcc          = zeros(nb_wells,1);
distance_nb_minloc_lcc                   = zeros(nb_wells,1);

if (compute_all_features==1)
  distance_ratio_white_gray                = zeros(nb_wells,1);
  distance_aspect_ratio_lcc                = zeros(nb_wells,1);
  distance_nb_max_white_area               = zeros(nb_wells,1);
  distance_nb_connected_comp               = zeros(nb_wells,1);
  distance_concavity_lcc                   = zeros(nb_wells,1);
  distance_texture_psX_2eros_lcc           = zeros(nb_wells,1);
  distance_texture_psX_3lastlevel_eros_lcc = zeros(nb_wells,1);
  distance_nb_maxloc_lcc                   = zeros(nb_wells,1);
  distance_intensity_ppmaxloc_lcc          = zeros(nb_wells,1);
end % compute all features

% All the computations
% 'normal' nuclei
the_area_lcc_normal                        = [];
the_elliptical_regularity_lcc_normal       = [];
the_texture_psX_1eros_lcc_normal           = [];
the_intensity_ppminloc_lcc_normal          = [];
the_nb_minloc_lcc_normal                   = [];
the_ratio_white_gray_normal                = [];
the_aspect_ratio_lcc_normal                = [];
the_nb_max_white_area_normal               = [];
the_nb_connected_comp_normal               = [];
the_concavity_lcc_normal                   = [];
the_texture_psX_2eros_lcc_normal           = [];
the_texture_psX_3lastlevel_eros_lcc_normal = [];
the_nb_maxloc_lcc_normal                   = [];
the_intensity_ppmaxloc_lcc_normal          = [];

for wellIndex=index_normal_wells

  plateName = DATA(wellIndex).plate;
  wellName = DATA(wellIndex).well;

  %disp(['normal : ' plateName wellName])
  disp(['reference : ' plateName wellName])

  for siteIndex = 1:nb_site_per_well

    siteName = ['s' num2str(siteIndex)];

    if (compute_all_features==1)
      load(['../results/features/' plateName '/res_' plateName wellName siteName '_features_all.mat']);
    else
      load(['../results/features/' plateName '/res_' plateName wellName siteName '_features.mat']);
    end

    the_area_lcc_normal                        = [the_area_lcc_normal area_lcc];
    the_elliptical_regularity_lcc_normal       = [the_elliptical_regularity_lcc_normal elliptical_regularity_lcc];
    the_texture_psX_1eros_lcc_normal           = [the_texture_psX_1eros_lcc_normal texture_psX_1eros_lcc];
    the_intensity_ppminloc_lcc_normal          = [the_intensity_ppminloc_lcc_normal intensity_ppminloc_lcc];
    the_nb_minloc_lcc_normal                   = [the_nb_minloc_lcc_normal nb_minloc_lcc];

    if (compute_all_features==1)
      the_ratio_white_gray_normal                = [the_ratio_white_gray_normal ratio_white_gray];
      the_aspect_ratio_lcc_normal                = [the_aspect_ratio_lcc_normal aspect_ratio_lcc];
      the_nb_max_white_area_normal               = [the_nb_max_white_area_normal nb_max_white_area];
      the_nb_connected_comp_normal               = [the_nb_connected_comp_normal nb_connected_comp];
      the_concavity_lcc_normal                   = [the_concavity_lcc_normal concavity_lcc];
      the_texture_psX_2eros_lcc_normal           = [the_texture_psX_2eros_lcc_normal texture_psX_2eros_lcc];
      the_texture_psX_3lastlevel_eros_lcc_normal = [the_texture_psX_3lastlevel_eros_lcc_normal texture_psX_3lastlevel_eros_lcc];
      the_nb_maxloc_lcc_normal                   = [the_nb_maxloc_lcc_normal nb_maxloc_lcc]; 
      the_intensity_ppmaxloc_lcc_normal          = [the_intensity_ppmaxloc_lcc_normal intensity_ppmaxloc_lcc];
    end % compute all features

  end
end


% the normal or abnormal wells
for wellIndex=1:nb_wells

  the_area_lcc_abnormal                        = [];
  the_elliptical_regularity_lcc_abnormal       = [];
  the_texture_psX_1eros_lcc_abnormal           = [];
  the_intensity_ppminloc_lcc_abnormal          = [];
  the_nb_minloc_lcc_abnormal                   = [];

  if (compute_all_features==1)
    the_ratio_white_gray_abnormal                = [];
    the_aspect_ratio_lcc_abnormal                = [];
    the_nb_max_white_area_abnormal               = [];
    the_nb_connected_comp_abnormal               = [];
    the_concavity_lcc_abnormal                   = [];
    the_texture_psX_2eros_lcc_abnormal           = [];
    the_texture_psX_3lastlevel_eros_lcc_abnormal = [];
    the_nb_maxloc_lcc_abnormal                   = [];
    the_intensity_ppmaxloc_lcc_abnormal          = [];
  end % compute all features

  plateName = DATA(wellIndex).plate;
  wellName = DATA(wellIndex).well;

  %disp(['abnormal : ' plateName wellName])
  disp(['test : ' plateName wellName])

  for siteIndex = 1:nb_site_per_well

    siteName = ['s' num2str(siteIndex)];

    if (compute_all_features==1)
      load(['../results/features/' plateName '/res_' plateName wellName siteName '_features_all.mat']);
    else
      load(['../results/features/' plateName '/res_' plateName wellName siteName '_features.mat']);
    end

    the_area_lcc_abnormal                        = [the_area_lcc_abnormal area_lcc];
    the_elliptical_regularity_lcc_abnormal       = [the_elliptical_regularity_lcc_abnormal elliptical_regularity_lcc];
    the_texture_psX_1eros_lcc_abnormal           = [the_texture_psX_1eros_lcc_abnormal texture_psX_1eros_lcc];
    the_intensity_ppminloc_lcc_abnormal          = [the_intensity_ppminloc_lcc_abnormal intensity_ppminloc_lcc];
    the_nb_minloc_lcc_abnormal                   = [the_nb_minloc_lcc_abnormal nb_minloc_lcc];

    if (compute_all_features==1)
      the_ratio_white_gray_abnormal                = [the_ratio_white_gray_abnormal ratio_white_gray];
      the_aspect_ratio_lcc_abnormal                = [the_aspect_ratio_lcc_abnormal aspect_ratio_lcc];
      the_nb_max_white_area_abnormal               = [the_nb_max_white_area_abnormal nb_max_white_area];
      the_nb_connected_comp_abnormal               = [the_nb_connected_comp_abnormal nb_connected_comp];
      the_concavity_lcc_abnormal                   = [the_concavity_lcc_abnormal concavity_lcc];
      the_texture_psX_2eros_lcc_abnormal           = [the_texture_psX_2eros_lcc_abnormal texture_psX_2eros_lcc];
      the_texture_psX_3lastlevel_eros_lcc_abnormal = [the_texture_psX_3lastlevel_eros_lcc_abnormal texture_psX_3lastlevel_eros_lcc];
      the_nb_maxloc_lcc_abnormal                   = [the_nb_maxloc_lcc_abnormal nb_maxloc_lcc];
      the_intensity_ppmaxloc_lcc_abnormal          = [the_intensity_ppmaxloc_lcc_abnormal intensity_ppmaxloc_lcc];
    end % compute all features

  end

  % Compute distances...
  distance_area_lcc(wellIndex)                        = f_get_fisher_criterion(the_area_lcc_normal, the_area_lcc_abnormal);
  distance_elliptical_regularity_lcc(wellIndex)       = f_get_fisher_criterion(the_elliptical_regularity_lcc_normal, the_elliptical_regularity_lcc_abnormal);
  distance_texture_psX_1eros_lcc(wellIndex)           = f_get_fisher_criterion(the_texture_psX_1eros_lcc_normal, the_texture_psX_1eros_lcc_abnormal);
  distance_intensity_ppminloc_lcc(wellIndex)          = f_get_fisher_criterion(the_intensity_ppminloc_lcc_normal, the_intensity_ppminloc_lcc_abnormal);
  distance_nb_minloc_lcc(wellIndex)                   = f_get_fisher_criterion(the_nb_minloc_lcc_normal, the_nb_minloc_lcc_abnormal);
  if (compute_all_features==1)
    distance_ratio_white_gray(wellIndex)                = f_get_fisher_criterion(the_ratio_white_gray_normal, the_ratio_white_gray_abnormal);
    distance_aspect_ratio_lcc(wellIndex)                = f_get_fisher_criterion(the_aspect_ratio_lcc_normal, the_aspect_ratio_lcc_abnormal);
    distance_nb_max_white_area(wellIndex)               = f_get_fisher_criterion(the_nb_max_white_area_normal, the_nb_max_white_area_abnormal);
    distance_nb_connected_comp(wellIndex)               = f_get_fisher_criterion(the_nb_connected_comp_normal, the_nb_connected_comp_abnormal);
    distance_concavity_lcc(wellIndex)                   = f_get_fisher_criterion(the_concavity_lcc_normal, the_concavity_lcc_abnormal);
    distance_texture_psX_2eros_lcc(wellIndex)           = f_get_fisher_criterion(the_texture_psX_2eros_lcc_normal, the_texture_psX_2eros_lcc_abnormal);
    distance_texture_psX_3lastlevel_eros_lcc(wellIndex) = f_get_fisher_criterion(the_texture_psX_3lastlevel_eros_lcc_normal, the_texture_psX_3lastlevel_eros_lcc_abnormal);
    distance_nb_maxloc_lcc(wellIndex)                   = f_get_fisher_criterion(the_nb_maxloc_lcc_normal, the_nb_maxloc_lcc_abnormal);
    distance_intensity_ppmaxloc_lcc(wellIndex)          = f_get_fisher_criterion(the_intensity_ppmaxloc_lcc_normal, the_intensity_ppmaxloc_lcc_abnormal);
  end % compute all distances

end

toc

if (compute_all_features==1)
  save(['../results/all_the_distances_all.mat'], 'distance_area_lcc', 'distance_elliptical_regularity_lcc', 'distance_texture_psX_1eros_lcc', 'distance_intensity_ppminloc_lcc', 'distance_nb_minloc_lcc', 'distance_ratio_white_gray', 'distance_aspect_ratio_lcc', 'distance_nb_max_white_area', 'distance_nb_connected_comp', 'distance_concavity_lcc', 'distance_texture_psX_2eros_lcc', 'distance_texture_psX_3lastlevel_eros_lcc', 'distance_nb_maxloc_lcc', 'distance_intensity_ppmaxloc_lcc');
else
  save(['../results/all_the_distances.mat'], 'distance_area_lcc', 'distance_elliptical_regularity_lcc', 'distance_texture_psX_1eros_lcc', 'distance_intensity_ppminloc_lcc', 'distance_nb_minloc_lcc');
end
