% Function that computes for each well, for each feature, 
% and for each threshold,
% the distance measures between the histogram of that feature 
% from normal nuclei well and the histogram of the same feature 
% from 'abnormal' nuclei well
% (see the descriptions of the outputs)

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

function [ distance_area_lcc_for_opti, ...
           distance_elliptical_regularity_lcc_for_opti, ...
           distance_texture_psX_1eros_lcc_for_opti, ...
           distance_intensity_ppminloc_lcc_for_opti, ...
           distance_nb_minloc_lcc_for_opti ...
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
                                                         data, index_wells_for_opti, ...
							 index_normal_wells_for_opti, ...
							 nb_site_per_well, ...
							 connexity_threshold_range, ...
							 intensity_threshold_range, ...
							 gray_threshold_range, ...
							 white_threshold_range, ...
							 white_area_threshold_range)


% Inputs:
% *******
% compute_all_features              : if 0 => compute the 5 features
%                                     if 1 => compute all the features
% data                              : cell of structures that contains the data
%                                     (plate name, well name, gene name)
%                                     to deduce the name of the files
% index_wells_for_opti              : array of wells index among the indexes
%                                     in the data structure
% index_normal_wells_for_opti       : array of "normal" wells index among the indexes
%                                     in the data structure
% nb_site_per_well                  : number of site per well  
%                                     (image number)
% connexity_threshold_range         : range of "c" thresholds for the connexity
% intensity_threshold_range         : range of "i" thresholds for the intensity value
%
%% if compute all the features the following parameters must be defined
%
% gray_threshold_range              : range of "g" thresholds for the gray value
% white_threshold_range             : range of "w" thresholds for the white value
% white_area_threshold_range        : range of "wa" thresholds for the definition of the white area
%
% Outputs:
% ********
% distance_area_lcc_for_opti                        : for each connexity threshold,
%                                                     distance for the size of the biggest 
%                                                     connected component (c x 1 vector)
% distance_elliptical_regularity_lcc_for_opti       : for each connexity threshold,
%                                                     distance for the elliptical regularity 
%                                                     measure of the biggest connected
%                                                     component (c x 1 vector)
% distance_texture_psX_1eros_lcc_for_opti           : for each connexity and intensity thresholds,
%                                                     distance for the percent of pixels with 
%                                                     an intensity smaller than X after one 
%                                                     erosion of the biggest connected component
%                                                     (c x i matrix)
% distance_intensity_ppminloc_lcc_for_opti          : for each connexity threshold,
%                                                     distance for the smallest intensity 
%                                                     of the local minima of the biggest 
%                                                     connected component (c x 1 vector)
% distance_nb_minloc_lcc_for_opti                   : for each connexity threshold,
%                                                     distance for the number of local minima 
%                                                     of the biggest connected component 
%                                                     (c x 1 vector)
%
%% if compute all the features the following parameters are computed
%
% distance_ratio_white_gray_for_opti                : for each gray and white thresholds,
%                                                     distance for the white/gray ratio
%                                                     = number of pixels with an intensity 
%                                                     above a threshold in the connected 
%                                                     components (g x w matrix)
% distance_aspect_ratio_lcc_for_opti                : for each connexity threshold,
%                                                     distance for the aspect ration measure 
%                                                     of the biggest connected component
%                                                     (c x 1 vector)
% distance_nb_max_white_area_for_opti               : for each connexity and white area thresholds,
%                                                     distance for the maximal number of 
%                                                     white areas in a connected component 
%                                                     (c x wa matrix)
% distance_nb_connected_comp_for_opti               : for each connexity threshold,
%                                                     distance for the number of nucleolus 
%                                                     connected components (c x 1 vector)
% distance_concavity_lcc_for_opti                   : for each connexity threshold,
%                                                     distance for the concavity measure 
%                                                     of the biggest connected component
%                                                     (c x 1 vector)
% distance_texture_psX_2eros_lcc_for_opti           : for each connexity and intensity thresholds,
%                                                     distance for the percent of pixels with 
%                                                     an intensity smaller than X after two
%                                                     erosions of the biggest connected
%                                                     component (c x i matrix) 
% distance_texture_psX_3lastlevel_eros_lcc_for_opti : for each connexity and intensity thresholds,
%                                                     distance for the percent of pixels with 
%                                                     an intensity smaller than X on the last 
%                                                     three levels of erosions of the biggest 
%                                                     connected component (c x i matrix)
% distance_nb_maxloc_lcc_for_opti                   : for each connexity threshold,
%                                                     distance for the number of local maxima 
%                                                     of the biggest connected component 
%                                                     (c x 1 vector)
% distance_intensity_ppmaxloc_lcc_for_opti          : for each connexity threshold,
%                                                     distance for the smallest intensity 
%                                                     of the local maxima of the biggest 
%                                                     connected component (c x 1 vector)


  % Check the input and output parameters
  if (compute_all_features==0) % compute the 5 features
    if ((nargin~=7) || (nargout~=5))
      error(['invalid number of input ' num2str(nargin) ' (should be 10) or output ' num2str(nargout) ' (should the 6) parameters']);
    end
  else % compute all the features
    if ((nargin~=10) || (nargout~=14))
      error(['invalid number of input ' num2str(nargin) ' (should be 17) or output ' num2str(nargout) ' (should the 15) parameters']);
    end
  end


  nb_connexity_threshold  = size(connexity_threshold_range, 2);
  nb_intensity_threshold  = size(intensity_threshold_range, 2);
  if (compute_all_features==1)
    nb_gray_threshold       = size(gray_threshold_range, 2);
    nb_white_threshold      = size(white_threshold_range, 2);
    nb_white_area_threshold = size(white_area_threshold_range, 2);
  end % compute all features

  % Initialization
  nb_well_index = size(index_wells_for_opti,2);

  distance_area_lcc_for_opti                        = cell(nb_well_index,1);
  distance_elliptical_regularity_lcc_for_opti       = cell(nb_well_index,1);
  distance_texture_psX_1eros_lcc_for_opti           = cell(nb_well_index,1);
  distance_intensity_ppminloc_lcc_for_opti          = cell(nb_well_index,1);
  distance_nb_minloc_lcc_for_opti                   = cell(nb_well_index,1);
  if (compute_all_features==1)
    distance_ratio_white_gray_for_opti                = cell(nb_well_index,1);
    distance_aspect_ratio_lcc_for_opti                = cell(nb_well_index,1);
    distance_nb_max_white_area_for_opti               = cell(nb_well_index,1);
    distance_nb_connected_comp_for_opti               = cell(nb_well_index,1);
    distance_concavity_lcc_for_opti                   = cell(nb_well_index,1);
    distance_texture_psX_2eros_lcc_for_opti           = cell(nb_well_index,1);
    distance_texture_psX_3lastlevel_eros_lcc_for_opti = cell(nb_well_index,1);
    distance_nb_maxloc_lcc_for_opti                   = cell(nb_well_index,1);
    distance_intensity_ppmaxloc_lcc_for_opti          = cell(nb_well_index,1);
  end % compute all features
  for index=1:nb_well_index
    distance_area_lcc_for_opti{index}                        = zeros(nb_connexity_threshold,1);
    distance_elliptical_regularity_lcc_for_opti{index}       = zeros(nb_connexity_threshold,1);
    distance_texture_psX_1eros_lcc_for_opti{index}           = zeros(nb_connexity_threshold,nb_intensity_threshold);
    distance_intensity_ppminloc_lcc_for_opti{index}          = zeros(nb_connexity_threshold,1);
    distance_nb_minloc_lcc_for_opti{index}                   = zeros(nb_connexity_threshold,1);
    if (compute_all_features==1)
      distance_ratio_white_gray_for_opti{index}                = zeros(nb_gray_threshold,nb_white_threshold);
      distance_aspect_ratio_lcc_for_opti{index}                = zeros(nb_connexity_threshold,1);
      distance_nb_max_white_area_for_opti{index}               = zeros(nb_connexity_threshold,nb_white_area_threshold);
      distance_nb_connected_comp_for_opti{index}               = zeros(nb_connexity_threshold,1);
      distance_concavity_lcc_for_opti{index}                   = zeros(nb_connexity_threshold,1);
      distance_texture_psX_2eros_lcc_for_opti{index}           = zeros(nb_connexity_threshold,nb_intensity_threshold);
      distance_texture_psX_3lastlevel_eros_lcc_for_opti{index} = zeros(nb_connexity_threshold,nb_intensity_threshold);
      distance_nb_maxloc_lcc_for_opti{index}                   = zeros(nb_connexity_threshold,1);
      distance_intensity_ppmaxloc_lcc_for_opti{index}          = zeros(nb_connexity_threshold,1);
    end % compute all features
  end

  % All the computations
  % 'normal' nuclei
  the_area_lcc_normal                        = [];
  the_elliptical_regularity_lcc_normal       = [];
  the_texture_psX_1eros_lcc_normal           = [];
  the_intensity_ppminloc_lcc_normal          = [];
  the_nb_minloc_lcc_normal                   = [];
  if (compute_all_features==1)
    the_ratio_white_gray_normal                = [];
    the_aspect_ratio_lcc_normal                = [];
    the_nb_max_white_area_normal               = [];
    the_nb_connected_comp_normal               = [];
    the_concavity_lcc_normal                   = [];
    the_texture_psX_2eros_lcc_normal           = [];
    the_texture_psX_3lastlevel_eros_lcc_normal = [];
    the_nb_maxloc_lcc_normal                   = [];
    the_intensity_ppmaxloc_lcc_normal          = [];
  end % compute all features

  for wellIndex=index_normal_wells_for_opti

    plateName = data(wellIndex).plate;
    wellName = data(wellIndex).well;

    disp(['normal : ' plateName wellName])

    for siteIndex = 1:nb_site_per_well

      siteName = ['s' num2str(siteIndex)];

      if (compute_all_features==1)
        load(['../results/features_for_parameters_optimization/' plateName '/res_' plateName wellName siteName '_features_for_parameters_optimization_all.mat']);
      else
        load(['../results/features_for_parameters_optimization/' plateName '/res_' plateName wellName siteName '_features_for_parameters_optimization.mat']);
      end

      the_area_lcc_normal                        = [the_area_lcc_normal area_lcc];
      the_elliptical_regularity_lcc_normal       = [the_elliptical_regularity_lcc_normal elliptical_regularity_lcc];
      the_texture_psX_1eros_lcc_normal           = cat(3, the_texture_psX_1eros_lcc_normal, texture_psX_1eros_lcc);
      the_intensity_ppminloc_lcc_normal          = [the_intensity_ppminloc_lcc_normal intensity_ppminloc_lcc];
      the_nb_minloc_lcc_normal                   = [the_nb_minloc_lcc_normal nb_minloc_lcc];
      if (compute_all_features==1)
        the_ratio_white_gray_normal                = cat(3, the_ratio_white_gray_normal, ratio_white_gray);
        the_aspect_ratio_lcc_normal                = [the_aspect_ratio_lcc_normal aspect_ratio_lcc];
        the_nb_max_white_area_normal               = cat(3, the_nb_max_white_area_normal, nb_max_white_area);
        the_nb_connected_comp_normal               = [the_nb_connected_comp_normal nb_connected_comp];
        the_concavity_lcc_normal                   = [the_concavity_lcc_normal concavity_lcc];
        the_texture_psX_2eros_lcc_normal           = cat(3, the_texture_psX_2eros_lcc_normal, texture_psX_2eros_lcc);
        the_texture_psX_3lastlevel_eros_lcc_normal = cat(3, the_texture_psX_3lastlevel_eros_lcc_normal, texture_psX_3lastlevel_eros_lcc);
        the_nb_maxloc_lcc_normal                   = [the_nb_maxloc_lcc_normal nb_maxloc_lcc]; 
        the_intensity_ppmaxloc_lcc_normal          = [the_intensity_ppmaxloc_lcc_normal intensity_ppmaxloc_lcc];
      end % compute all features
    end
  end % index normal wells for opti


  % the normal or abnormal wells
  index = 0;
  for wellIndex=index_wells_for_opti

    index = index + 1;


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

    plateName = data(wellIndex).plate;
    wellName = data(wellIndex).well;

    disp(['abnormal : ' plateName wellName])

    for siteIndex = 1:nb_site_per_well

      siteName = ['s' num2str(siteIndex)];
      if (compute_all_features==1)
        load(['../results/features_for_parameters_optimization/' plateName '/res_' plateName wellName siteName '_features_for_parameters_optimization_all.mat']);
      else
        load(['../results/features_for_parameters_optimization/' plateName '/res_' plateName wellName siteName '_features_for_parameters_optimization.mat']);
      end

      the_area_lcc_abnormal                        = [the_area_lcc_abnormal area_lcc];
      the_elliptical_regularity_lcc_abnormal       = [the_elliptical_regularity_lcc_abnormal elliptical_regularity_lcc];
      the_texture_psX_1eros_lcc_abnormal           = cat(3, the_texture_psX_1eros_lcc_abnormal, texture_psX_1eros_lcc);
      the_intensity_ppminloc_lcc_abnormal          = [the_intensity_ppminloc_lcc_abnormal intensity_ppminloc_lcc];
      the_nb_minloc_lcc_abnormal                   = [the_nb_minloc_lcc_abnormal nb_minloc_lcc];
      if (compute_all_features==1)
        the_ratio_white_gray_abnormal                = cat(3, the_ratio_white_gray_abnormal, ratio_white_gray);
        the_aspect_ratio_lcc_abnormal                = [the_aspect_ratio_lcc_abnormal aspect_ratio_lcc];
        the_nb_max_white_area_abnormal               = cat(3, the_nb_max_white_area_abnormal, nb_max_white_area);
        the_nb_connected_comp_abnormal               = [the_nb_connected_comp_abnormal nb_connected_comp];
        the_concavity_lcc_abnormal                   = [the_concavity_lcc_abnormal concavity_lcc];
        the_texture_psX_2eros_lcc_abnormal           = cat(3, the_texture_psX_2eros_lcc_abnormal, texture_psX_2eros_lcc);
        the_texture_psX_3lastlevel_eros_lcc_abnormal = cat(3, the_texture_psX_3lastlevel_eros_lcc_abnormal, texture_psX_3lastlevel_eros_lcc);
        the_nb_maxloc_lcc_abnormal                   = [the_nb_maxloc_lcc_abnormal nb_maxloc_lcc];
        the_intensity_ppmaxloc_lcc_abnormal          = [the_intensity_ppmaxloc_lcc_abnormal intensity_ppmaxloc_lcc];
      end % compute all features
    end

    % Compute distances...
    for ind_connexity_threshold=1:nb_connexity_threshold

      distance_area_lcc_for_opti{index}(ind_connexity_threshold)                        = f_get_fisher_criterion(the_area_lcc_normal(ind_connexity_threshold,:), the_area_lcc_abnormal(ind_connexity_threshold,:));
      distance_elliptical_regularity_lcc_for_opti{index}(ind_connexity_threshold)       = f_get_fisher_criterion(the_elliptical_regularity_lcc_normal(ind_connexity_threshold,:), the_elliptical_regularity_lcc_abnormal(ind_connexity_threshold,:));

      for ind_intensity_threshold=1:nb_intensity_threshold
        if (intensity_threshold_range(ind_intensity_threshold)>connexity_threshold_range(ind_connexity_threshold))
          normal   = the_texture_psX_1eros_lcc_normal(ind_connexity_threshold,ind_intensity_threshold,:);
          abnormal = the_texture_psX_1eros_lcc_abnormal(ind_connexity_threshold,ind_intensity_threshold,:);
          distance_texture_psX_1eros_lcc_for_opti{index}(ind_connexity_threshold, ind_intensity_threshold)           = f_get_fisher_criterion(normal(:)', abnormal(:)');
        end
      end % intensity threshold
      distance_intensity_ppminloc_lcc_for_opti{index}(ind_connexity_threshold)          = f_get_fisher_criterion(the_intensity_ppminloc_lcc_normal(ind_connexity_threshold,:), the_intensity_ppminloc_lcc_abnormal(ind_connexity_threshold,:));
      distance_nb_minloc_lcc_for_opti{index}(ind_connexity_threshold)                   = f_get_fisher_criterion(the_nb_minloc_lcc_normal(ind_connexity_threshold,:), the_nb_minloc_lcc_abnormal(ind_connexity_threshold,:));

      if (compute_all_features==1)
        distance_aspect_ratio_lcc_for_opti{index}(ind_connexity_threshold)                = f_get_fisher_criterion(the_aspect_ratio_lcc_normal(ind_connexity_threshold,:), the_aspect_ratio_lcc_abnormal(ind_connexity_threshold,:));

        for ind_white_area_threshold=1:nb_white_area_threshold
          normal   = the_nb_max_white_area_normal(ind_connexity_threshold, ind_white_area_threshold,:);
          abnormal = the_nb_max_white_area_abnormal(ind_connexity_threshold, ind_white_area_threshold,:);
          distance_nb_max_white_area_for_opti{index}(ind_connexity_threshold, ind_white_area_threshold) = f_get_fisher_criterion(normal(:)', abnormal(:)');
        end % intensity threshold

        distance_nb_connected_comp_for_opti{index}(ind_connexity_threshold)               = f_get_fisher_criterion(the_nb_connected_comp_normal(ind_connexity_threshold,:), the_nb_connected_comp_abnormal(ind_connexity_threshold,:));
        distance_concavity_lcc_for_opti{index}(ind_connexity_threshold)                   = f_get_fisher_criterion(the_concavity_lcc_normal(ind_connexity_threshold,:), the_concavity_lcc_abnormal(ind_connexity_threshold,:));


        for ind_intensity_threshold=1:nb_intensity_threshold
          if (intensity_threshold_range(ind_intensity_threshold)>connexity_threshold_range(ind_connexity_threshold))
            normal   = the_texture_psX_2eros_lcc_normal(ind_connexity_threshold,ind_intensity_threshold,:);
            abnormal = the_texture_psX_2eros_lcc_abnormal(ind_connexity_threshold,ind_intensity_threshold,:);
            distance_texture_psX_2eros_lcc_for_opti{index}(ind_connexity_threshold, ind_intensity_threshold)           = f_get_fisher_criterion(normal(:)', abnormal(:)');

            normal   = the_texture_psX_3lastlevel_eros_lcc_normal(ind_connexity_threshold,ind_intensity_threshold,:);
            abnormal = the_texture_psX_3lastlevel_eros_lcc_abnormal(ind_connexity_threshold,ind_intensity_threshold,:);
            distance_texture_psX_3lastlevel_eros_lcc_for_opti{index}(ind_connexity_threshold, ind_intensity_threshold) = f_get_fisher_criterion(normal(:)', abnormal(:)');
          end
        end % intensity threshold

        distance_nb_maxloc_lcc_for_opti{index}(ind_connexity_threshold)                   = f_get_fisher_criterion(the_nb_maxloc_lcc_normal(ind_connexity_threshold,:), the_nb_maxloc_lcc_abnormal(ind_connexity_threshold,:));
        distance_intensity_ppmaxloc_lcc_for_opti{index}(ind_connexity_threshold)          = f_get_fisher_criterion(the_intensity_ppmaxloc_lcc_normal(ind_connexity_threshold,:), the_intensity_ppmaxloc_lcc_abnormal(ind_connexity_threshold,:));
      end % compute all features
    end % connexity threshold

    if (compute_all_features==1)
      for ind_white_threshold=1:nb_white_threshold
        for ind_gray_threshold=1:nb_gray_threshold
  	  normal   = the_ratio_white_gray_normal(ind_gray_threshold,ind_white_threshold,:);
          abnormal = the_ratio_white_gray_abnormal(ind_gray_threshold,ind_white_threshold,:);
          distance_ratio_white_gray_for_opti{index}(ind_gray_threshold, ind_white_threshold)                = f_get_fisher_criterion(normal(:)', abnormal(:)');
        end % gray threshold
      end % white threshold
    end % compute all features

  end % index wells for opti


end
