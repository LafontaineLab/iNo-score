% Function that computes for each valid nuclei of the image a set of features
% for a range of parameter threshold
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

function [ nb_nuclei_valid, area_lcc, ...
           elliptical_regularity_lcc, texture_psX_1eros_lcc, ...
           intensity_ppminloc_lcc, nb_minloc_lcc, ...
           ratio_white_gray, aspect_ratio_lcc, ...
           nb_max_white_area, nb_connected_comp, concavity_lcc, ...
           texture_psX_2eros_lcc, texture_psX_3lastlevel_eros_lcc, ...
           nb_maxloc_lcc, intensity_ppmaxloc_lcc ...
         ] = f_get_nucleus_features_for_parameters_optimization(compute_all_features, ...
                                                                name1, name2, name3, ...
                                                                percentile_for_shape, ...
                                                                connexity_threshold_range, ...
								intensity_threshold_range, ...
								gray_threshold_range, ...
								white_threshold_range, ...
								white_area_threshold_range)

% Inputs:
% *******
% compute_all_features              : if 0 => compute the 5 features
%                                     if 1 => compute all the features
% name1                             : the first part of the image name 
%                                     (plate name)
% name2                             : the second part of the image name  
%                                     (well name)
% name3                             : the third part of the image name  
%                                     (image number)
% percentile_for_shape              : percentile value (value between 0 and 100)
%                                     of the nucleolus intensity normalization 
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
% nb_nuclei_valid                   : number "n" of valid nuclei in the image
% area_lcc                          : for each nucleus, size of the largest 
%                                     connected component (c x n matrix)
% elliptical_regularity_lcc         : for each nucleus, elliptical regularity 
%                                     measure of the largest connected
%                                     component (c x n matrix)
% texture_psX_1eros_lcc             : for each nucleus, percent of pixels with 
%                                     an intensity smaller than X after one 
%                                     erosion of the largest connected component
%                                     (c x i x n matrix)
% intensity_ppminloc_lcc            : for each nucleus, smallest intensity 
%                                     of the local minima of the largest 
%                                     connected component (c x n matrix)
% nb_minloc_lcc                     : for each nucleus, number of local minima 
%                                     of the largest connected component 
%                                     (c x n matrix)
%
%% if compute all the features the following parameters are computed
%
% ratio_white_gray                  : for each nucleus, the white/gray ratio
%                                     = number of pixels with an intensity 
%                                       above a threshold in the connected 
%                                       components (g x w x n matrix)
% aspect_ratio_lcc                  : for each nucleus, aspect ration measure 
%                                     of the largest connected component
%                                     (c x n matrix)
% nb_max_white_area                 : for each nucleus, maximal number of 
%                                     white areas in a connected component 
%                                     (c x wa x n matrix)
% nb_connected_comp                 : for each nucleus, number of nucleolus 
%                                     connected components (c x n matrix)
% concavity_lcc                     : for each nucleus, concavity measure 
%                                     of the largest connected component
%                                     (c x n matrix)
% texture_psX_2eros_lcc             : for each nucleus, percent of pixels with 
%                                     an intensity smaller than X after two
%                                     erosions of the largest connected
%                                     component (c x i x n matrix) 
% texture_psX_3lastlevel_eros_lcc   : for each nucleus, percent of pixels with 
%                                     an intensity smaller than X on the last 
%                                     three levels of erosions of the largest 
%                                     connected component (c x i x n matrix)
% nb_maxloc_lcc                     : for each nucleus, number of local maxima 
%                                     of the largest connected component 
%                                     (c x n matrix)
% intensity_ppmaxloc_lcc            : for each nucleus, smallest intensity 
%                                     of the local maxima of the largest 
%                                     connected component (c x n matrix)


  % Check the input and output parameters
  if (compute_all_features==0) % compute the 5 features
    if ((nargin~=7) || (nargout~=6))
      error(['invalid number of input ' num2str(nargin) ' (should be 10) or output ' num2str(nargout) ' (should the 6) parameters']);
    end
  else % compute all the features
    if ((nargin~=10) || (nargout~=15))
      error(['invalid number of input ' num2str(nargin) ' (should be 17) or output ' num2str(nargout) ' (should the 15) parameters']);
    end
  end


  nb_connexity_threshold = size(connexity_threshold_range, 2);
  nb_intensity_threshold = size(intensity_threshold_range, 2);

  if (compute_all_features==1)
    nb_gray_threshold = size(gray_threshold_range, 2);
    nb_white_threshold = size(white_threshold_range, 2);
    nb_white_area_threshold = size(white_area_threshold_range, 2);
  end % compute all features

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Labelling of segmented nuclei
  [I_nuclei, I_nucleoli, BW] = f_read_image(name1, name2, name3);

  [ht, lg] = size(I_nucleoli);

  [L, nblabels] = bwlabel(BW, 8);

  perimeter = bwperim(BW);
  [y x] = find(perimeter==1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute information on the nuclei

  stats_nucleus = regionprops(L, I_nucleoli, 'Area', 'BoundingBox', ...
		    'PixelIdxList', 'PixelList', 'MaxIntensity', ...
		    'MeanIntensity', 'MinIntensity', 'PixelValues', ...
		    'WeightedCentroid', 'MinorAxisLength', 'MajorAxisLength');

  area_lcc                        = zeros(nb_connexity_threshold, nblabels);
  elliptical_regularity_lcc       = zeros(nb_connexity_threshold, nblabels);
  texture_psX_1eros_lcc           = zeros(nb_connexity_threshold, nb_intensity_threshold, nblabels);
  intensity_ppminloc_lcc          = zeros(nb_connexity_threshold, nblabels);
  nb_minloc_lcc                   = zeros(nb_connexity_threshold, nblabels);

  if (compute_all_features==1)
    ratio_white_gray                = zeros(nb_gray_threshold, nb_white_threshold, nblabels);
    aspect_ratio_lcc                = zeros(nb_connexity_threshold, nblabels);
    nb_max_white_area               = zeros(nb_connexity_threshold, nb_white_area_threshold, nblabels);
    nb_connected_comp               = zeros(nb_connexity_threshold, nblabels);
    concavity_lcc                   = zeros(nb_connexity_threshold, nblabels);
    texture_psX_2eros_lcc           = zeros(nb_connexity_threshold, nb_intensity_threshold, nblabels);
    texture_psX_3lastlevel_eros_lcc = zeros(nb_connexity_threshold, nb_intensity_threshold, nblabels);
    nb_maxloc_lcc                   = zeros(nb_connexity_threshold, nblabels);
    intensity_ppmaxloc_lcc          = zeros(nb_connexity_threshold, nblabels);
  end

  nb_nuclei_valid = 0;

  for label_nucleus = 1:nblabels

    disp([name1 name2 name3 ' : nucleus ' num2str(label_nucleus)])

    bb = floor(stats_nucleus(label_nucleus).BoundingBox);
    idx_pts = stats_nucleus(label_nucleus).PixelIdxList;
    pts = stats_nucleus(label_nucleus).PixelList;
    nb_pts = size(pts,1);
    pts_nucleus = pts-ones(nb_pts,1)*[bb(1) bb(2)]+ones(nb_pts,1)*[10 10];

    I_nucleus = zeros(bb(4)+20,bb(3)+20);
    I_nucleus(pts_nucleus(:,2)+(pts_nucleus(:,1)-1)*(bb(4)+20)) = I_nuclei(idx_pts);

    I_nucleo = zeros(bb(4)+20,bb(3)+20);
    I_nucleo(pts_nucleus(:,2)+(pts_nucleus(:,1)-1)*(bb(4)+20)) = I_nucleoli(idx_pts);

    I_BW = zeros(bb(4)+20,bb(3)+20);
    I_BW(pts_nucleus(:,2)+(pts_nucleus(:,1)-1)*(bb(4)+20)) = BW(idx_pts);

    type_nucleus = f_get_type_nucleus(I_nucleus, I_nucleo, I_BW);

    % If it is a simple nucleus, i.e. no multiple nucleus, nor strange one...
    if (type_nucleus==1)

      nb_nuclei_valid = nb_nuclei_valid + 1;

      % The nucleus is enlarged to be sure to have the totality of the nucleus
      I_BW = imdilate(I_BW, ones(13,13));
      I_BW = imerode(I_BW, ones(3,3));

      % pixel inside a connected component
      [y_d x_d] = find(I_BW==1);
      y_d_im = y_d-10+bb(2);
      x_d_im = x_d-10+bb(1);
      i_in = find(  (1<=y_d_im) & (y_d_im<=ht)   ...
                  & (1<=x_d_im) & (x_d_im<=lg) );

      y_d_loc = y_d(i_in);
      x_d_loc = x_d(i_in);
      ind_loc = y_d_loc+(x_d_loc-1)*(bb(4)+20);
      y_d_im  = y_d_im(i_in);
      x_d_im  = x_d_im(i_in);
      ind_im  = y_d_im+(x_d_im-1)*ht;

      I_BW = zeros(bb(4)+20,bb(3)+20);
      I_BW(ind_loc) = 1;
      I_nucleus(ind_loc) = I_nuclei(ind_im);
      I_nucleo(ind_loc) = I_nucleoli(ind_im);

      %%%%%%%%%%%%%%%%%%%%%%%%%%
      % Normalization of the nuclei image (nucleus by nucleus)
      % The intensity is scaled between 1 and 255 on the basis of the 
      % percentile for shape
      % percentile_for_shape = 99.9

      I_nucleo_normalized = I_nucleo;
      intensity = I_nucleo_normalized(ind_loc);
      intensity_percentile = prctile(intensity,percentile_for_shape);
      min_intensity = min(intensity);
      a = (255-1)/(intensity_percentile-min_intensity);
      b = (intensity_percentile-255*min_intensity)/(intensity_percentile-min_intensity);
      nv_intensity = a*intensity + b;
      nv_intensity(nv_intensity>255) = 255;
      I_nucleo_normalized(ind_loc) = nv_intensity;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Compute all the feature of one nucleus
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ind_connexity_threshold = 0;

      for connexity_thres=connexity_threshold_range

        ind_connexity_threshold = ind_connexity_threshold + 1;

        %%% AREA LCC
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cut the nuclei in 2 gray levels
        % 0 => background
        % 1 => high enough intensity => connexity (connected component)
        I_nucleo_connexity = I_nucleo_normalized>connexity_thres;

        I_BW_connexity = I_nucleo_connexity>0;
        [L_nucleus_connexity, nblabels_nucleus_connexity] = bwlabel(I_BW_connexity, 8);

        stats_nucleus_connexity = regionprops(L_nucleus_connexity, 'Area');

        % index of the largest connected component
        [area_max, ind_max_size] = max(cat(1,stats_nucleus_connexity.Area));

        %% area of the largest connected component
        area_lcc(ind_connexity_threshold, nb_nuclei_valid) = area_max;

        %%% ASPECT RATIO, ELLIPTICAL REGULARITY, CONCAVITY LCC
        I_nucleo_connexity = I_nucleo_normalized>connexity_thres;

        I_BW_connexity = I_nucleo_connexity>0;
        [L_nucleus_connexity, nblabels_nucleus_connexity] = bwlabel(I_BW_connexity, 8);

        stats_nucleus_connexity = regionprops(L_nucleus_connexity, 'MajorAxisLength', 'MinorAxisLength', 'Area');

        % index of the largest connected component (lcc)
        [area_max, ind_max_size] = max(cat(1,stats_nucleus_connexity.Area));

        if (compute_all_features==1)

          %% aspect ratio of the largest connected component
          aspect_ratio_lcc(ind_connexity_threshold, nb_nuclei_valid) = stats_nucleus_connexity(ind_max_size).MinorAxisLength/stats_nucleus_connexity(ind_max_size).MajorAxisLength;

          %% concavity of the largest connected component
          concavity = f_get_concavity(I_BW_connexity);
          concavity_lcc(ind_connexity_threshold, nb_nuclei_valid) = concavity(ind_max_size);

        end % compute all features

        %% elliptical regularity of the largest connected component
        [area_inertie_pixel, area_inertie_analytique] = f_get_area_circumscribed_ellipse(I_BW_connexity, I_nucleo, 0);
        elliptical_regularity_lcc(ind_connexity_threshold, nb_nuclei_valid) = stats_nucleus_connexity(ind_max_size).Area/area_inertie_pixel(ind_max_size);

        if (compute_all_features==1)

          %%% NUMBER OF CONNECTED COMPONENTS
          I_nucleo_connexity = I_nucleo_normalized>connexity_thres;
          I_BW_connexity = I_nucleo_connexity>0;
          [L_nucleus_connexity, nblabels_nucleus_connexity] = bwlabel(I_BW_connexity, 8);

          %% number of connected components (= of nucleoli) of the nucleus
          nb_connected_comp(ind_connexity_threshold, nb_nuclei_valid) = nblabels_nucleus_connexity;


          %%% MAXIMAL NUMBER OF WHITE AREA PER NUCLEUS
          I_nucleo_connexity = I_nucleo_normalized>connexity_thres;

          I_BW_connexity = I_nucleo_connexity>0;
          [L_nucleus_connexity, nblabels_nucleus_connexity] = bwlabel(I_BW_connexity, 8);

          stats_nucleus_connexity = regionprops(L_nucleus_connexity, 'PixelIdxList');

          ind_white_area_threshold = 0;
          for white_area_thres=white_area_threshold_range
            ind_white_area_threshold = ind_white_area_threshold + 1;

            if (white_area_thres>=connexity_thres)
              %% Maximal number of white area per nucleus
              I_BW_peak = I_nucleo_normalized>white_area_thres;
              [L_nucleus_peak, nblabels_nucleus_peak] = bwlabel(I_BW_peak, 8);

              % Number of white peaks by connexity (connected component)
              nb_peak_per_comp = zeros(1,nblabels_nucleus_connexity);
              for lab=1:nblabels_nucleus_connexity
                ll = setdiff(unique(L_nucleus_peak(stats_nucleus_connexity(lab).PixelIdxList)),0);
                if isempty(ll)
                  nb_peak_per_comp(lab) = 0;
                else
                  nb_peak_per_comp(lab) = size(ll,1);
                end
              end

              nb_max_white_area(ind_connexity_threshold, ind_white_area_threshold, nb_nuclei_valid) = max(nb_peak_per_comp);
            end
          end % for white area threshold
        end % compute all features

        %%% LOCAL MINMAX
        I_nucleo_connexity = I_nucleo_normalized>connexity_thres;

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the distance of the pixels of a connected component 
        % to its boundary
        I_nucleo_connexity_erode_niv = f_compute_distance_to_border(I_nucleo_connexity);

        %% local minmax
        [xy_min, xy_max, nb_min, nb_max, intensity_ppmin, intensity_ppmax] = ...
          f_get_local_min_max_lcc(I_nucleo, I_nucleo_normalized, ...
                                  I_nucleo_connexity, ...
                                  I_nucleo_connexity_erode_niv, 1);

        intensity_ppminloc_lcc(ind_connexity_threshold, nb_nuclei_valid) = intensity_ppmin;
        nb_minloc_lcc(ind_connexity_threshold, nb_nuclei_valid) = nb_min;

        if (compute_all_features==1)
          intensity_ppmaxloc_lcc(ind_connexity_threshold, nb_nuclei_valid) = intensity_ppmax;
          nb_maxloc_lcc(ind_connexity_threshold, nb_nuclei_valid) = nb_max;
        end % compute all features

        %%% TEXTURE PSX
        I_nucleo_connexity = I_nucleo_normalized>connexity_thres;

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the distance of the pixels of a connected component 
        % to its boundary
        I_nucleo_connexity_erode_niv = f_compute_distance_to_border(I_nucleo_connexity);

        I_BW_connexity = I_nucleo_connexity>0;
        [L_nucleus_connexity, nblabels_nucleus_connexity] = bwlabel(I_BW_connexity, 8);

        stats_nucleus_connexity = regionprops(L_nucleus_connexity, 'Area');

        % index of the largest connected component
        [area_max, ind_max_size] = max(cat(1,stats_nucleus_connexity.Area));

        %% texture percentage under X erosions lcc (largest connected component)
        intensity_pixel_1eros           = sort(I_nucleo_normalized(I_nucleo_connexity_erode_niv>1  & L_nucleus_connexity==ind_max_size)','ascend');
        nb_intensity_pixel_1eros           = size(intensity_pixel_1eros,2);

        if (compute_all_features==1)
          intensity_pixel_2eros           = sort(I_nucleo_normalized(I_nucleo_connexity_erode_niv>2  & L_nucleus_connexity==ind_max_size)','ascend');
          nb_intensity_pixel_2eros           = size(intensity_pixel_2eros,2);

          intensity_pixel_3lastlevel_eros = sort(I_nucleo_normalized(I_nucleo_connexity_erode_niv>max(0,max(I_nucleo_connexity_erode_niv(L_nucleus_connexity==ind_max_size))-3) & L_nucleus_connexity==ind_max_size)','ascend');
          nb_intensity_pixel_3lastlevel_eros = size(intensity_pixel_3lastlevel_eros,2);
        end % compute all features

        ind_intensity_threshold = 0;
        for intensity_thres=intensity_threshold_range
          ind_intensity_threshold = ind_intensity_threshold + 1;
          if (intensity_thres>=connexity_thres)
            texture_psX_1eros_lcc(ind_connexity_threshold, ind_intensity_threshold, nb_nuclei_valid)           = sum(intensity_pixel_1eros<=intensity_thres)/max(1,nb_intensity_pixel_1eros);
            if (compute_all_features==1)
              texture_psX_2eros_lcc(ind_connexity_threshold, ind_intensity_threshold, nb_nuclei_valid)           = sum(intensity_pixel_2eros<=intensity_thres)/max(1,nb_intensity_pixel_2eros); 
              texture_psX_3lastlevel_eros_lcc(ind_connexity_threshold, ind_intensity_threshold, nb_nuclei_valid) = sum(intensity_pixel_3lastlevel_eros<=intensity_thres)/max(1,nb_intensity_pixel_3lastlevel_eros);
            end % compute all features
          end
        end % for intensity threshold

      end % for connexity threshold

      if (compute_all_features==1)
        %%% RATIO WHITE GRAY
        %% ratio white gray
        ind_gray_threshold = 0;
        for gray_thres=gray_threshold_range

          ind_gray_threshold = ind_gray_threshold + 1;
          ind_white_threshold = 0;

          for white_thres=white_threshold_range
            ind_white_threshold = ind_white_threshold + 1;
            if (white_thres>=gray_thres)
              ratio_white_gray(ind_gray_threshold, ind_white_threshold, nb_nuclei_valid) = sum(sum(I_nucleo_normalized>white_thres))/sum(sum(I_nucleo_normalized>gray_thres));
            end
          end % for white threshold
        end % for gray threshold
      end % compute all features

    end % if type_nucleus==normal

  end % for all segmented nuclei

  % in case there were non valid nuclei
  area_lcc                        = area_lcc(:,1:nb_nuclei_valid);
  elliptical_regularity_lcc       = elliptical_regularity_lcc(:,1:nb_nuclei_valid);
  texture_psX_1eros_lcc           = texture_psX_1eros_lcc(:,:,1:nb_nuclei_valid);
  intensity_ppminloc_lcc          = intensity_ppminloc_lcc(:,1:nb_nuclei_valid);
  nb_minloc_lcc                   = nb_minloc_lcc(:,1:nb_nuclei_valid);

  if (compute_all_features==1)
    ratio_white_gray                = ratio_white_gray(:,:,1:nb_nuclei_valid);
    aspect_ratio_lcc                = aspect_ratio_lcc(:,1:nb_nuclei_valid);
    nb_max_white_area               = nb_max_white_area(:,:,1:nb_nuclei_valid);
    nb_connected_comp               = nb_connected_comp(:,1:nb_nuclei_valid);
    concavity_lcc                   = concavity_lcc(:,1:nb_nuclei_valid);
    texture_psX_2eros_lcc           = texture_psX_2eros_lcc(:,:,1:nb_nuclei_valid);
    texture_psX_3lastlevel_eros_lcc = texture_psX_3lastlevel_eros_lcc(:,:,1:nb_nuclei_valid);
    nb_maxloc_lcc                   = nb_maxloc_lcc(:,1:nb_nuclei_valid);
    intensity_ppmaxloc_lcc          = intensity_ppmaxloc_lcc(:,1:nb_nuclei_valid);
  end

end % function
