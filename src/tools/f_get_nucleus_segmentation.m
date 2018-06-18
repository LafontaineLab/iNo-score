% Function that segments the nuclei image with an adaptative threshold

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

function BW = f_get_nucleus_segmentation(I_nuclei, ...
                                         small_thres, big_thres, ...
                                         the_bin_thres)

% Inputs:
% *******
% I_nuclei      : nuclei image
% small_thres   : size of the smallest nucleus (in pixel)
% big_thres     : size of the biggest nucleus (in pixel)
% the_bin_thres : the vector of binarization thresholds (nx1 vector)
%
% Output:
% *******
% BW          : nuclei segmented image (1: nucleus, 0: background)

  [ht, lg] = size(I_nuclei);

  BW_final = zeros(ht, lg);
  deb = 1;

  for binarization_thres = the_bin_thres

    I_BW = I_nuclei>binarization_thres;
    I_BW = imclose(I_BW,ones(3,3));
    BW = I_BW;
    stat = regionprops(I_BW,'FilledImage','BoundingBox');

    % Holes filling
    for i =1:size(stat,1)
      top = stat(i).BoundingBox(2)+0.5;
      bottom = stat(i).BoundingBox(2)+0.5+stat(i).BoundingBox(4)-1;
      left = stat(i).BoundingBox(1)+0.5;
      right = stat(i).BoundingBox(1)+0.5+stat(i).BoundingBox(3)-1;
      BW(top:bottom,left:right) =  BW(top:bottom,left:right) ...
                                 | stat(i).FilledImage;
    end

    [L, nblabels] = bwlabel(BW, 8);
    stat = regionprops(L, I_nuclei, 'Area', 'MeanIntensity');
    area = cat(1,stat.Area);
    mean_intensity = cat(1,stat.MeanIntensity); 

    if (deb==0)

      for i=1:nb_labels_prec
        % Check the size of the connected components inside 
        % a previous big connected component
        ind_cell = setdiff(unique(L(L_prec==i)),0);
        nb_big_cell = sum(area(ind_cell)>small_thres);

        % If a big connected component is divided into at least two big 
        % enough connected components, keep those connected components
        % Else if the big connected component reminds big, keep the 
        % previous segmentation except if there is a relevent improvement
        % in the segmentation
        % Else if the big connected component is divided into one big 
        % connected and small connected components, keep the previous 
        % segmentation.

        if (nb_big_cell>1)
          BW_final(L_prec==i) = 0;
          BW_final(L_prec==i & L>0) = 1;
        elseif (nb_big_cell==1 & size(ind_cell,1)==1 & area(ind_cell)>1400)
          if (mean_intensity(ind_cell)/mean_intensity_prec(i)>1.15)
            BW_final(L_prec==i) = 0;
            BW_final(L_prec==i & L>0) = 1;
          end
        end
      end % for previous connected components

    else
      BW_final = BW;
      deb=0;
    end

    I_BW = BW;

    [L, nb_labels] = bwlabel(BW_final, 8);
    % Remove too small connected components
    stat = regionprops(L, 'Area', 'PixelIdxList');
    for i = 1:nb_labels
      if (stat(i).Area<small_thres)
        BW_final(stat(i).PixelIdxList) = 0;
      end
    end

    [L_prec, nb_labels_prec] = bwlabel(BW_final, 8);

    stat = regionprops(L_prec, I_nuclei, 'MeanIntensity');
    mean_intensity_prec = cat(1,stat.MeanIntensity);

  end % for binarization threshold

  BW = BW_final;

  [L, nblabels] = bwlabel(BW, 8);

  % Remove the connected component that touch the image boundary
  % => incomplete nuclei
  ind_label = unique([L(:,1); L(:,end); L(1,:)'; L(end,:)']);

  for i=ind_label'
    BW(L==i) = 0;
  end
  [L, nblabels] = bwlabel(BW, 8);

  % Remove too big nucleus
  stat = regionprops(L,'Area','PixelIdxList');
  for i = 1:size(stat,1)
    if stat(i).Area>big_thres
      BW(stat(i).PixelIdxList) = 0;
    end
  end

end % function
