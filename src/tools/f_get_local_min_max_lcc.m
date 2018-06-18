% Function that compute the local maxima and minima on the 
% biggest connected component of the nucleoli signal of a nucleus

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

function [ xy_min, xy_max, nb_min, nb_max, ...
           intensity_ppmin, intensity_ppmax ...
         ] = f_get_local_min_max_lcc(I, I_normalized, I_segmented, ...
                                     I_erode, nb_erosion)

% Inputs:
% *******
% I            : nucleolus image for the local maxima and minima search
% I_normalized : normalized image of the nucleoli
% I_segmented  : image of the segmented nucleoli
% I_erode      : image that contains for each pixel of a segmented area 
%                its index of erosion
% nb_erosion   : number of erosions to remove the local maxima and minima 
%                that are close to the boundary of a connected component
%
% Outputs:
% ********
% xy_min          : coordinates of the local minima in the biggest 
%                   connected component
% xy_max          : coordinates of the local maxima in the biggest 
%                   connected component
% nb_min          : number of local minima in the biggest connected component
% nb_max          : number of local maxima in the biggest connected component
% intensity_ppmin : intensity of the smallest local minima in the biggest 
%                   connected component (=255 if no local minimum)
% intensity_ppmax : intensity of the smallest local maxima in the biggest 
%                   connected component (=0 if no local maximum)

  I_BW = I>0;
  I_BW = imerode(I_BW,strel('disk',2));

  [L_segmented_BW, nb_labels] = bwlabel(I_segmented, 8);
  stats = regionprops(L_segmented_BW, 'Area');

  I_erode_BW = I_erode>nb_erosion; % removed nb_erosions erosions

  [ht, lg] = size(I);

  I_vd = zeros(ht,lg);
  I_vd(:,1:end-1) = sign(I(:,1:end-1)-I(:,2:end));
  I_vg = zeros(ht,lg);
  I_vg(:,2:end) = sign(I(:,2:end)-I(:,1:end-1));
  I_vh = zeros(ht,lg);
  I_vh(2:end,:) = sign(I(2:end,:)-I(1:end-1,:));
  I_vb = zeros(ht,lg);
  I_vb(1:end-1,:) = sign(I(1:end-1,:)-I(2:end,:));

  I_vhd = zeros(ht,lg);
  I_vhd(2:end,1:end-1) = sign(I(2:end,1:end-1)-I(1:end-1,2:end));
  I_vbd = zeros(ht,lg);
  I_vbd(1:end-1,1:end-1) = sign(I(1:end-1,1:end-1)-I(2:end,2:end));

  I_vhg = zeros(ht,lg);
  I_vhg(2:end,2:end) = sign(I(2:end,2:end)-I(1:end-1,1:end-1));
  I_vbg = zeros(ht,lg);
  I_vbg(1:end-1,2:end) = sign(I(1:end-1,2:end)-I(2:end,1:end-1));

  sI = I_vd+I_vg+I_vh+I_vb+I_vhd+I_vhg+I_vbg+I_vbd;

  [marea lab] = max(cat(1,stats.Area));

  % local maxima
  [y x] = find(sI==8 & L_segmented_BW==lab & I_erode_BW==1);

  if isempty(x)
    nb_max = 0;
    xy_max = [];
    intensity_ppmax = 0;
  else
    nb_max = size(x,1);
    xy_max = [x y];
    intensity_ppmax = min(I_normalized(xy_max(:,2)+(xy_max(:,1)-1)*ht));
  end

  % local minima
  [y x] = find(sI==-8 & L_segmented_BW==lab & I_erode_BW==1);

  if isempty(x)
    nb_min = 0;
    xy_min = [];
    intensity_ppmin = 255;
  else
    nb_min = size(x,1);
    xy_min = [x y];
    intensity_ppmin = min(I_normalized(xy_min(:,2)+(xy_min(:,1)-1)*ht));
  end

end % function
