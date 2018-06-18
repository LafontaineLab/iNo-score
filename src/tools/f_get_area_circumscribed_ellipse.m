% Function that determines and computes, for each connected component, 
% the area of the circumscribed ellipse of this connected component.

% We call "circumscribed ellipse" the smallest ellipse such that all 
% the pixels of the the connected component is included in that ellipse.

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

function [area_pixel, area_analytical] = f_get_area_circumscribed_ellipse(I_BW, I, type_ellipse)

% Inputs:
% *******
% I_BW         : black and white image
% I            : image
% type_ellipse : 1 => we are interested by the shape of the connected component
%                0 => we are interested by the intensity distribution 
%                     of the connected component
%
% Outputs:
% *******
% area_pixel      : area of the ellipse determine by the number of pixels
%                   inside the ellipse
% area_analytical : area of the ellipse determined by the analytical equation

  [ht, lg] = size(I_BW);

  [L_BW, nb_labels] = bwlabel(I_BW, 8);
  stats = regionprops(L_BW, I, 'PixelList', 'PixelValues');

  % Boundary of the connected components
  perim = bwperim(I_BW);

  % Pixel area of the circumscribed ellipse
  % (count the number of pixels inside the ellipse)
  area_pixel = zeros(1,nb_labels);

  % Analytical area of the circumscribed ellipse: A = pi*a*b 
  % (where a and b are the half-axes of the ellipse)
  area_analytical = zeros(1,nb_labels);

  % List of pixel coordinates inside the image
  [x_im,y_im] = meshgrid(1:lg, 1:ht);
  x_im = x_im(:);
  y_im = y_im(:);

  for l = 1:nb_labels

    % Pixels coordinates of the connected component
    x_cc = stats(l).PixelList(:,1);
    y_cc = stats(l).PixelList(:,2);

    % Pixels coordinates of the boundary of the connected component
    [y_contour x_contour] = find(L_BW==l & perim==1);

    if (type_ellipse==1)

      % Principal component analysis of the points of the connected component
      % => principal axes of the connected component
      % => length of the half-axes (cf. eigenvalues)
      [changement_repere, xy_cc_nv_repere, lambda] = pca([x_cc y_cc]);

      origine = [mean(x_cc) mean(y_cc)];

    else

      % Analysis of the points intensity of the connected component

      %% moment of ordre 2 of the connected component
      pts = stats(l).PixelList;
      nb_pts = size(pts,1);

      intensity = stats(l).PixelValues;
      normalisation = sum(intensity);

      % gravity center
      gravity_center = intensity'*pts/normalisation;

      pts_center = pts - ones(nb_pts,1)*gravity_center;

      % Calculate normalized second central moments for the region. 1/12 is
      % the normalized second central moment of a pixel with unit length.
      uxx = intensity'*(pts_center(:,1).^2)/normalisation + 1/12;
      uyy = intensity'*(pts_center(:,2).^2)/normalisation + 1/12;
      uxy = intensity'*(pts_center(:,1).*pts_center(:,2))/normalisation;

      A = [ uxx uxy
            uxy uyy];

      [changement_repere, lambda] = eig(A);
      lambda = diag(lambda);

      origine = gravity_center;

    end

    % Lengths of the principal axes of the ellipse
    a = sqrt(lambda(1));
    b = sqrt(lambda(2));

    if (a>0 & b>0) % real ellipse

      % Project the image pixels in the new basis
      xy_im_nv = [x_im-origine(1) y_im-origine(2)]*changement_repere;

      % Project the boundary pixels in the new basis
      xy_contour_nv = [x_contour-origine(1) y_contour-origine(2)]*changement_repere;
  
      % For each boundary point of the ellipse, compute (x/a)^2 + (y/b)^2
      dist = sqrt((xy_contour_nv(:,1)/a).^2 + (xy_contour_nv(:,2)/b).^2);

      % Half_axes of the circumscribed ellipse
      ra = max(dist)*a;
      rb = max(dist)*b;

      % Number of pixels in the ellipse of radius ra and rb
      dist_in = sqrt((xy_im_nv(:,1)/ra).^2 + (xy_im_nv(:,2)/rb).^2);

      area_analytical(l) = pi*ra*rb;

      area_pixel(l) = sum(dist_in<=1);

    else

      area_analytical(l) = size(x_cc,1);

      area_pixel(l) = size(x_cc,1);

    end

  end % for connected component

end % function
