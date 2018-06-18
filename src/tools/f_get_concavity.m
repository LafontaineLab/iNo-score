% Function that computes the concavity measure of each connected component
%
% concavity = 1 - connected_component_area / convex_hull_area
%
% Remark : concavity value in [0,1] : 0 => convex shape
%                                     1 => highly concave shape

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

function concavity = f_get_concavity(I_BW)

% Input:
% ******
% I_BW : black and white image with "n" white connected components
%
% Output: 
% *******
% concavity : concavity measure for each connected component (1xn vector)

  [L_BW, nb_labels] = bwlabel(I_BW, 8);
  stats = regionprops(L_BW, 'Area','PixelList');

  [ht lg] = size(I_BW);

  % List of pixels coordinates in the image
  [x_im,y_im] = meshgrid(1:lg, 1:ht);
  x_im = x_im(:);
  y_im = y_im(:);

  % Boundary of the connected components
  perim = bwperim(I_BW);

  concavity = zeros(1,nb_labels);

  for l = 1:nb_labels

    % Pixels on the boundary of the connected component
    [y x] = find(perim==1 & L_BW==l);

    if (size(y,1)>2) % At least 2 points

      % Are the points colinear ?
      co = [x-mean(x) y-mean(y)];
      co = co'*co/size(y,1);

      if (rank(co)==2) % All the points are not colinear
   
        % Convex hull
        k = convhull(x,y);

        % Pixels inside the convex hull
        in = inpolygon(x_im,y_im,x(k),y(k));

        % Concavity measure
        concavity(l) = 1 - stats(l).Area/sum(in);

      else % Points are colinear

        concavity(l) = 0;

      end
    else % less than 2 points

      concavite(l) = 0;

    end

  end % for connected components

end % function
