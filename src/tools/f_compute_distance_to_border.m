% Function that associates to each pixel of a connected component, 
% its distance to the boundary of the connected component.

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

function I_distance = f_compute_distance_to_border(I_BW_comp)

% Input:
% ******
% I_BW_comp  : black and white image containing in white the connected 
%             components relative to the nucleoli
%
% Output:
% *******
% I_distance : image containing for each pixel its ditance to the boundary 
%              of its connected component

  I_distance = zeros(size(I_BW_comp));

  element = ones(3,3);
  I_BW_comp_new_prec = I_BW_comp;
  i = 1;
  while sum(sum(I_BW_comp_new_prec))>0
    I_BW_comp_new = imerode(I_BW_comp_new_prec,element);
    I_distance(I_BW_comp_new==0&I_BW_comp_new_prec==1) = i;
    I_BW_comp_new_prec = I_BW_comp_new;
    i = i + 1;
  end

end % function
