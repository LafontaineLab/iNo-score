% Function that reads the nuclei, nucleoli and segmentation 
% images from the files

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

function [I_nuclei, I_nucleoli, I_BW] = f_read_image(name1, name2, name3)

% Inputs:
% *******
% name1 : the first part of the image name (plate name)
% name2 : the second part of the image name (well name)
% name3 : the third part of the image name (image number)
%
% Outputs:
% ********
% I_nuclei   : the nuclei image
% I_nucleoli : the nucleoli image
% I_BW       : the segmentation image

  I_nuclei = double(imread(['../data/' name1 '/' name1 '_' name2 '_' name3 '_dapi.TIF']));

  if (nargout>1)
    I_nucleoli = double(imread(['../data/' name1 '/' name1 '_' name2 '_' name3 '_gfp.TIF']));
  end

  if (nargout>2)
    I_BW = double(imread(['../results/bw/' name1 '/' name1 '_' name2 '_' name3 '_bw.png']));
  end

end % function
