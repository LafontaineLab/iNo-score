% Script that computes the nuclei segmentation for all the images
% and save the BW mask in a file

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
dataconfig

% Creation of the subdirectories to save the results
if ~exist(['../results/bw'])
  mkdir('../results/', 'bw')
end
for plateNameIndex=1:nb_wells
  if ~exist(['../results/bw/' DATA(plateNameIndex).plate])
    mkdir('../results/bw/', DATA(plateNameIndex).plate)
  end
end

% Initialization of the thresholds
small_thres = 500; % size of the smallest nucleus (in pixel)
big_thres = 5000;  % size of the biggest nucleus (in pixel)
the_binarization_thres=[250:20:450 500:100:1300]; % the tested thresholds 
                                                  % for binarization/segmentation

for wellIndex=1:nb_wells

  plateName = DATA(wellIndex).plate;
  wellName = DATA(wellIndex).well;

  for siteIndex = 1:nb_site_per_well

    siteName = ['s' num2str(siteIndex)];

    disp(['Compute segmentation: ' plateName ' ' wellName ' ' siteName]);

    I_nuclei = f_read_image(plateName, wellName, siteName);
tic
    BW = f_get_nucleus_segmentation(I_nuclei, small_thres, ...
                                    big_thres, the_binarization_thres);
toc
    imwrite(BW, ['../results/bw/' plateName '/' plateName '_' wellName '_' siteName '_bw.png']);

  end % for site indexes

end % for well indexes
