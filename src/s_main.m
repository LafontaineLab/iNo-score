% Main script that launch all the scripts in order to analyze
% the degree of disruption of the nucleolus.

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

% Compute the nuclei segmentation
s_compute_segmentation

% Compute all the features on all the nuclei
s_compute_all_features

% Compute all the distances between the set of normal nuclei and abnormal ones
s_compute_all_distances

% Analyze the feature and sort the target gene by degree of disruption 
% of the nucleolus
s_analyze_features_vectors

