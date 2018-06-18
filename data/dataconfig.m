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

% File that gives the configuration of the data.

% The names of the image file are the following
%  - X_Y_sZ_gfp.TIF (for the GFP image)
%  - X_Y_sZ_dapi.TIF (for the DAPI image)
% where:
%  - X is the name of the plate
%  - Y is the name of the well (ex: A01, C05, H12)
%  - Z is the index of the site/image taken in the well

% Number of features to compute (either 5 or all)
compute_all_features = 0; % 0 => compute the 5 features
                          % 1 => compute all the features

% Number of sites/images per well
nb_site_per_well = 16;

% List of all the gene names
gene_name = {...
'eL24'
'eS8'
'uL5'
'uL18'
'uS8'
'SCR'
};

nb_gene_name = size(gene_name,1);

% List of the gene names that will be considered as normal (reference gene)
normal_gene_name = {...
'SCR'
};

% The structure contains:
%  - the name of the plate
%  - the name of the well (A05 means row A, column 5 of the plate)
%  - the name of the gene
%  - the index of the gene initialized to 0, it will be automatically filled
DATA = [...
   struct('plate', 'plate1', 'well', 'A05', 'gene_name', 'eL24', 'gene_index', 0);
   struct('plate', 'plate1', 'well', 'A06', 'gene_name', 'eL24', 'gene_index', 0);
   struct('plate', 'plate1', 'well', 'A07', 'gene_name', 'eL24', 'gene_index', 0);
   struct('plate', 'plate1', 'well', 'B12', 'gene_name', 'SCR',  'gene_index', 0);
   struct('plate', 'plate1', 'well', 'C03', 'gene_name', 'eS8',  'gene_index', 0);
   struct('plate', 'plate1', 'well', 'C04', 'gene_name', 'eS8',  'gene_index', 0);
   struct('plate', 'plate1', 'well', 'C05', 'gene_name', 'eS8',  'gene_index', 0);
   struct('plate', 'plate1', 'well', 'C06', 'gene_name', 'uL18', 'gene_index', 0);
   struct('plate', 'plate1', 'well', 'C07', 'gene_name', 'uL18', 'gene_index', 0);
   struct('plate', 'plate1', 'well', 'C08', 'gene_name', 'uL18', 'gene_index', 0);
   struct('plate', 'plate1', 'well', 'E05', 'gene_name', 'uS8',  'gene_index', 0);
   struct('plate', 'plate1', 'well', 'E06', 'gene_name', 'uS8',  'gene_index', 0);
   struct('plate', 'plate1', 'well', 'E07', 'gene_name', 'uS8',  'gene_index', 0);
   struct('plate', 'plate1', 'well', 'E08', 'gene_name', 'uL5',  'gene_index', 0);
   struct('plate', 'plate1', 'well', 'E09', 'gene_name', 'uL5',  'gene_index', 0);
   struct('plate', 'plate1', 'well', 'E10', 'gene_name', 'uL5',  'gene_index', 0);
];

% Total number of wells
nb_wells = size(DATA,1);

% gene_index settings
for i=1:nb_wells
  not_find=1;
  j = 1;
  while (not_find==1)
%[ i j strcmp(gene_name{j}, DATA(i).gene_name)]
    if (j>nb_gene_name)
      error(['gene name "' DATA(i).gene_name '" not in the list of gene name'])
    end
    if (strcmp(gene_name{j}, DATA(i).gene_name)==1)
      DATA(i).gene_index = j;
      not_find = 0;
    end
    j = j+1;
  end
end

% Set the index of the normal_gene_name
index_normal_wells = [];
gene_index = cat(1,DATA.gene_index);
for i=1:size(normal_gene_name,1)
  not_find=1;
  j = 1;
  while (not_find==1)
    if (j>nb_gene_name)
      error(['gene name "' normal_gene_name{i} '" not in the list of gene name'])
    end
    if (strcmp(gene_name{j}, normal_gene_name{i})==1)
      index_normal_wells = [index_normal_wells; find(gene_index==j)];
      not_find = 0;
    end
    j = j+1;
  end
end

% For the legends in Matlab
the_name_legend = cell(nb_wells,1);
for ind = 1:nb_wells
  the_name_legend{ind} = [DATA(ind).plate ' ' DATA(ind).well ' ' DATA(ind).gene_name ];
end
