```
copyright 2017 All rights reserved

Pascaline Parisot (pascaline.parisot.ucl@gmail.com) 
Christophe De Vleeschouwer (christophe.devleeschouwer@uclouvain.be)
ISPGroup, Universite catholique de Louvain (Belgium)
http://sites.uclouvain.be/ispgroup/

Denis L.J. Lafontaine (denis.lafontaine@ulb.ac.be)
RNA Molecular Biology, Universite Libre de Bruxelles (Belgium)
http://www.LafontainLab.com
http://www.RibosomalProteins.com
http://www.RibosomeSynthesis.com
```

This archive contains codes to sort the target gene by index of nucleolar disruption (iNo).

- data should contain the dataset of images named as X_Y_sZ_w2.TIF (for GFP image) and X_Y_sZ_w2.TIF for (DAPI image) (where X is the name of the plate, Y the name of the well, Z the index of the images/sites in a well):
  - one directory by plate (the name of the directory is the name of the plate)
  - a "dataconfig.m" file (look at the example file to build your own file)
- results will contain all the result files
- src contains the source codes
  - Launch the "s_main" script in order to segment the nuclei, compute the features and analyze them.
