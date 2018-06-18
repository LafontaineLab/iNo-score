% Function that determine the type of the nucleus 
% (simple, multiple or 'strange')

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

function type_nucleus = f_get_type_nucleus(I_nucleus, I_nucleolus, I_BW)

% Inputs:
% *******
% I_nucleus   : image of one nucleus
% I_nucleolus : image of the nucleolus associated to the nucleus
% I_BW        : segmentation of the nucleus
%
% Output:
% *******
% type_nucleus : 1 => simple nucleus
%                2 => double/multiple (its boundary includes at least 
%                                      one strong concavity)
%                3 => strange (the nucleolus or nucleus signal is completely 
%                              too strong or too low)

  % normalization of the images
  I_nucleus_n = I_nucleus/max(max(I_nucleus));
  I_nucleo_n = I_nucleolus/max(max(I_nucleolus));

  % size of the images
  [ht, lg] = size(I_nucleus);

  % size of the nucleus
  nb_pixels = sum(sum(I_BW));

  type_nucleus = 1;

  % Test if the nucleus is considered as 'strange' or not
  % i.e. : the histogram of the intensity of the nucleus is not white
  %     or the histogram of the intensities of the nucleolus is not dark
  if (  sum(I_nucleus_n(I_BW==1)>=0.47)<nb_pixels/2 ...
      | sum(I_nucleo_n(I_BW==1)<=0.53)<nb_pixels/2)

    type_nucleus = 3;

  else

    % Test if the nucleus is a simple nucleus or a double/multiple nucleus
    % A nucleus is double/multiple if it exists at least one significant
    % concavity on its boundary
    
    perimetre = bwperim(I_BW);

    [y x] = find(perimetre==1);

    [r, s, t] = pca([x y]);

    origine = [mean(x) mean(y)]';

    % if the barycenter of the boundary point is not inside the nucleus,
    % then the nucleus is double/multiple
    if (sum(sum(I_BW(floor(origine(2)):floor(origine(2))+1, ...
  		     floor(origine(1)):floor(origine(1)+1))))~=4)

      type_nucleus = 2;

    else

      % line parallel of the new X axis
      % go through all the lines until finding a concavity
      % How many changes 'inside nucleus' <-> 'outside nucleus' ?
      trouve = 0;
      j = floor(min(s(:,2)));

      while ((j<=floor(max(s(:,2)))) & (trouve==0)) % parcours on Y

        % Points on the line of equation Y = cste
        ind = find(s(:,2)>j & s(:,2)<=j+1);
        [ssort indsort] = sort(s(ind,1));
        ssort = s(ind(indsort),:);

        if (size(ind,1)>2)
          % Distance between two consecutive points on the line Y = cste
          dd = sqrt((ssort(1:end-1,1)-ssort(2:end,1)).^2 ...
                  + (ssort(1:end-1,2)-ssort(2:end,2)).^2 );
          
          % We consider that two points are normally close,
          % if their distance is smaller than 5 pixels
          ii = find(dd>6);

          % We consider that there is a turnaround point,
          % if they are at least 3 groupes of points distant of more than 5 pixels
          if (size(ii,1)>=3)
            trouve = 1;
          end
        end

        j = j+1;

      end % while j

      if (trouve==1)
 
        type_nucleus = 2;

      else 

        % line parallel of the new Y axis
        % go through all the lines until finding a turnaround point
        % How many changes 'inside nucleus' <-> 'outside nucleus' ?
        j=floor(min(s(:,1)));

        while ((j<=floor(max(s(:,1)))) & (trouve==0)) % parcours on X

          % Points on the line of equation X = cste
          ind = find(s(:,1)>j & s(:,1)<=j+1);
          [ssort indsort] = sort(s(ind,2));
          ssort = s(ind(indsort),:);

          if size(ind,1)>2
            % Distance between two consecutive points on the line X = cste
            dd = sqrt((ssort(1:end-1,1)-ssort(2:end,1)).^2 ...
                    + (ssort(1:end-1,2)-ssort(2:end,2)).^2 );
           
            % We consider that two points are normally close,
            % if their distance is smaller than 5 pixels
            ii = find(dd>6);


            % We consider that there is a turnaround point,
            % if they are at least 3 groupes of points distant of more than 5 pixels
            if (size(ii,1)>=3)
              trouve = 1;
            end
          end

          j = j + 1;

        end % while j

        if trouve==1
          type_nucleus = 2;
        end

      end % if trouve==1 after analysis on X

    end % if center outside nucleus

  end % if 'strange'

end % function
