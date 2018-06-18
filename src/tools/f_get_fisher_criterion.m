% Function that computes the Fisher criterion between two sets 
% of unidimensional data (p et q)
% C = (m_p - m_q) / sqrt(s_p + s_q)
% where: m_p and m_q are the means of p and q
%        s_p and s_q are the variances of p et q (non normalized)

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

function [crit, crit2] = f_get_fisher_criterion(p, q)

% Inputs:
% p : one distribution (1xn vector)
% q : an other distribution (1xm vector) 
%
% Ouputs:
% crit  : Fisher criterion
% crit2 : square of Fisher criterion

  m_p = mean(p);
  m_q = mean(q);

  s_p = (p-m_p)*(p-m_p)';
  s_q = (q-m_q)*(q-m_q)';

  crit = (m_p - m_q)/sqrt(s_p + s_q);
  crit2 = (m_p - m_q)^2/(s_p + s_q);

end
