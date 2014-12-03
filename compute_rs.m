function [rs, norm_terms] = compute_rs(I, bs, var_b, pi_n, B)
%INPUT
%  matrix I: [h x w]
%  matrix bs:
%    Stores the bead locations: bs(b,:) -> [x,y]. image-frame.
%  float var_b:
%    Variance of the beads.
%  float pi_n:
%    Mixing coef. between uniform and gaussian mixture.
%OUTPUT
%  matrix rs:
%    Stores the responsibilities.
%      rs(b,i,j) -> resp. for bead b for pixel at (i,j)
%  matrix norm_terms:
%    Stores the norm_terms for each pixel:
%      norm_terms(i,j) -> float norm_term
% See pg. 27 (Eq. 3.5)
rs = zeros(B, size(I,1), size(I, 2));
norm_terms = zeros(size(I,1), size(I,2));
a = size(I,1)*size(I,2); % Area of image
for i=1:size(I,1)
    for j=1:size(I,2)
        if (~I(i,j))
            continue
        end
        %% Compute normalization term
        norm_term = (pi_n*B)/((1-pi_n)*a);
        for b=1:B
            mu_b = bs(:,b); % col vec
            norm_term = norm_term + mvnpdf([j;i], mu_b,var_b*eye(2));
        end
        norm_terms(i,j) = norm_term;            
        %% Compute rs_b(i,j)
        for b=1:B
            mu_b = bs(:, b); % col vec            
            rs(b,i,j) = mvnpdf([j;i], mu_b, var_b*eye(2)) / norm_term;
        end
    end
end
end
