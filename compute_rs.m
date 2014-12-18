function [rs, norm_terms] = compute_rs(I, inked, bs, var_b, pi_n, B)
%INPUT
%  matrix I: [h x w]
%  matrix inked: [N x 2]
%    Stores the row and column of all inked pixels in I
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
% see pg. 23 for f_i(y), f(y) terms
rs = zeros(B, size(I,1), size(I, 2));
norm_terms = zeros(size(I,1), size(I,2));
%a = size(I,1)*size(I,2); % Area of image
a = size(inked,1);
inked0 = inked; % rows are (x,y) coords
inked0(:,1) = inked(:,2);
inked0(:,2) = inked(:,1);
inds = sub2ind(size(norm_terms), inked(:,1), inked(:,2));

if 1 % This is buggy for some reason! :(
%% Compute normalization terms (denominator)
norm_terms(inds) = (pi_n*B)/((1-pi_n)*a);
%norm_terms(inds) = (pi_n/a);
for b=1:B
    mu_b = bs(:,b); % col vec
    sig_b = var_b*eye(2);
    f_is = mvnpdf(inked0, mu_b', sig_b);
    %norm_terms(inds) = norm_terms(inds) + ((1-pi_n)/B)*f_is;
    norm_terms(inds) = norm_terms(inds) + f_is;
end
if 0
    % Visualize norm_terms
    figure; imshow(imresize(norm_terms, 5.0), []);
    fprintf('Press any key to continue.\n'); pause; close all;
end

%% Compute rs_b(i,j), ie gaussian responsibilities
for b=1:B
    mu_b = bs(:,b); % col vec
    sig_b = var_b*eye(2);
    terms = mvnpdf(inked0, mu_b', sig_b) ./ norm_terms(inds);
    if max(terms(:)) > 1.0
        mm = mvnpdf(inked0, mu_b', sig_b);
        nn = norm_terms(inds);
        ii = find(mm > nn);
        keyboard
    end
    %rs(b, inked(:,1), inked(:,2)) = reshape(terms, [1,size(inked,1), size(inked,2)]);
    ii = sub2ind(size(rs), ones([length(terms), 1])*b, inked(:,1), inked(:,2));
    rs(ii) = terms;
end

if 0
    % Visualize resps
    hfig = figure;
    II(:,:,1) = I;
    II(:,:,2) = I;
    II(:,:,3) = I;
    % Draw ctrl pts on image
    for b=1:B
        mu_b = bs(:,b);
        xc = round(mu_b(1)); yc = round(mu_b(2));
        II(yc-2:yc+2, xc, 1) = 255;
        II(yc, xc-2:xc+2, 1) = 255;
    end
    for b=1:B
        IIb = II;
        mu_b = bs(:,b);
        xc = round(mu_b(1)); yc = round(mu_b(2));
        % Highlight this ctrl pt
        IIb(yc-2:yc+2, xc, 2) = 255;
        IIb(yc, xc-2:xc+2, 2) = 255;
        % Overlay image and resp
        Ib = imresize(squeeze(rs(b,:,:)), 4.0);
        If = imfuse(imresize(IIb, 4.0), Ib, 'blend');
        imshow(If,[]);
        title(sprintf('Bead %d/%d', b, B));
        pause;
    end
    close(hfig);
end

if (max(rs(:)) > 1.0)
    keyboard
end

return;
end
if 1
    %% Old way of computing
    for ndx = 1 : size(inked,1)
        i = inked(ndx,1);
        j = inked(ndx,2);
            %% Compute normalization term
            norm_term = (pi_n*B)/((1-pi_n)*a);
            for b=1:B
                mu_b = bs(:,b); % col vec
                norm_term = norm_term + mvnpdf([j;i], mu_b,var_b*eye(2)); % 51.7%
            end
            norm_terms(i,j) = norm_term;            
            %% Compute rs_b(i,j)
            for b=1:B
                mu_b = bs(:, b); % col vec            
                rs(b,i,j) = mvnpdf([j;i], mu_b, var_b*eye(2)) / norm_term; % 46.2%
            end
    end
end

if (max(rs(:)) > 1.0)
    keyboard
end

end
