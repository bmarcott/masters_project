% Test script to sanity check min_E_def.m
c = [1; 2; 3; 4; 5; 6;];
c_home = [0; 0.1; 3; 3.1; 7; 6.9];
[A, t] = min_E_def(c, c_home);

n = length(c) / 2;
c1 = reshape(c, [2, n]);
c2 = reshape(c_home, [2, n]);
c2_est = A*c1 + repmat(t, [1, n]);

c2_est
c2

%% Compute estimate errors
err = 0.0;
for i=1:size(c1,2)
    err = err + norm(c2_est(:,i) - c2(:,i))^2;
end
fprintf('Err: %.8f\n', err);

err = sum(sum(((A*c1+repmat(t,[1,n])) - c2).^2, 1));
fprintf('    Err: %.8f    (should be equiv)\n', err);
