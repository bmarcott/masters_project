% Test script to sanity check min_E_def.m
rng(42);
c = rand(10, 1);
c_home = rand(10, 1);
n = length(c) / 2;
c1 = reshape(c, [2, n]);
c2 = reshape(c_home, [2, n]);

[A, t] = min_E_def(c, c_home); % (A,t) map obj frame -> img frame
c1_est = A*c2 + repmat(t, [1, n]);

%% Compute estimate errors
err = 0.0;
for i=1:size(c1,2)
    err = err + norm(c1_est(:,i) - c1(:,i))^2;
end
fprintf('Err: %.8f\n', err);

err = compute_E_def(c1, c2, A, t);
fprintf('    Err: %.8f    (should be equiv)\n', err);
cond(A)
