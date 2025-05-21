function [p, mu, k] = waterfilling(G, P, sigma2)
[U, S, V] = svd(G, 'econ');
lambda = diag(S);
n = sigma2 ./ (lambda.^2);
[n_sorted, sort_idx] = sort(n, 'ascend');
k = length(n_sorted);
found = false;
while ~found && k >= 1
    total_noise = sum(n_sorted(1:k));
    mu_trial = (P + total_noise) / k;
    if mu_trial > n_sorted(k)
        found = true;
    else
        k = k - 1;
    end
end
if k < 1
    error('无法找到可行解：总功率不足或信道条件太差');
end
mu = (P + sum(n_sorted(1:k))) / k;
p_sorted = zeros(size(n_sorted));
p_sorted(1:k) = mu - n_sorted(1:k);
p_sorted = max(p_sorted, 0);
p = zeros(size(p_sorted));
p(sort_idx) = p_sorted;
assert(abs(sum(p) - P) < 1e-10, '功率分配误差过大');
end