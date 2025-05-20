function [p, mu, k] = waterfilling(G, P, sigma2)
% 输入参数：
%   G: 128x8 信道矩阵
%   P: 总发射功率
%   sigma2: 接收端噪声方差
% 输出参数：
%   p: 分配的功率向量（8x1），对应8个子信道
%   mu: 注水水位
%   k: 激活的子信道数量

% 1. 奇异值分解 (SVD)
[U, S, V] = svd(G, 'econ');    % 经济型SVD分解
lambda = diag(S);               % 提取奇异值（8x1向量）

% 2. 计算等效噪声功率
n = sigma2 ./ (lambda.^2);      % 等效噪声功率（8x1向量）

% 3. 对等效噪声功率进行升序排序
[n_sorted, sort_idx] = sort(n, 'ascend');

% 4. 确定注水水位μ和激活的子信道数k
k = length(n_sorted);           % 初始尝试激活全部子信道
found = false;

while ~found && k >= 1
    total_noise = sum(n_sorted(1:k));
    mu_trial = (P + total_noise) / k;
    
    if mu_trial > n_sorted(k)
        found = true;           % 找到满足条件的k
    else
        k = k - 1;              % 减少激活的子信道数
    end
end

if k < 1
    error('无法找到可行解：总功率不足或信道条件太差');
end

% 5. 计算最终水位μ
mu = (P + sum(n_sorted(1:k))) / k;

% 6. 功率分配
p_sorted = zeros(size(n_sorted));
p_sorted(1:k) = mu - n_sorted(1:k);
p_sorted = max(p_sorted, 0);    % 确保功率非负

% 恢复原始子信道顺序
p = zeros(size(p_sorted));
p(sort_idx) = p_sorted;         % 根据排序索引恢复原始顺序

% 验证总功率约束（理论上精确满足）
assert(abs(sum(p) - P) < 1e-10, '功率分配误差过大');
end