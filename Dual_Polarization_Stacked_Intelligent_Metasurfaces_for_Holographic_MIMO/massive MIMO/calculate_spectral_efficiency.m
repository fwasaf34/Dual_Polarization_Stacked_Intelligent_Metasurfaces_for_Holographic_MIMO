function C = calculate_spectral_efficiency(G, P, noise)
    % 1. 执行注水算法分配功率
    [p, ~, ~] = waterfilling(G, P, noise);  % 调用之前的注水函数
    
    % 2. 对信道矩阵G进行SVD分解，获取奇异值
    [~, S, ~] = svd(G, 'econ');
    lambda = diag(S);  % 奇异值向量（8x1）
    
    % 3. 计算每个子信道的容量贡献
    C_i = log2(1 + (p .* lambda.^2) / noise);
    
    % 4. 总频谱效率
    C = sum(C_i);
end