Ps = zeros(S, 1);
lam_sq = diag(lam).^2;
inverse_gains = noise ./ lam_sq;
finite_inverse_gains = inverse_gains(isfinite(inverse_gains));

% 处理所有逆增益为无穷大的情况（信道增益全为零）
if isempty(finite_inverse_gains)
    disp('所有信道增益为零，无法分配功率。');
    Ps = zeros(size(inverse_gains));
    return;
end

tau_low = min(finite_inverse_gains);
tau_high = max(finite_inverse_gains) + P;
tolerance = 1e-6;
max_iterations = 100;
current_iteration = 0;
optimal_tau = 0;

while (tau_high - tau_low) > tolerance && current_iteration < max_iterations
    tau_mid = (tau_low + tau_high) / 2;
    current_allocated_ps_temp = max(0, tau_mid - inverse_gains);
    current_P_sum = sum(current_allocated_ps_temp);
    
    if abs(current_P_sum - P) < tolerance
        optimal_tau = tau_mid;
        break;
    elseif current_P_sum < P 
        tau_low = tau_mid;
    else
        tau_high = tau_mid;
    end
    
    current_iteration = current_iteration + 1;
end

% 确保即使未收敛也设置optimal_tau
if optimal_tau == 0
    optimal_tau = (tau_low + tau_high) / 2;
end

Ps = max(0, optimal_tau - inverse_gains);

% 确保总功率不超过P（处理浮点误差）
Ps = Ps * min(1, P / sum(Ps));

disp('--- 注水法功率分配结果 ---');
disp(['噪声功率 (sigma^2): ', num2str(noise)]);
disp('信道奇异值 (lambda_s):');
disp(diag(lam));
disp(['总可用功率 (Pt): ', num2str(P)]);
disp(['计算得到的注水水位 (tau): ', num2str(optimal_tau)]);
disp('分配功率 (Ps):');
disp(Ps);
disp(['总分配功率之和: ', num2str(sum(Ps))]);