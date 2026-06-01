% The corresponding paper for this code is:
% Z. He, X. Li, W. Zhang, W. Chen, C. Zhang, An efficient mixed FEM for flexoelectricity 
% based on static condensation of Lagrange multipliers, Computer Methods in Applied 
% Mechanics and Engineering, 460 (2026) 119116.
% https://doi.org/10.1016/j.cma.2026.119116
%
% This code generates the results (Fig. 7(b)and(c), classical) discussed in Section 4.4.
% 
% If you have any questions, please feel free to contact me at: hezz@zju.edu.cn

clc; clear; close all;

%% 1. 定义参数并封装到结构体中 (替代 global，提升效率)
p = struct();

% 几何参数
p.h = 5e-9;
p.a = 2 * p.h;

% 材料参数 (弹性矩阵)
p.c11 = 173E9;
p.c13 = 82E9;
p.c33 = 173E9;
p.c44 = 108E9;

% 材料尺度参数
p.l = p.h / 5*0;

% 密度
p.rho = 5996;


%% 2. 绘制 FEM 数据背景
if isfile('MFEM_StrainGradient_DispersionCurve.txt')
     femData = readmatrix('MFEM_StrainGradient_DispersionCurve.txt');
    plot(femData(:, 1), femData(:, 2:end), '.', ...
        'MarkerSize', 10, 'Color', [0.85 0.85 0.85]);
    hold on;
end

%% 3. 计算主体
data = []; 
i = 1;

for k = [0.001 0.02:0.02:1]
    tic;
    disp(["开始计算 k = " num2str(k)]);
    
    % 调用优化后的求根函数
    % 使用匿名函数将参数 p 传递给 compute_M
    da = Searoots(@(f, kv) compute_M(f, kv, p), 1E8, 1e9/2, 6e11, k);
    
    if ~isempty(da)
        L2 = plot(k * ones(size(da,1),1), da(:,1), 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r'); 
        hold on;
        drawnow;
        % 动态记录数据（自动补零适应不同数量的根）
        data(i, 1:size(da,1)+1) = [k, da(:,1)']; 
    else
        data(i, 1) = k;
    end
    i = i + 1;
    toc;
end

% 格式化图表
xlim([0, 1]);
% ylim([0, 7e10]);
if exist('L1', 'var') && exist('L2', 'var')
    legend([L1(1) L2(1)], 'FEM', '解析解', 'fontsize', 15, 'Location', "northwest");
end

title('Classical Wave Dispersion');
xlim([0, 1]);
xlabel('Wave number (\pi/(2h))');
ylabel('Frequency (Hz)');
% save('AnalyticSolution_Classic.txt',"data",'-ascii')

%% ================= 子函数区域 ================= %%

function detM = compute_M(f, k, p)
    % 高度向量化计算行列式的值，不再使用缓慢的双重 for 循环和 switch
    omega = f * 2 * pi;
    k1 = pi / p.a * k;

    % k3p 为 8x1 的列向量
    k3p = double(getK3(omega, k1, p));
    
    % 预计算重复项
    k3p2 = k3p.^2;
    k3p3 = k3p.^3;
    k3p4 = k3p.^4;

    % 批量计算 k11, k12 和 bj (全为 8x1 向量)
    k11 = -p.c11*k1^2 - p.c44.*k3p2 + p.rho*omega^2 - p.l^2*p.c11*k1^4 - p.l^2*(p.c11+p.c44)*k1^2.*k3p2 - p.l^2*p.c44.*k3p4;
    k12 = -(p.c13+p.c44)*k1.*k3p - p.l^2*(p.c13+p.c44)*k1^3.*k3p - p.l^2*(p.c13+p.c44)*k1.*k3p3;
    
    bj = -k11 ./ k12;
    % 注：原代码中 cj=0，因此所有乘以 cj 的项全部省略以提升速度
    
    % 构建 M 矩阵 (8x8)，按行直接向量化赋值
    M = zeros(4, 4);
    
    % 指数项 (1x8 行向量)
    exp_plus  = exp(1i .* k3p.' .* (+p.h/2));
    exp_minus = exp(1i .* k3p.' .* (-p.h/2));

    % Case 1 & 2 公共部分 (1x8 行向量)
    term12 = p.c44*1i.*k3p + p.c44*1i*k1.*bj + p.l^2*p.c44*(1i*k1^2.*k3p + 1i*k1^3.*bj) ...
             + p.l^2*p.c44*(1i.*k3p3 + 1i*k1.*k3p2.*bj) ...
             + p.l^2*p.c11*1i*k1^2.*k3p + p.l^2*p.c13*1i*k1.*k3p2.*bj;
    
    M(1, :) = term12.' .* exp_plus;
    M(2, :) = term12.' .* exp_minus;

    % Case 3 & 4 公共部分
    term34 = p.c13*1i*k1 + p.c33*1i.*k3p.*bj + p.l^2*p.c13*1i*k1^3 + p.l^2*p.c33*1i*k1^2.*k3p.*bj ...
             + p.l^2*p.c13*1i*k1.*k3p2 + p.l^2*p.c33*1i.*k3p3.*bj ...
             + p.l^2*p.c44*(1i*k1.*k3p2 + 1i*k1^2.*k3p.*bj);
             
    M(3, :) = term34.' .* exp_plus;
    M(4, :) = term34.' .* exp_minus;

    % 计算行列式的绝对值
    detM = abs(det(M));
end

function result = Searoots(fun, x_start, x_step, x_end, k)
    % 优化后的搜根程序：粗网格扫描 -> 局部区间细化
    
    % 1. 粗网格评估
    % 原代码相当于在 x_start 到 x_end 内每三分之一个 x_step 评估一次
    dx_grid = x_start : (x_step/3) : (x_end + x_step*4/3); 
    n_points = length(dx_grid);
    y_grid = zeros(1, n_points);
    
    for idx = 1:n_points
        y_grid(idx) = fun(dx_grid(idx), k);
    end
    
    result = [];
    
    % 2. 寻找极小值并做区间缩小(二分)迭代
    for idx = 2:(n_points - 1)
        if y_grid(idx) < y_grid(idx-1) && y_grid(idx) < y_grid(idx+1)
            % 发现局部极小值中心，开始局部细化
            nwf = dx_grid(idx);
            tx = (x_step/3) / 2;
            s1 = y_grid(idx);
            s1_initial = s1; % 记录初始极小值以便最终验证
            
            % 迭代 20 次逼近真实谷底
            for kk = 1:20
                ax = nwf - tx;
                bx = nwf + tx;
                
                y_ax = fun(ax, k);
                y_bx = fun(bx, k);
                
                if y_ax < s1 && y_ax < y_bx
                    s1 = y_ax;
                    nwf = ax;
                elseif y_bx < s1
                    s1 = y_bx;
                    nwf = bx;
                end
                tx = tx / 2;
            end
            
            % 零点判断：如果细化后的谷底值比初始值小三个量级以上，认为是根
            if (s1 * 1e3) < s1_initial
                result = [result; nwf, s1];
            end
        end
    end
end

function k3p = getK3(omega, k1, p)
    % 计算多项式系数
    Coeffs = [
        p.c33*p.c44,...
        p.c11*p.c33*k1^2 - p.c33*omega^2*p.rho - p.c44*omega^2*p.rho - p.c13^2*k1^2 - 2*p.c13*p.c44*k1^2,...
        omega^4*p.rho^2 + p.c11*p.c44*k1^4 - p.c11*k1^2*omega^2*p.rho - p.c44*k1^2*omega^2*p.rho];

    % 求解 k3 平方并开方，返回包含正负值的列向量
    % 求解 k3 平方并开方，返回包含正负值的列向量
    k3p = sqrt(roots(Coeffs));
    k3p = [k3p; -k3p];
end