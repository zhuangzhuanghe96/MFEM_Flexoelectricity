% =========================================================================
% 二维挠曲电混合有限元 (M-FEM) 求解器
% 基于位移-电势-辅助应变的多物理场强耦合计算
% =========================================================================

clc; clear;close all;

% 打印求解器初始化信息
fprintf('============================================================\n');
fprintf('   2D FlexoFEM: Flexible Element Types                      \n');
fprintf('============================================================\n');

%% 1. 全局参数与仿真配置
conf.filename      = 'mesh_Example2_7d9k.mphtxt'; 
conf.el_type       = 'Q8';       % 位移场插值单元类型 ('Q8' 或 'Q9')
conf.lam_type      = 'P0';       % 拉格朗日乘子单元类型 ('Q4' 或 'P0')
conf.condense      = true;       % 是否在单元级别执行静态凝聚
conf.analysis_type = 'static';   % 求解类型：静态分析 ('static') 

% 输出当前单元和求解策略配置
fprintf('[Config] Displacement: %s\n', conf.el_type);
fprintf('[Config] Multiplier:   %s\n', conf.lam_type);
fprintf('[Config] Condensation: %s\n', string(conf.condense));

% 设定材料物理参数
mat.E       = 139e9;     % 杨氏模量 (Pa)
mat.nu      = 0.3;       % 泊松比
mat.kappa   = 1e-9;      % 介电常数 (F/m)
mat.l_scale = 2e-6;      % 内部特征长度尺度 (m)

% 挠曲电系数 (C/m)
mat.mu_L = 3e-6;         % 纵向挠曲电系数
mat.mu_T = 1e-6;         % 横向挠曲电系数
mat.mu_S = 1e-6;         % 剪切挠曲电系数

mat.rho = 2300;          % 材料密度 (kg/m^3)

%% 2. 模型前处理与自由度映射
% 读取并解析网格文件
mesh = import_mesh_adapter(conf.filename);
if mesh.n_elems == 0
    error('未读取到网格，请检查文件!'); 
end

% 计算材料属性矩阵与多物理场数值缩放因子
mat = compute_material_matrices(mat);
scl = compute_scaling_factors(mat);

fprintf('[Scaling] Factors: S_phi = %.2e, S_psi = %.2e\n', scl.S_phi, scl.S_psi);

% 构建全局自由度 (DOF) 映射表
[dof, mesh, info] = build_dof_map_flexible(mesh, conf);

%% 3. 全局刚度矩阵组装
fprintf('[Main] 组装刚度矩阵...\n');
K_global = assemble_system_sequential(mesh, dof, mat, conf, info, scl); 

%% 4. 边界条件与外部载荷施加
% 定义并施加 Dirichlet 边界条件 (位移与电势)
fprintf('[Main] 应用 Dirichlet 边界条件...\n');
BC_List = {}; 
BC_List{end+1} = struct('tag', [4, 6], 'type', 'u', 'value', {{@(x,y) x./(x.^2+y.^2).^0.5*0.045E-6, @(x,y) y./(x.^2+y.^2).^0.5*0.045E-6}});
BC_List{end+1} = struct('tag', [5, 7], 'type', 'u', 'value', {{@(x,y) x./(x.^2+y.^2).^0.5*0.05E-6, @(x,y) y./(x.^2+y.^2).^0.5*0.05E-6}});
BC_List{end+1} = struct('tag', 3, 'type', 'u', 'value', {{NaN, 0}});
BC_List{end+1} = struct('tag', 1, 'type', 'u', 'value', {{0, NaN}});
BC_List{end+1} = struct('tag', [4, 6], 'type', 'phi', 'value', 0);
BC_List{end+1} = struct('tag', [5, 7], 'type', 'phi', 'value', 1);
BC_List{end+1} = struct('tag', [1, 3], 'type', 'psi', 'value', {{NaN, NaN, 0}}); 

[fixed_dofs, U_fixed_vals] = apply_dirichlet_bc(mesh, dof, BC_List, scl);

% 定义并施加 Neumann 边界条件 (面力或电荷)
fprintf('[Main] 应用 Neumann 载荷...\n');
Load_List = {};
F_global = apply_neumann_load(mesh, dof, Load_List);

%% 5. 线性方程组求解
free_dofs = setdiff(1:dof.n_total, fixed_dofs);
fprintf('[Main] 求解线性方程组 (DOF: %d)...\n', length(free_dofs));

% 求解系统主方程
U_sol = solve_system_dirichlet(K_global, F_global, fixed_dofs, U_fixed_vals, dof.n_total);

% 提取并逆缩放物理节点解
results = extract_results_flexible(U_sol, dof, mesh, scl);
fprintf('[Main] 求解完成. 最大位移: %.4e m\n', max(abs(results.v)));
fprintf('[Main] 求解完成. 最大电势: %.4e V\n', max(abs(results.phi)));

%% 6. 结果后处理与可视化
visualize_results(mesh, results);

fprintf('\n[Post] 绘制边界曲线...\n');
plot_boundary_result(mesh, results, 3, 'U_mag', 'arc', "k-", "markerSize", 15, 'LineWidth', 1.5);
plot_boundary_result(mesh, results, 3, 'phi', 'arc', "r-", "markerSize", 15, 'LineWidth', 1.5);

% 提取特定边界 (Tag 3) 的数据以进行验证
fprintf('\n[Data] 提取上表面电势数据...\n');
boundary_data = extract_boundary_data(mesh, results, 3, 'phi');

% 加载解析解系数进行理论对比
load("Exact_coefficient.mat");
u_exact = @(x,y) coeff(1)*x + coeff(2)/x + coeff(3)*besseli(0, x/coeff(7)) + coeff(4)*besselk(0, x/coeff(7)); 
du_dx   = @(x,y) coeff(1) - coeff(2)/x.^2 + (coeff(3)/coeff(7)) * besseli(1, x/coeff(7)) - (coeff(4)/coeff(7)) * besselk(1, x/coeff(7));

f_sum = 3E-6;
k = 1e-9;
phi_exact = @(x,y) coeff(5)*log(x) + coeff(6) + (f_sum / k) * (du_dx(x,y) + u_exact(x,y)./x);

% 计算数值解与解析解的 L2 相对误差
err  = compute_boundary_error(mesh, results, 3, 'u', u_exact);
err2 = compute_boundary_error(mesh, results, 3, 'phi', phi_exact);

% 聚合边界数据并保存验证文件
U = extract_boundary_data(mesh, results, 3, "u");
data = [U.x, U.values];
U = extract_boundary_data(mesh, results, 3, "phi");
data = [data, U.values];
U = extract_boundary_data(mesh, results, 3, "psi11");
data = [data, U.values];
U = extract_boundary_data(mesh, results, 3, "psi22");
data = [data, U.values];
% save("MFEM_staicvalidation.txt", "data", '-ascii');




%以下为模型所需调用的子函数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 前处理========================================================================

function mesh = import_mesh_adapter(filename)
    % 适配器函数：解析并构建有限元网格拓扑结构与边界信息
    [node_coords, elem_struct, boundaries] = import_mesh_mixed_v3(filename);
    plot_mesh_check_v3(node_coords, elem_struct, boundaries);
    
    if isempty(elem_struct.quad2)
        error('网格文件中未找到 quad2 (Q9) 单元！');
    end
    
    mesh.elems = elem_struct.quad2(:, 1:9);      
    mesh.quad_tags = elem_struct.quad2(:, 10);   
    
    if ~isempty(boundaries)
        mesh.edge_elems = boundaries(:, 1:end-1); 
        mesh.edge_tags  = boundaries(:, end);
    else
        mesh.edge_elems = [];
        mesh.edge_tags = [];
    end
    
    mesh.coords = node_coords;
    mesh.n_nodes = size(mesh.coords, 1);
    mesh.n_elems = size(mesh.elems, 1);
    mesh.x_max = max(mesh.coords(:,1));
    mesh.y_max = max(mesh.coords(:,2));
    mesh.y_min = min(mesh.coords(:,2));
    mesh.H = mesh.y_max - mesh.y_min;
    
    fprintf('[Mesh Adapter] 已适配: %d 节点, %d 单元 (包含边界信息)\n', mesh.n_nodes, mesh.n_elems);
end

function mat = compute_material_matrices(mat)
    % 依据材料基本参数构建弹性、介电和挠曲电耦合矩阵
    if isfield(mat, 'C') && ~isempty(mat.C)
        fprintf('[Material] 使用用户输入的弹性矩阵 mat.C (3x3).\n');
        if ~all(size(mat.C) == [3, 3])
            error('自定义 mat.C 必须是 3x3 矩阵 (对应 Voigt 顺序: 11, 22, 12)。');
        end
    else
        if ~isfield(mat, 'E') || ~isfield(mat, 'nu')
            error('未找到 mat.C，且缺少 mat.E 或 mat.nu，无法计算弹性矩阵。');
        end
        const = mat.E / ((1+mat.nu)*(1-2*mat.nu));
        mat.C = const * [1-mat.nu, mat.nu, 0; 
                         mat.nu, 1-mat.nu, 0; 
                         0,       0,      (1-2*mat.nu)/2];
    end

    if isscalar(mat.kappa)
        mat.kappa = mat.kappa * eye(2); 
    else
        if ~all(size(mat.kappa) == [2, 2])
            error('Kappa must be scalar or 2x2 matrix.');
        end
    end
    
    if isfield(mat, 'mu_mat') && ~isempty(mat.mu_mat)
        if ~all(size(mat.mu_mat) == [6, 2])
            error('自定义 mat.mu 必须是 6x2 矩阵。\n');
        end
    else
        mat.mu_mat(1, 1) = mat.mu_L; 
        mat.mu_mat(5, 2) = mat.mu_L;      
        mat.mu_mat(2, 1) = mat.mu_T; 
        mat.mu_mat(4, 2) = mat.mu_T;        
        mat.mu_mat(6, 1) = mat.mu_S; 
        mat.mu_mat(3, 2) = mat.mu_S; 
    end

    mat.D_sge = mat.l_scale^2 * blkdiag(mat.C, mat.C);
end

function scl = compute_scaling_factors(mat)
    % 评估不同物理场的量级并生成预处理缩放因子以改善矩阵条件数
    C_ref = max(diag(mat.C));
    if isscalar(mat.kappa)
        k_ref = mat.kappa;
    else
        k_ref = max(abs(diag(mat.kappa)));
    end
    
    scl.S_phi = sqrt(C_ref / k_ref);
    scl.S_psi = 1.0 / mat.l_scale; 
    scl.S_lam = C_ref;
    
    fprintf('[Scaling] S_phi = %.2e (Voltage Scaling)\n', scl.S_phi);
    fprintf('[Scaling] S_psi = %.2e (Gradient Scaling)\n', scl.S_psi);
end

function [dof, mesh, info] = build_dof_map_flexible(mesh, conf)
    % 动态分配单元与节点的自由度映射 (支持 Q8/Q9 及气泡函数扩展)
    if strcmp(conf.el_type, 'Q9')
        info.idx_u_nodes = 1:9; 
    elseif strcmp(conf.el_type, 'Q8')
        info.idx_u_nodes = 1:8; 
    else
        error('Unsupported el_type: %s', conf.el_type);
    end
    
    info.n_bub_dof = 3; 
    
    if strcmp(conf.lam_type, 'P0')
        info.n_lam_dof = 3;  
    elseif strcmp(conf.lam_type, 'Q4')
        info.n_lam_dof = 12; 
    else
        error('Unsupported lam_type: %s', conf.lam_type);
    end
    
    mesh.is_corner = false(mesh.n_nodes, 1);
    unique_corners = unique(mesh.elems(:, 1:4)); 
    mesh.is_corner(unique_corners) = true;
    
    is_u_node = false(mesh.n_nodes, 1);
    used_nodes = mesh.elems(:, info.idx_u_nodes); 
    is_u_node(unique(used_nodes(:))) = true;
    
    dof.map = zeros(mesh.n_nodes, 6); 
    curr = 0;
    
    for i = 1:mesh.n_nodes
        if is_u_node(i)
            dof.map(i, 1:2) = curr + [1, 2];
            curr = curr + 2;
        end
        
        if mesh.is_corner(i)
            dof.map(i, 3) = curr + 1;
            dof.map(i, 4:6) = curr + (2:4);
            curr = curr + 4;
        end
    end
    
    dof.elem_map = [];
    if ~conf.condense
        n_internal = info.n_bub_dof + info.n_lam_dof;
        dof.elem_map = zeros(mesh.n_elems, n_internal);
        for e = 1:mesh.n_elems
            dof.elem_map(e, :) = curr + (1:n_internal);
            curr = curr + n_internal;
        end
        fprintf('[DOF] 全耦合模式: 已分配 Bubble/Lambda 全局自由度。\n');
    else
        fprintf('[DOF] 静态凝聚模式: Bubble/Lambda 仅在单元内处理。\n');
    end
    
    dof.n_total = curr;
    
    fprintf('[DOF] 统计:\n');
    fprintf('      Total DOFs: %d\n', dof.n_total);
    fprintf('      Displacement Nodes: %d (Type: %s)\n', sum(is_u_node), conf.el_type);
    fprintf('      Corner Nodes:       %d\n', sum(mesh.is_corner));
    if ~conf.condense
        fprintf('      Internal DOFs:      %d (Bubble+Lambda)\n', mesh.n_elems * (info.n_bub_dof + info.n_lam_dof));
    end
end
%% 刚度组装========================================================================
function K = assemble_system_sequential(mesh, dof, mat, conf, info, scl)
    % 遍历每个单元并组装全局刚度稀疏矩阵，包含静态凝聚判断逻辑
    est_nz = mesh.n_elems * 64^2; 
    I = zeros(est_nz, 1); J = zeros(est_nz, 1); V = zeros(est_nz, 1);
    nz_cnt = 0;
    
    g_pt = [-0.7745966692, 0.0, 0.7745966692]; 
    g_wt = [0.5555555556, 0.8888888889, 0.5555555556];
    
    fprintf('[Assembly] 开始组装 (Condensation: %s)...\n', string(conf.condense));

    for e = 1:mesh.n_elems
        nodes = mesh.elems(e, :);
        u_nodes = nodes(info.idx_u_nodes); 
        corners = nodes(1:4);
        
        dofs_u   = reshape(dof.map(u_nodes, 1:2)', [], 1);
        dofs_phi = dof.map(corners, 3);
        dofs_psi = reshape(dof.map(corners, 4:6)', [], 1); 
        
        retained_dofs = [dofs_u; dofs_phi; dofs_psi];
        
        if conf.condense
            ele_global_dofs = retained_dofs;
        else
            if isfield(dof, 'elem_map') && ~isempty(dof.elem_map)
                local_dofs = dof.elem_map(e, :)';
                ele_global_dofs = [retained_dofs; local_dofs];
            else
                error('配置要求全耦合求解，但未找到 dof.elem_map。请检查 build_dof_map_flexible。');
            end
        end
        
        Ke = element_stiffness_flexible(mesh.coords(nodes, :), mat, g_pt, g_wt, conf, info, scl);
        
        n_edof = length(ele_global_dofs);
        if size(Ke, 1) ~= n_edof
            error('单元 %d 刚度矩阵尺寸与全局自由度映射数不匹配。', e);
        end
        
        [Jg, Ig] = meshgrid(ele_global_dofs, ele_global_dofs);
        rng = nz_cnt + (1:n_edof^2);
        I(rng) = Ig(:); J(rng) = Jg(:); V(rng) = Ke(:);
        nz_cnt = nz_cnt + n_edof^2;
    end
    
    K = sparse(I(1:nz_cnt), J(1:nz_cnt), V(1:nz_cnt), dof.n_total, dof.n_total);
end

function Ke_final = element_stiffness_flexible(xy_all, mat, g_pt, g_wt, conf, info, scl)
    % 基于高斯积分评估局部单元刚度矩阵的所有子块(包括力学、电学及其耦合块)
    n_u_node = length(info.idx_u_nodes);
    n_ret = n_u_node*2 + 4 + 12; 
    n_total = n_ret + info.n_bub_dof + info.n_lam_dof; 
    
    Ke_aug = zeros(n_total, n_total);
    
    idx_u   = 1:n_u_node*2;
    idx_phi = n_u_node*2 + (1:4);
    idx_psi = n_u_node*2 + 4 + (1:12);
    idx_bub = n_ret + (1:info.n_bub_dof);
    idx_lam = n_ret + info.n_bub_dof + (1:info.n_lam_dof);
    
    xy_u = xy_all(info.idx_u_nodes, :);
    S_phi = scl.S_phi;
    S_psi = scl.S_psi;
    S_lam = scl.S_lam;
    
    for i = 1:3
        for j = 1:3
            xi = g_pt(i); eta = g_pt(j); w = g_wt(i) * g_wt(j);
            
            if strcmp(conf.el_type, 'Q9')
                [Nu, dNu_dxi] = shape_Q9(xi, eta); J_jac = dNu_dxi * xy_u; 
            else
                [Nu, dNu_dxi] = shape_Q8(xi, eta); J_jac = dNu_dxi * xy_u; 
            end
            detJ = det(J_jac); dNu_dx = J_jac \ dNu_dxi;
            
            [N4, dN4_dxi] = shape_Q4(xi, eta); dN4_dx = J_jac \ dN4_dxi;
            [Nb, dNb_dxi] = shape_Bubble(xi, eta); dNb_dx = J_jac \ dNb_dxi;

            Bu = build_Bu(dNu_dx, n_u_node); 
            Bphi = build_Bphi(dN4_dx)*S_phi; 
            
            Bpsi_n = build_Bpsi(N4)* S_psi;       
            B_gpsi_n = build_B_gpsi(dN4_dx)* S_psi;

            Bg = Bu; 
            
            Bpsi_b = (Nb * eye(3))* S_psi;         
            B_gpsi_b = build_B_gpsi_bubble(dNb_dx)* S_psi;
            
            if strcmp(conf.lam_type, 'P0')
                B_lam = eye(3); 
            else
                B_lam = build_Bpsi(N4)* S_psi; 
            end
            
            Ke_aug(idx_u, idx_u)     = Ke_aug(idx_u, idx_u)     + (Bu' * mat.C * Bu) * detJ * w;
            Ke_aug(idx_phi, idx_phi) = Ke_aug(idx_phi, idx_phi) - (Bphi' * mat.kappa * Bphi) * detJ * w;
            k_flx_n = (B_gpsi_n' * mat.mu_mat * Bphi) * detJ * w;
            Ke_aug(idx_psi, idx_phi) = Ke_aug(idx_psi, idx_phi) - k_flx_n; 
            Ke_aug(idx_phi, idx_psi) = Ke_aug(idx_phi, idx_psi) - k_flx_n';
            Ke_aug(idx_psi, idx_psi) = Ke_aug(idx_psi, idx_psi) + (B_gpsi_n' * mat.D_sge * B_gpsi_n) * detJ * w;
            
            k_sge_bb = (B_gpsi_b' * mat.D_sge * B_gpsi_b) * detJ * w;
            Ke_aug(idx_bub, idx_bub) = Ke_aug(idx_bub, idx_bub) + k_sge_bb;
            
            k_flx_b = (B_gpsi_b' * mat.mu_mat * Bphi) * detJ * w;
            Ke_aug(idx_bub, idx_phi) = Ke_aug(idx_bub, idx_phi) - k_flx_b; 
            Ke_aug(idx_phi, idx_bub) = Ke_aug(idx_phi, idx_bub) - k_flx_b';
            
            k_sge_bn = (B_gpsi_b' * mat.D_sge * B_gpsi_n) * detJ * w;
            Ke_aug(idx_bub, idx_psi) = Ke_aug(idx_bub, idx_psi) + k_sge_bn; 
            Ke_aug(idx_psi, idx_bub) = Ke_aug(idx_psi, idx_bub) + k_sge_bn';
            
            q_lp = (B_lam' * eye(3) * Bpsi_n) * detJ * w * S_lam;
            Ke_aug(idx_lam, idx_psi) = Ke_aug(idx_lam, idx_psi) + q_lp; 
            Ke_aug(idx_psi, idx_lam) = Ke_aug(idx_psi, idx_lam) + q_lp';
            
            q_lb = (B_lam' * eye(3) * Bpsi_b) * detJ * w * S_lam;
            Ke_aug(idx_lam, idx_bub) = Ke_aug(idx_lam, idx_bub) + q_lb; 
            Ke_aug(idx_bub, idx_lam) = Ke_aug(idx_bub, idx_lam) + q_lb';
            
            q_lu = - (B_lam' * eye(3) * Bg) * detJ * w * S_lam;
            Ke_aug(idx_lam, idx_u)   = Ke_aug(idx_lam, idx_u)   + q_lu; 
            Ke_aug(idx_u, idx_lam)   = Ke_aug(idx_u, idx_lam)   + q_lu';
        end
    end
    
    Ke_final = process_condensation_sequential(Ke_aug, conf, info, n_ret);
end

function Ke_out = process_condensation_sequential(Ke_aug, conf, info, n_ret)
    % 依据 Schur 补定理对气泡自由度和乘子进行顺序消除（静态凝聚）
    if ~conf.condense
        Ke_out = Ke_aug;
        return;
    end

    n_bub = info.n_bub_dof;
    n_lam = info.n_lam_dof;
    
    idx_r = 1:n_ret;                
    idx_b = n_ret + (1:n_bub);      
    idx_l = n_ret + n_bub + (1:n_lam);
    
    K_rr = Ke_aug(idx_r, idx_r);
    K_rb = Ke_aug(idx_r, idx_b); K_br = Ke_aug(idx_b, idx_r);
    K_rl = Ke_aug(idx_r, idx_l); K_lr = Ke_aug(idx_l, idx_r);
    
    K_bb = Ke_aug(idx_b, idx_b);
    K_bl = Ke_aug(idx_b, idx_l); K_lb = Ke_aug(idx_l, idx_b);
    
    K_ll = Ke_aug(idx_l, idx_l); 
    
    invK_bb = eye(3)/(K_bb); 
    
    K_rr_s = K_rr - K_rb * invK_bb * K_br;
    K_rl_s = K_rl - K_rb * invK_bb * K_bl;
    K_lr_s = K_lr - K_lb * invK_bb * K_br;
    K_ll_s = K_ll - K_lb * invK_bb * K_bl; 
    
    if rcond(K_ll_s) < 1e-16
        disp("乘子矩阵条件数过大，请检查原因!!")
        K_ll_s = K_ll_s + 1e-12 * mean(diag(K_rr_s)) * eye(size(K_ll_s));
    end
    
    Ke_out = K_rr_s - K_rl_s * (K_ll_s \ K_lr_s);
end

function M = assemble_mass_matrix(mesh, dof, mat, conf, info)
    % 提取并组装整体系统的位移自由度一致质量矩阵
    est_nz = mesh.n_elems * (18)^2; 
    I = zeros(est_nz, 1); J = zeros(est_nz, 1); V = zeros(est_nz, 1);
    nz_cnt = 0;
    
    g_pt = [-0.7745966692, 0.0, 0.7745966692]; 
    g_wt = [0.5555555556, 0.8888888889, 0.5555555556];
    
    rho = mat.rho;
    
    fprintf('[Assembly] 组装质量矩阵 (Density: %.1f)...\n', rho);

    for e = 1:mesh.n_elems
        nodes = mesh.elems(e, :);
        u_nodes = nodes(info.idx_u_nodes);
        n_u_node = length(u_nodes);
        dofs_u = reshape(dof.map(u_nodes, 1:2)', [], 1);
        xy_u = mesh.coords(u_nodes, :);
        Me = zeros(n_u_node*2, n_u_node*2);
        
        for i = 1:3
            for j = 1:3
                xi = g_pt(i); eta = g_pt(j); w = g_wt(i) * g_wt(j);
                
                if strcmp(conf.el_type, 'Q9')
                    [Nu, dNu_dxi] = shape_Q9(xi, eta);
                else
                    [Nu, dNu_dxi] = shape_Q8(xi, eta);
                end
                
                J_jac = dNu_dxi * xy_u;
                detJ = det(J_jac);
                
                Nu_mat = zeros(2, n_u_node*2);
                for k = 1:n_u_node
                    col = (k-1)*2 + [1,2];
                    Nu_mat(1, col(1)) = Nu(k);
                    Nu_mat(2, col(2)) = Nu(k);
                end
                
                Me = Me + rho * (Nu_mat' * Nu_mat) * detJ * w;
            end
        end
        
        n_edof = length(dofs_u);
        [Jg, Ig] = meshgrid(dofs_u, dofs_u);
        rng = nz_cnt + (1:n_edof^2);
        I(rng) = Ig(:); J(rng) = Jg(:); V(rng) = Me(:);
        nz_cnt = nz_cnt + n_edof^2;
    end
    
    M = sparse(I(1:nz_cnt), J(1:nz_cnt), V(1:nz_cnt), dof.n_total, dof.n_total);
end
%% 边界约束========================================================================
function [fixed_dofs, U_fixed_vals] = apply_dirichlet_bc(mesh, dof, BC_List, scl)
    % 评估边界几何并应用包含物理量逆向缩放逻辑的狄利克雷边界约束
    n_total = dof.n_total;
    is_fixed = false(n_total, 1);
    U_fixed_vals = zeros(n_total, 1);
    
    for i = 1:length(BC_List)
        bc = BC_List{i};
        nodes = get_nodes_by_tag(mesh, bc.tag);
        
        if isempty(nodes)
            fprintf('     [Warn] Tag %d 未找到节点，跳过。\n', bc.tag); continue; 
        end
        
        node_coords = mesh.coords(nodes, :); 
        
        switch lower(bc.type)
            case 'u',   target_nodes = nodes; cols = [1, 2];
            case 'phi', target_nodes = nodes(mesh.is_corner(nodes)); cols = 3;
            case 'psi', target_nodes = nodes(mesh.is_corner(nodes)); cols = [4, 5, 6];
            otherwise,  error('未知 BC 类型');
        end
        
        vals_input = bc.value;
        if ~iscell(vals_input)
            vals_input = num2cell(vals_input);
        end
        
        for k = 1:length(cols)
            if k > length(vals_input), continue; end
            val_item = vals_input{k};
            
            if isnumeric(val_item) && isnan(val_item), continue; end
            
            if length(target_nodes) ~= length(nodes)
                current_coords = mesh.coords(target_nodes, :);
            else
                current_coords = node_coords;
            end
            
            dof_indices = dof.map(target_nodes, cols(k));
            valid_mask = dof_indices > 0;
            dof_indices = dof_indices(valid_mask);
            current_coords = current_coords(valid_mask, :);
            
            is_fixed(dof_indices) = true;

             scale_val = 1.0;
             if cols(k) == 3        
                 scale_val = 1.0 / scl.S_phi;
             elseif cols(k) >= 4    
                 scale_val = 1.0 / scl.S_psi;
             end
            
            if isa(val_item, 'function_handle')
                calc_vals = val_item(current_coords(:,1), current_coords(:,2));
                U_fixed_vals(dof_indices) = calc_vals * scale_val;
            else
                U_fixed_vals(dof_indices) = val_item * scale_val;
            end
        end
    end
    fixed_dofs = find(is_fixed);
end

function F = apply_neumann_load(mesh, dof, Load_List)
    % 评估边界标签分布并执行边缘积分以构建系统载荷向量
    F = sparse(dof.n_total, 1);
    for i = 1:length(Load_List)
        load = Load_List{i};
        mask = ismember(mesh.edge_tags, load.tag);
        edge_indices = find(mask);
        
        if isempty(edge_indices)
            if ~isempty(load.tag)
                fprintf('     [Warn] Tags %s 未找到关联的边界单元，跳过。\n', mat2str(load.tag));
            end
            continue; 
        end
        
        edges = mesh.edge_elems(edge_indices, :);
        fprintf('     [Neumann] Tags %s (%s): 在 %d 条边上积分...\n', mat2str(load.tag), load.type, length(edge_indices));
        
        if ~iscell(load.value)
            load.value = num2cell(load.value); 
        end
        F = integrate_boundary_load(F, edges, mesh.coords, dof, load);
    end
end

function F = integrate_boundary_load(F, edges, coords, dof, load)
    % 具体实施沿曲线的三点高斯边界积分计算面力与面电荷
    g_pt = [-0.7745966692, 0.0, 0.7745966692];
    g_wt = [0.5555555556, 0.8888888889, 0.5555555556];
    n_edges = size(edges, 1);
    
    for e = 1:n_edges
        nodes = edges(e, :); 
        xy_nodes = coords(nodes, :); 
        
        for k = 1:3
            xi = g_pt(k); w = g_wt(k);
            [N_quad, dN_dxi] = shape_1D_quad(xi);
            
            dxdxi = dN_dxi * xy_nodes(:,1);
            dydxi = dN_dxi * xy_nodes(:,2);
            J_line = sqrt(dxdxi^2 + dydxi^2);
            
            x_g = N_quad * xy_nodes(:,1);
            y_g = N_quad * xy_nodes(:,2);
            
            switch lower(load.type)
                case 'traction'
                    val_tx = load.value{1};
                    if isa(val_tx, 'function_handle'), tx = val_tx(x_g, y_g); else, tx = val_tx; end
                    
                    val_ty = load.value{2};
                    if isa(val_ty, 'function_handle'), ty = val_ty(x_g, y_g); else, ty = val_ty; end
                    
                    for n = 1:3
                        global_node = nodes(n);
                        dof_u = dof.map(global_node, 1);
                        dof_v = dof.map(global_node, 2);
                        F(dof_u) = F(dof_u) + N_quad(n) * tx * J_line * w;
                        F(dof_v) = F(dof_v) + N_quad(n) * ty * J_line * w;
                    end
                    
                case 'charge'
                    val_q = load.value{1};
                    if isa(val_q, 'function_handle'), q = val_q(x_g, y_g); else, q = val_q; end
                    
                    N_lin = shape_1D_linear(xi);
                    local_idx = [1, 2]; 
                    
                    for n = 1:2
                        global_node = nodes(local_idx(n));
                        dof_phi = dof.map(global_node, 3);
                        if dof_phi > 0
                            F(dof_phi) = F(dof_phi) + N_lin(n) * q * J_line * w;
                        end
                    end
            end
        end
    end
end

function [T, slave_dofs] = apply_floquet_bc(mesh, dof, periodic_pairs, k_vec)
    % 基于波矢 k 建立系统周期性边界相位关系矩阵（Floquet 理论）
    n_total = dof.n_total;
    constraints = cell(n_total, 1);
    is_slave = false(n_total, 1);
    
    fprintf('[Floquet] 构建周期性约束 (k = [%.2f, %.2f])...\n', k_vec(1), k_vec(2));
    
    for p = 1:length(periodic_pairs)
        pair = periodic_pairs(p);
        [node_pairs, dist_vec] = match_periodic_nodes(mesh, pair.master_tag, pair.slave_tag);
        
        if isempty(node_pairs)
            warning('Tags %d -> %d 未找到匹配的周期性节点对。', pair.master_tag, pair.slave_tag);
            continue;
        end
        
        phase = exp(-1i * dot(k_vec, dist_vec)); 
        fprintf('     Pair %d->%d: %d nodes matched. Phase = %.2f + %.2fi\n', pair.master_tag, pair.slave_tag, size(node_pairs, 1), real(phase), imag(phase));
            
        for i = 1:size(node_pairs, 1)
            m_node = node_pairs(i, 1);
            s_node = node_pairs(i, 2);
            for col = 1:size(dof.map, 2)
                m_dof = dof.map(m_node, col);
                s_dof = dof.map(s_node, col);
                if m_dof > 0 && s_dof > 0
                    if ~is_slave(s_dof)
                        constraints{s_dof} = struct('master', m_dof, 'val', phase);
                        is_slave(s_dof) = true;
                    end
                end
            end
        end
    end
    
    independent_dofs = find(~is_slave);
    n_reduced = length(independent_dofs);
    old_2_new = zeros(n_total, 1);
    old_2_new(independent_dofs) = 1:n_reduced;
    
    T_rows = zeros(n_total, 1);
    T_cols = zeros(n_total, 1);
    T_vals = zeros(n_total, 1);
    count = 0;
    
    for row = 1:n_total
        count = count + 1;
        T_rows(count) = row;
        
        if ~is_slave(row)
            T_cols(count) = old_2_new(row);
            T_vals(count) = 1.0;
        else
            c = constraints{row};
            master_global = c.master;
            chain_limit = 5; chain_iter = 0;
            current_phase = c.val;
            
            while is_slave(master_global) && chain_iter < chain_limit
                c_next = constraints{master_global};
                master_global = c_next.master;
                current_phase = current_phase * c_next.val;
                chain_iter = chain_iter + 1;
            end
            
            if is_slave(master_global)
                error('检测到循环约束或过深的约束链，请检查周期性定义。');
            end
            
            T_cols(count) = old_2_new(master_global);
            T_vals(count) = current_phase;
        end
    end
    
    T = sparse(T_rows, T_cols, T_vals, n_total, n_reduced);
    slave_dofs = find(is_slave);
end

function [pairs, dist_vec] = match_periodic_nodes(mesh, tag_m, tag_s)
    % 对齐并匹配周期边界的主从节点坐标空间位置
    nodes_m = get_nodes_by_tag(mesh, tag_m);
    nodes_s = get_nodes_by_tag(mesh, tag_s);
    
    if length(nodes_m) ~= length(nodes_s)
        warning('主从边界节点数量不一致 (%d vs %d)，网格可能非周期性。', length(nodes_m), length(nodes_s));
    end
    
    xy_m = mesh.coords(nodes_m, :);
    xy_s = mesh.coords(nodes_s, :);
    
    center_m = mean(xy_m, 1);
    center_s = mean(xy_s, 1);
    dist_vec = center_s - center_m;
    
    if abs(dist_vec(1)) > abs(dist_vec(2)), sort_axis = 2; else, sort_axis = 1; end
    
    [~, idx_m] = sort(xy_m(:, sort_axis));
    [~, idx_s] = sort(xy_s(:, sort_axis));
    pairs = [nodes_m(idx_m), nodes_s(idx_s)];
    
    sorted_xy_m = xy_m(idx_m, :);
    sorted_xy_s = xy_s(idx_s, :);
    diffs = (sorted_xy_s - sorted_xy_m) - dist_vec;
    err = max(sqrt(sum(diffs.^2, 2)));
    
    if err > 1e-4 * norm(dist_vec)
        warning('周期性节点匹配误差较大 (%.2e)，网格可能未对齐。', err);
    end
end

function T_final = merge_floquet_dirichlet(T_floquet, fixed_dofs)
    % 融合约束：将位移锁定限制更新至广义系统变换矩阵
    if isempty(fixed_dofs)
        T_final = T_floquet;
        return;
    end
    fprintf('[BC] 正在合并 Dirichlet 约束到 Floquet 矩阵...\n');
    rows_to_fix = fixed_dofs;
    [~, cols_impacted] = find(T_floquet(rows_to_fix, :));
    cols_to_remove = unique(cols_impacted);
    
    if ~isempty(cols_to_remove)
        fprintf('     检测到 %d 个全局自由度被固定，移除了 T 矩阵的 %d 列。\n', length(rows_to_fix), length(cols_to_remove));
    end
    
    T_final = T_floquet;
    T_final(:, cols_to_remove) = [];
    if size(T_final, 2) == 0, warning('所有自由度都被 Dirichlet 或 Floquet 约束消除了，系统为空！'); end
end
%% 系统求解===============================================================
function U = solve_system_dirichlet(K, F, fixed_dofs, U_fixed_vals, n_total)
    % 通用稳态解算核心，移除已知主自由度并利用高斯消去求解缩减系统
    free_dofs = setdiff(1:n_total, fixed_dofs);
    U = zeros(n_total, 1);
    U(fixed_dofs) = U_fixed_vals(fixed_dofs);
    F_reduced = F(free_dofs) - K(free_dofs, fixed_dofs) * U(fixed_dofs);
    U(free_dofs) = K(free_dofs, free_dofs) \ F_reduced;
end

function [V_modes, freqs] = solve_eigen_system(K, M, fixed_dofs, n_total, num_modes)
    % 消除固定自由度后提取广义结构特征值与对应的基频振型
    is_fixed = false(n_total, 1);
    is_fixed(fixed_dofs) = true;
    free_dofs = find(~is_fixed);
    
    K_free = K(free_dofs, free_dofs);
    M_free = M(free_dofs, free_dofs);
    
    sigma = 0.1; 
    opts.issym = true; 
    opts.isreal = true;
    
    fprintf('     [Eigen] Solving generalized eigenvalue problem (Size: %d)...\n', length(free_dofs));
    t_start = tic;
    
    try
        [V_free, D] = eigs(K_free, M_free, num_modes, 'sm', opts);
    catch ME
        warning('eigs 求解失败，尝试使用带 shift 的求解...');
        [V_free, D] = eigs(K_free, M_free, num_modes, sigma, opts);
    end
    fprintf('     [Eigen] Solved in %.2f s.\n', toc(t_start));
    
    eig_vals = diag(D);
    eig_vals(eig_vals < 0) = 0;
    freqs = sqrt(eig_vals) / (2*pi);
    
    [freqs, sort_idx] = sort(freqs);
    V_free = V_free(:, sort_idx);
    V_modes = zeros(n_total, num_modes);
    V_modes(free_dofs, :) = V_free;
end

function [V_full, freqs] = solve_eigen_floquet(K, M, T, num_modes, use_condensation)
    % 带 Floquet 相位转换的模态解算器，支持频散特性提取
    if nargin < 5, use_condensation = false; end

    fprintf('[Floquet] Step 1: 边界投影 (Boundary Reduction)... '); 
    t0 = tic;
    K_red = T' * K * T;
    M_red = T' * M * T;
    
    K_red = (K_red + K_red') / 2;
    M_red = (M_red + M_red') / 2;
    fprintf('耗时 %.2f s (Size: %d, Condensation: %s)\n', toc(t0), size(K_red, 1), string(use_condensation));
    
    opts.issym = false; 
    opts.isreal = false;
    opts.disp = 0;

    if use_condensation
        diag_M = real(diag(M_red)); 
        max_mass = max(diag_M);
        mass_tol = 1e-12 * max_mass;
        
        idx_d = find(diag_M > mass_tol);  
        idx_s = find(diag_M <= mass_tol); 
        n_d = length(idx_d);
        n_s = length(idx_s);
        
        fprintf('[Floquet] Step 2: 静态凝聚 (Guyan Reduction)...\n');
        if n_d < num_modes
            warning('动态自由度 (%d) 少于请求的模态数 (%d)，已自动调整。', n_d, num_modes);
            num_modes = n_d;
        end
        
        t0 = tic;
        K_dd = K_red(idx_d, idx_d); K_ds = K_red(idx_d, idx_s);
        K_sd = K_red(idx_s, idx_d); K_ss = K_red(idx_s, idx_s);
        M_cond = M_red(idx_d, idx_d); 
        K_cond = K_dd - K_ds * (K_ss \ K_sd);
        K_cond = (K_cond + K_cond') / 2; 
        
        fprintf('     凝聚完成，耗时 %.2f s. Final Size: %d\n', toc(t0), n_d);
        
        try
            [V_d, D] = eigs(K_cond, M_cond, num_modes, 'smallestabs', opts);
        catch
            warning('凝聚系统求解困难，尝试添加 shift...');
            [V_d, D] = eigs(K_cond, M_cond, num_modes, 0.1, opts);
        end
        
        lambda = diag(D);
        freqs = sqrt(abs(real(lambda))) / (2*pi);
        [freqs, sort_idx] = sort(freqs);
        V_d = V_d(:, sort_idx);
        V_s = - K_ss \ (K_sd * V_d);
        V_red = zeros(size(K_red, 1), num_modes);
        V_red(idx_d, :) = V_d;
        V_red(idx_s, :) = V_s;
    else
        fprintf('[Floquet] Step 2: 直接求解全系统 (Skip Condensation)...\n');
        sigma = 1; 
        t0 = tic;
        try
            [V_red, D] = eigs(K_red, M_red, num_modes, sigma, opts);
        catch ME
            error('全系统求解失败。请尝试开启 use_condensation=true 或调整 sigma 值。');
        end
        fprintf('     求解完成，耗时 %.2f s.\n', toc(t0));
        
        lambda = diag(D);
        abs_lambda = abs(real(lambda));
        freqs = sqrt(abs_lambda) / (2*pi);
        [freqs, sort_idx] = sort(freqs);
        V_red = V_red(:, sort_idx);
    end
    V_full = T * V_red;
end

function [history, time_vec] = solve_transient_newmark(K, M, C, F_static, fixed_dofs, fixed_vals_amp, dof, conf)
    % 在时域内显式更新系统状态量 (采用 Newmark-beta 动力学方法)
    dt = conf.time.dt;
    T_end = conf.time.T_end;
    beta = conf.time.beta;
    gamma = conf.time.gamma;
    
    n_total = dof.n_total;
    n_steps = ceil(T_end / dt);
    
    is_global_format = (length(fixed_vals_amp) == n_total);
    if is_global_format
        fprintf('[Transient] 检测到固定值采用全局向量格式 (Size: %d)\n', n_total);
    end
    
    a0 = 1 / (beta * dt^2);
    a1 = gamma / (beta * dt);
    a2 = 1 / (beta * dt);
    a3 = 1 / (2*beta) - 1;
    a4 = gamma/beta - 1;
    a5 = dt/2 * (gamma/beta - 2);
    
    u = zeros(n_total, 1);
    u_dot = zeros(n_total, 1);
    u_ddot = zeros(n_total, 1);
    
    fprintf('[Transient] 构建有效刚度矩阵 K_eff...\n');
    K_eff = K + a0 * M + a1 * C;
    K_eff_coupling = K_eff(:, fixed_dofs); 
    [K_eff_bc, ~, fixed_indices] = apply_bc_matrix_method(K_eff, fixed_dofs);
    
    fprintf('[Transient] 分解 K_eff (Pre-factorization)...\n');
    [L_eff, U_eff, P_eff, Q_eff] = lu(K_eff_bc); 
    
    time_vec = zeros(n_steps + 1, 1); 
    time_vec(1) = 0;
    history = struct('t', {}, 'u_vec', {}, 'max_v', {});
    
    save_cnt = 1;
    history(save_cnt).t = 0;
    history(save_cnt).u_vec = u; 
    
    mask_v = (dof.map(:, 2) > 0); 
    if any(mask_v)
        history(save_cnt).max_v = 0;
        history(save_cnt).min_v = 0;
    end
    
    fprintf('[Transient] 开始时间步循环 (%d steps, Load: %s)...\n', n_steps, conf.load.type);
    t_start = tic;
    
    for step = 1:n_steps
        t = step * dt;
        switch lower(conf.load.type)
            case 'sine'
                freq = conf.load.freq;
                omega = 2 * pi * freq;
                time_factor = 1.0 * sin(omega * t);
            case 'step'
                ramp_time = 10 * dt;
                if t < ramp_time, time_factor = t / ramp_time; else, time_factor = 1.0; end
            otherwise
                error('Unknown load type: %s', conf.load.type);
        end
        
        F_t = F_static * time_factor;
        if is_global_format
            vals_current_step = fixed_vals_amp(fixed_indices) * time_factor;
        else
            vals_current_step = fixed_vals_amp * time_factor;
        end
        
        term_M = a0*u + a2*u_dot + a3*u_ddot;
        term_C = a1*u + a4*u_dot + a5*u_ddot;
        F_eff = F_t + M * term_M + C * term_C;
        F_eff = F_eff - K_eff_coupling * vals_current_step;
        F_eff(fixed_indices) = vals_current_step;
        
        u_next = Q_eff * (U_eff \ (L_eff \ (P_eff * F_eff)));
        a_next = a0 * (u_next - u) - a2 * u_dot - a3 * u_ddot;
        v_next = u_dot + dt * ((1-gamma)*u_ddot + gamma*a_next);
        
        u = u_next; u_dot = v_next; u_ddot = a_next;
        time_vec(step + 1) = t; 
        
        if mod(step, conf.time.save_freq) == 0
            save_cnt = save_cnt + 1;
            history(save_cnt).t = t;
            if any(mask_v)
                v_disp = u(mask_v);
                history(save_cnt).max_v = max(v_disp);
                history(save_cnt).min_v = min(v_disp);
            end
            history(save_cnt).u_vec = u; 
        end
        if mod(step, 100) == 0
            fprintf('   Step %d/%d (t=%.2e s) | Factor=%.2f\n', step, n_steps, t, time_factor);
        end
    end
    fprintf('[Transient] 完成. 耗时 %.2f s.\n', toc(t_start));
end

function [K_mod, F_mod, fixed_idx] = apply_bc_matrix_method(K, fixed_dofs)
    % 使用对角罚函数法强制应用动态矩阵边界约束
    n = size(K, 1);
    fixed_idx = fixed_dofs;
    diag_01 = ones(n, 1);
    diag_01(fixed_idx) = 0;
    
    D = spdiags(diag_01, 0, n, n);
    I_minus_D = spdiags(1-diag_01, 0, n, n);
    K_mod = D * K * D + I_minus_D;
    F_mod = []; 
end
%% 后处理==========================================================================
function res = extract_results_flexible(U, dof, mesh, scl)
    % 分配并还原缩放的各类物理响应到空间节点字典映射
    n_nodes = mesh.n_nodes;
    res.u     = nan(n_nodes, 1);
    res.v     = nan(n_nodes, 1);
    res.phi   = nan(n_nodes, 1);
    res.psi11 = nan(n_nodes, 1);
    res.psi22 = nan(n_nodes, 1);
    res.psi12 = nan(n_nodes, 1);
    
    mask_u = (dof.map(:, 1) > 0); 
    if any(mask_u)
        res.u(mask_u) = U(dof.map(mask_u, 1));
        res.v(mask_u) = U(dof.map(mask_u, 2));
    end
    
    mask_phi = (dof.map(:, 3) > 0); 
    if any(mask_phi)
        res.phi(mask_phi) = U(dof.map(mask_phi, 3))* scl.S_phi;
        res.psi11(mask_phi) = U(dof.map(mask_phi, 4))* scl.S_psi;
        res.psi22(mask_phi) = U(dof.map(mask_phi, 5))* scl.S_psi;
        res.psi12(mask_phi) = U(dof.map(mask_phi, 6))* scl.S_psi;
    end
    
    if all(isnan(res.u)), warning('提取结果为空. 请检查 dof.map'); end
end

function result_matrix = export_full_results(mesh, results, filename)
    % 将所有响应参量组装至结构化数组并输出为文件存储格式
    x = mesh.coords(:, 1);
    y = mesh.coords(:, 2);
    result_matrix = [x, y, results.u, results.v, results.psi11, results.psi22, results.psi12, results.phi];
    result_matrix = rmmissing(result_matrix);
    fprintf('[Export] 成功生成全场结果矩阵，尺寸: %d x 8\n', size(result_matrix, 1));
    
    if nargin > 2 && ~isempty(filename)
        save(filename, 'result_matrix', '-ascii');
        fprintf('[Export] 矩阵已保存至文件: %s\n', filename);
    end
end

function data = extract_boundary_data(mesh, results, target_tags, var_name, sort_mode)
    % 指定标签路径并依照弧长插值获取曲线结果趋势
    if nargin < 5, sort_mode = 'auto'; end
    if ~isfield(mesh, 'edge_tags'), error('Mesh 缺少 edge_tags'); end
    
    mask = ismember(mesh.edge_tags, target_tags);
    edge_idx = find(mask);
    if isempty(edge_idx), data = []; return; end
    
    current_edges = mesh.edge_elems(edge_idx, :);
    node_indices = unique(current_edges(:));
    xy = mesh.coords(node_indices, :);
    
    x_span = max(xy(:,1)) - min(xy(:,1));
    y_span = max(xy(:,2)) - min(xy(:,2));
    switch lower(sort_mode)
        case 'x', base = 1; case 'y', base = 2;
        case 'auto', if x_span >= y_span, base = 1; else, base = 2; end
    end
    [~, sort_order] = sort(xy(:, base));
    sorted_nodes = node_indices(sort_order);
    sorted_xy = xy(sort_order, :);
    
    d_vec = diff(sorted_xy, 1, 1);
    arc_len = [0; cumsum(sqrt(sum(d_vec.^2, 2)))];
    
    switch lower(var_name)
        case 'u',      vals = results.u(sorted_nodes);
        case 'v',      vals = results.v(sorted_nodes);
        case 'phi',    vals = results.phi(sorted_nodes);
        case {'psi11', 'e11', 'exx', 'psixx'}, vals = results.psi11(sorted_nodes);
        case {'psi22', 'e22', 'eyy', 'psiyy'}, vals = results.psi22(sorted_nodes);
        case {'psi12', 'e12', 'exy', 'psixy'}, vals = results.psi12(sorted_nodes);
        case 'u_mag',  vals = sqrt(results.u(sorted_nodes).^2 + results.v(sorted_nodes).^2);
        otherwise, error('未知变量: %s', var_name);
    end
    
    data.nodes = sorted_nodes;
    data.x = sorted_xy(:,1); data.y = sorted_xy(:,2); data.arc_len = arc_len;
    data.values = vals;
    valid = ~isnan(vals);
    data.clean_x = data.x(valid); data.clean_y = data.y(valid);
    data.clean_arc_len = data.arc_len(valid); data.clean_values = vals(valid);
end

function visualize_results(mesh, res)
    % 生成场图渲染系统全场特性
    hFig = figure('Color','w', 'Name', 'FlexoFEM Results');
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
    title(t, 'Flexoelectric Simulation Results (\psi_{ij} formulation)', 'FontSize', 14);
    colormap(hFig, 'jet');
    
    nexttile;
    max_disp = max(abs(res.v));
    scale = 1.0; if max_disp > 1e-12, scale = 0.1 * mesh.x_max / max_disp; end
    deformed_coords = mesh.coords + scale * [res.u, res.v];
    patch('Faces', mesh.elems(:, 1:4), 'Vertices', deformed_coords, 'FaceVertexCData', res.v, 'FaceColor', 'interp', 'EdgeColor', 'none');
    title(sprintf('Disp v (Scale: %.1f)', scale)); axis equal; axis tight; colorbar;
    
    nexttile; plot_field_patch(mesh, res.phi, 'Electric Potential \phi (V)');
    nexttile; plot_field_patch(mesh, res.psi11, 'Auxiliary Strain \psi_{11}');
    nexttile; plot_field_patch(mesh, res.psi12, 'Auxiliary Strain \psi_{12}');
end

function plot_field_patch(mesh, data, title_str)
    % 辅助画图模块，跳过非几何关联的失效点显示图层
    if all(isnan(data))
        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center'); axis off;
    else
        patch('Faces', mesh.elems(:, 1:4), 'Vertices', mesh.coords, 'FaceVertexCData', data, 'FaceColor', 'interp', 'EdgeColor', 'none');
        title(title_str); axis equal; axis tight; colorbar;
    end
end

function plot_boundary_result(mesh, results, target_tags, var_name, sort_mode, varargin)
    % 高亮呈现沿着边界方向提取的值，增加极点标记辅助读图
    if nargin < 5 || isempty(sort_mode), sort_mode = 'auto'; end
    if nargin < 6 || isempty(varargin), varargin = {'k-', 'LineWidth', 1.5}; end

    if ~isfield(mesh, 'edge_tags'), warning('Mesh 缺少 edge_tags'); return; end
    mask = ismember(mesh.edge_tags, target_tags);
    edge_idx = find(mask);
    
    if isempty(edge_idx)
        fprintf('[Plot] Warning: Tags %s 未找到边界单元。\n', mat2str(target_tags)); return;
    end
    
    current_edges = mesh.edge_elems(edge_idx, :);
    node_indices = unique(current_edges(:));
    xy = mesh.coords(node_indices, :);
    
    x_span = max(xy(:,1)) - min(xy(:,1)); y_span = max(xy(:,2)) - min(xy(:,2));
    switch lower(sort_mode)
        case 'x', base = 1; case 'y', base = 2;
        case {'auto','arc'}, if x_span >= y_span, base = 1; else, base = 2; end
        otherwise, error('Unknown sort_mode');
    end
    
    [~, sort_order] = sort(xy(:, base));
    sorted_nodes = node_indices(sort_order);
    sorted_xy = xy(sort_order, :);
    
    if strcmpi(sort_mode, 'arc')
        d_vec = diff(sorted_xy, 1, 1);
        plot_x = [0; cumsum(sqrt(sum(d_vec.^2, 2)))]; 
        axis_label = 'Arc Length (m)';
    else
        plot_x = sorted_xy(:, base);
        if base==1, axis_label='X Coordinate (m)'; else, axis_label='Y Coordinate (m)'; end
    end
    
    switch lower(var_name)
        case 'u',      vals = results.u(sorted_nodes); y_lab = 'u (m)';
        case 'v',      vals = results.v(sorted_nodes); y_lab = 'v (m)';
        case 'phi',    vals = results.phi(sorted_nodes); y_lab = '\phi (V)';
        case 'u_mag',  vals = sqrt(results.u(sorted_nodes).^2 + results.v(sorted_nodes).^2); y_lab = '|U| (m)';
        case {'psi11', 'e11', 'psixx'}, vals = results.psi11(sorted_nodes); y_lab = '\psi_{11}';
        case {'psi22', 'e22', 'psiyy'}, vals = results.psi22(sorted_nodes); y_lab = '\psi_{22}';
        case {'psi12', 'e12', 'psixy'}, vals = results.psi12(sorted_nodes); y_lab = '\psi_{12}';
        case {"eta111", "grad11x"}    , vals = results.eta111(sorted_nodes); y_lab = '\eta_{111}';
        otherwise, error('变量 %s 未知', var_name);
    end
    
    valid_mask = ~isnan(vals);
    if sum(valid_mask) < 2, return; end
    
    x_clean = plot_x(valid_mask); y_clean = vals(valid_mask);
    figure('Color','w', 'Name', sprintf('Boundary %s', mat2str(target_tags)));
    plot(x_clean, y_clean, varargin{:}); hold on; grid on; box on;
    
    [max_v, max_idx] = max(y_clean); [min_v, min_idx] = min(y_clean);
    plot(x_clean(max_idx), max_v, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    plot(x_clean(min_idx), min_v, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    
    text(x_clean(max_idx), max_v, sprintf(' Max: %.4e', max_v), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9);
    text(x_clean(min_idx), min_v, sprintf(' Min: %.4e', min_v), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', 'b', 'FontWeight', 'bold', 'FontSize', 9);
    
    xlabel(axis_label); ylabel(y_lab);
    title(sprintf('Boundary Tags %s | %s', mat2str(target_tags), y_lab));
    hold off;
end

function err_rel = compute_boundary_error(mesh, results, target_tags, var_name, exact_func)
    % 执行给定解析解函数并在 L2 范数下量化数值逼近误差
    if ~isfield(mesh, 'edge_tags'), error('Mesh 缺少 edge_tags'); end
    mask = ismember(mesh.edge_tags, target_tags);
    edge_indices = find(mask);
    
    if isempty(edge_indices), err_rel = NaN; return; end
    
    edges = mesh.edge_elems(edge_indices, :);
    n_edges = size(edges, 1);
    
    g_pt = [-0.7745966692, 0.0, 0.7745966692];
    g_wt = [0.5555555556, 0.8888888889, 0.5555555556];
    
    numerator_sq = 0.0; denominator_sq = 0.0; 
    fprintf('[Error] 正在计算边界 Tags %s 上的 L2 误差...\n', mat2str(target_tags));
    
    for e = 1:n_edges
        nodes = edges(e, :);      
        xy_nodes = mesh.coords(nodes, :); 
        
        switch lower(var_name)
            case 'u', vals_elem = results.u(nodes);
            case 'v', vals_elem = results.v(nodes);
            case 'phi', vals_elem = results.phi(nodes([1,2])); is_linear = true;
            case {'psi11','psixx'}, vals_elem = results.psi11(nodes([1,2])); is_linear = true;
            otherwise
                if isfield(results, var_name)
                    vals_elem = results.(var_name)(nodes); is_linear = false;
                else
                    error('变量 %s 未知', var_name);
                end
        end
        
        if strcmpi(var_name, 'u') || strcmpi(var_name, 'v'), is_linear = false; end
        
        for k = 1:3
            xi = g_pt(k); w = g_wt(k);
            [N_geom, dN_dxi] = shape_1D_quad(xi);
            
            x_g = N_geom * xy_nodes(:,1); y_g = N_geom * xy_nodes(:,2);
            dxdxi = dN_dxi * xy_nodes(:,1); dydxi = dN_dxi * xy_nodes(:,2);
            J_line = sqrt(dxdxi^2 + dydxi^2);
            
            if is_linear
                N_field = shape_1D_linear(xi); u_h = N_field * vals_elem;
            else
                u_h = N_geom * vals_elem;
            end
            
            u_ex = exact_func(x_g, y_g);
            numerator_sq   = numerator_sq   + (u_h - u_ex)^2 * J_line * w;
            denominator_sq = denominator_sq + (u_ex)^2       * J_line * w;
        end
    end
    
    if denominator_sq < 1e-20
        err_rel = sqrt(numerator_sq);
    else
        err_rel = sqrt(numerator_sq) / sqrt(denominator_sq);
    end
    fprintf('     L2 Relative Error (%s): %.4e (%.2f%%)\n', var_name, err_rel, err_rel*100);
end

%%========================================================================
%  MODULE: 导出场量计算 (应变、电场、应变梯度)
%  ========================================================================
function fields = compute_derived_fields(mesh, dof, U_global, conf, info, scl)
% COMPUTE_DERIVED_FIELDS 计算位移产生的应变、电场和应变梯度
%
% 输入:
%   mesh, dof, U_global, conf, info
%   scl:  缩放因子结构体 (包含 .S_phi, .S_psi)
%
% 输出: fields 结构体, 包含:
%   .strain_u : [n_nodes x 3] (exx, eyy, gam_xy)
%   .E_field  : [n_nodes x 2] (Ex, Ey)
%   .grad_psi : [n_nodes x 6] (dpsi11_dx, dpsi22_dx, ...)

    fprintf('[Post] Computing derived fields (Strain, E-field, Gradients) for %s...\n', conf.el_type);

    n_nodes = mesh.n_nodes;
    
    % --- 初始化累加器 (用于节点平均) ---
    sum_strain = zeros(n_nodes, 3); 
    cnt_strain = zeros(n_nodes, 1);
    
    sum_E = zeros(n_nodes, 2);
    sum_GP = zeros(n_nodes, 6); 
    cnt_Q4 = zeros(n_nodes, 1);
    
    % 局部坐标点 (包含所有9个可能的节点位置)
    local_coords = [-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0; 0 0]; 
    local_Q4 = [-1 -1; 1 -1; 1 1; -1 1]; % 仅角点
    
    % 判断位移场单元类型
    is_Q9 = strcmpi(conf.el_type, 'Q9');
    u_node_indices = info.idx_u_nodes; % Q9 为 1:9, Q8 为 1:8
    n_u_nodes = length(u_node_indices);
    
    % --- 遍历单元 ---
    for e = 1:mesh.n_elems
        nodes = mesh.elems(e, :);    % 全部节点
        corners = nodes(1:4);        % 4个角点
        
        % 提取当前单元有效的位移节点及坐标
        u_nodes = nodes(u_node_indices);
        xy_u = mesh.coords(u_nodes, :);
        
        % 提取单元节点解 (位移通常不缩放，S_u = 1.0)
        u_ele = U_global(dof.map(u_nodes, 1));
        v_ele = U_global(dof.map(u_nodes, 2));

        %考虑系统缩放
        phi_ele   = U_global(dof.map(corners, 3)) * scl.S_phi;
        psi11_ele = U_global(dof.map(corners, 4)) * scl.S_psi;
        psi22_ele = U_global(dof.map(corners, 5)) * scl.S_psi;
        psi12_ele = U_global(dof.map(corners, 6)) * scl.S_psi;
        
        % -------------------------------------------------------------
        % 1. 计算位移应变 (在有效的 u_nodes 上循环)
        % -------------------------------------------------------------
        for i = 1:n_u_nodes
            local_idx = u_node_indices(i);
            xi = local_coords(local_idx, 1); 
            eta = local_coords(local_idx, 2);
            
            % 动态选择形函数
            if is_Q9
                [~, dNu_dxi] = shape_Q9(xi, eta);
            else
                [~, dNu_dxi] = shape_Q8(xi, eta);
            end
            
            % 雅可比矩阵
            J = dNu_dxi * xy_u; 
            dNu_dx = J \ dNu_dxi; % [dN/dx; dN/dy]
            
            % Strain = [du/dx, dv/dy, du/dy + dv/dx]
            exx = dNu_dx(1,:) * u_ele;
            eyy = dNu_dx(2,:) * v_ele;
            gxy = dNu_dx(2,:) * u_ele + dNu_dx(1,:) * v_ele;
            
            global_idx = u_nodes(i);
            sum_strain(global_idx, :) = sum_strain(global_idx, :) + [exx, eyy, gxy];
            cnt_strain(global_idx) = cnt_strain(global_idx) + 1;
        end
        
        % -------------------------------------------------------------
        % 2. 计算电场和应变梯度 (仅在 4 个角点 Q4)
        % -------------------------------------------------------------
        for i = 1:4
            xi = local_Q4(i, 1); 
            eta = local_Q4(i, 2);
            
            [~, dN4_dxi] = shape_Q4(xi, eta);
            
            % 使用当前配置的等参映射计算几何雅可比
            if is_Q9
                [~, dNu_dxi_c] = shape_Q9(xi, eta);
            else
                [~, dNu_dxi_c] = shape_Q8(xi, eta);
            end
            J_geom = dNu_dxi_c * xy_u;
            dN4_dx = J_geom \ dN4_dxi;
            
            % E = - Grad(phi)  (此时 phi_ele 已是真实电压，求导结果为真实电场 V/m)
            Ex = - (dN4_dx(1,:) * phi_ele);
            Ey = - (dN4_dx(2,:) * phi_ele);
            
            % Grad(Psi) = [d11/dx, d22/dx, d12/dx, d11/dy, d22/dy, d12/dy]
            % (此时 psi_ele 已是真实辅助应变，求导结果为真实应变梯度 1/m)
            gp = zeros(1, 6);
            gp(1) = dN4_dx(1,:) * psi11_ele; 
            gp(2) = dN4_dx(1,:) * psi22_ele; 
            gp(3) = dN4_dx(1,:) * psi12_ele; 
            gp(4) = dN4_dx(2,:) * psi11_ele; 
            gp(5) = dN4_dx(2,:) * psi22_ele; 
            gp(6) = dN4_dx(2,:) * psi12_ele; 
            
            global_idx = corners(i);
            sum_E(global_idx, :) = sum_E(global_idx, :) + [Ex, Ey];
            sum_GP(global_idx, :) = sum_GP(global_idx, :) + gp;
            cnt_Q4(global_idx) = cnt_Q4(global_idx) + 1;
        end
    end
    
    % --- 节点平均 ---
    fields.strain_u = sum_strain ./ max(1, cnt_strain);
    fields.E_field  = sum_E ./ max(1, cnt_Q4);
    fields.grad_psi = sum_GP ./ max(1, cnt_Q4);
    
    % 清理未分配物理量的节点 (例如 Q8 的中心点或 Q4 的中点)
    % 设为 NaN 以便绘图函数 (Patch/Surf) 自动过滤
    mask_Q4 = (cnt_Q4 == 0);
    fields.E_field(mask_Q4, :) = NaN;
    fields.grad_psi(mask_Q4, :) = NaN;
    
    mask_strain = (cnt_strain == 0);
    fields.strain_u(mask_strain, :) = NaN;
end

function [time_vec, avg_curve, max_curve, raw_data] = extract_history_by_tag(history, mesh, dof, tag, var_name)
    % 为响应绘制提取瞬态模拟中时间步的步态变量
    nodes = get_nodes_by_tag(mesh, tag);
    if isempty(nodes), error('Tag %d 未包含任何节点。', tag); end
    
    switch lower(var_name)
        case 'u', col_idx = 1; case 'v', col_idx = 2; case 'phi', col_idx = 3;
        case 'psi11', col_idx = 4; case 'psi22', col_idx = 5; case 'psi12', col_idx = 6;
        otherwise, error('未知的变量名: %s', var_name);
    end
    
    global_dofs = dof.map(nodes, col_idx);
    valid_mask = (global_dofs > 0);
    
    if ~any(valid_mask)
        warning('Tag %d 上没有关于变量 %s 的有效自由度。', tag, var_name);
        n_steps = length(history);
        if n_steps > 0, time_vec = [history.t]'; else, time_vec = []; end
        avg_curve = nan(n_steps, 1); max_curve = nan(n_steps, 1); raw_data = []; return;
    end
    
    target_dofs = global_dofs(valid_mask); n_nodes_valid = length(target_dofs);
    n_steps = length(history);
    time_vec  = zeros(n_steps, 1); avg_curve = zeros(n_steps, 1); max_curve = zeros(n_steps, 1); raw_data  = zeros(n_steps, n_nodes_valid);
    
    for i = 1:n_steps
        time_vec(i) = history(i).t;
        if isfield(history(i), 'u_vec') && ~isempty(history(i).u_vec)
            u_current = history(i).u_vec; vals = u_current(target_dofs);
            raw_data(i, :) = vals'; avg_curve(i) = mean(vals); max_curve(i) = max(abs(vals));
        else
            avg_curve(i) = NaN; max_curve(i) = NaN; raw_data(i, :) = NaN;
        end
    end
end

function [time_vec, Q_curve] = extract_charge_history_by_tag(history, mesh, dof, K_global, tag, scl)
    % 通过恢复并集成边界上反应力解析动态电荷积累量
    fprintf('[Post] 正在计算 Tag %d 上的感应电荷 (Scaling Corrected)...\n', tag);
    if nargin < 6 || isempty(scl), S_phi = 1.0; else, S_phi = scl.S_phi; end

    nodes = get_nodes_by_tag(mesh, tag);
    if isempty(nodes), error('Tag %d 未包含任何节点。', tag); end
    
    phi_dofs = dof.map(nodes, 3); target_dofs = phi_dofs(phi_dofs > 0);
    if isempty(target_dofs), error('Tag %d 上没有有效的电势自由度。', tag); end
    
    K_sub = K_global(target_dofs, :);
    n_steps = length(history); time_vec = zeros(n_steps, 1); Q_curve  = zeros(n_steps, 1);
    
    for i = 1:n_steps
        time_vec(i) = history(i).t;
        if isfield(history(i), 'u_vec') && ~isempty(history(i).u_vec)
            U_curr = history(i).u_vec; 
            reaction_forces_scaled = K_sub * U_curr;
            Q_total_scaled = sum(reaction_forces_scaled);
            Q_curve(i) = Q_total_scaled / S_phi;
        else
            Q_curve(i) = NaN;
        end
    end
    max_Q = max(abs(Q_curve));
    fprintf('     Tag %d Peak Charge = %.4e C (Scaled by 1/%.0e)\n', tag, max_Q, S_phi);
end

%% 基础数学辅助组件 ========================================================================
%  包括组装微分操作阵及各类多项式有限元插值函数
%  ========================================================================

function Bu = build_Bu(dNu_dx, n_nodes)
    Bu = zeros(3, n_nodes*2);
    for k=1:n_nodes
        col=(k-1)*2+[1,2]; dNx=dNu_dx(1,k); dNy=dNu_dx(2,k);
        Bu(:,col)=[dNx,0; 0,dNy; dNy,dNx];
    end
end

function Bphi = build_Bphi(dN4_dx)
    Bphi = zeros(2, 4); for k=1:4, Bphi(:, k) = -dN4_dx(:, k); end
end

function Bpsi = build_Bpsi(N4)
    Bpsi = zeros(3, 12); for k=1:4, col=(k-1)*3+(1:3); Bpsi(:,col) = N4(k)*eye(3); end
end

function B_gpsi = build_B_gpsi(dN4_dx)
    B_gpsi = zeros(6, 12);
    for k=1:4
        col=(k-1)*3+(1:3); dNx=dN4_dx(1,k); dNy=dN4_dx(2,k);
        B_gpsi(1,col(1))=dNx; B_gpsi(2,col(2))=dNx; B_gpsi(3,col(3))=1*dNx;
        B_gpsi(4,col(1))=dNy; B_gpsi(5,col(2))=dNy; B_gpsi(6,col(3))=1*dNy;
    end
end

function B_gpsi_b = build_B_gpsi_bubble(dNb_dx)
    B_gpsi_b = zeros(6, 3); dNx = dNb_dx(1); dNy = dNb_dx(2);
    B_gpsi_b(1,1)=dNx; B_gpsi_b(2,2)=dNx; B_gpsi_b(3,3)=1*dNx;
    B_gpsi_b(4,1)=dNy; B_gpsi_b(5,2)=dNy; B_gpsi_b(6,3)=1*dNy;
end

function [N, dN] = shape_1D_quad(xi)
    N = zeros(1,3); dN = zeros(1,3);
    N(1)=0.5*xi*(xi-1); dN(1)=xi-0.5;
    N(2)=0.5*xi*(xi+1); dN(2)=xi+0.5;
    N(3)=1-xi^2;        dN(3)=-2*xi;
end

function N = shape_1D_linear(xi)
    N = [0.5*(1-xi), 0.5*(1+xi)];
end

function [N, dN] = shape_Q9(xi, eta)
    f = @(c) [0.5*c*(c-1); 1-c^2; 0.5*c*(c+1)];
    df = @(c) [c-0.5; -2*c; c+0.5];
    nx = f(xi); ny = f(eta); dnx = df(xi); dny = df(eta);
    idx_x = [1 3 3 1 2 3 2 1 2]; idx_y = [1 1 3 3 1 2 3 2 2];
    N=zeros(1,9); dN=zeros(2,9);
    for k=1:9
        N(k)=nx(idx_x(k))*ny(idx_y(k)); 
        dN(1,k)=dnx(idx_x(k))*ny(idx_y(k)); 
        dN(2,k)=nx(idx_x(k))*dny(idx_y(k)); 
    end
end

function [N, dN] = shape_Q8(xi, eta)
    N = zeros(1, 8); dN = zeros(2, 8);
    nodes_xi  = [-1,  1,  1, -1]; nodes_eta = [-1, -1,  1,  1];
    
    for i = 1:4
        xi_i = nodes_xi(i); et_i = nodes_eta(i);
        term_xi  = (1 + xi*xi_i); term_eta = (1 + eta*et_i); term_sum = (xi*xi_i + eta*et_i - 1);
        N(i) = 0.25 * term_xi * term_eta * term_sum;
        dN(1,i) = 0.25 * xi_i * term_eta * (2*xi*xi_i + eta*et_i);
        dN(2,i) = 0.25 * et_i * term_xi * (xi*xi_i + 2*eta*et_i);
    end
    
    N(5)    =  0.5 * (1 - xi^2) * (1 - eta); dN(1,5) = -xi  * (1 - eta);       dN(2,5) = -0.5 * (1 - xi^2);
    N(6)    =  0.5 * (1 + xi) * (1 - eta^2); dN(1,6) =  0.5 * (1 - eta^2);      dN(2,6) = -eta * (1 + xi);
    N(7)    =  0.5 * (1 - xi^2) * (1 + eta); dN(1,7) = -xi  * (1 + eta);        dN(2,7) =  0.5 * (1 - xi^2);
    N(8)    =  0.5 * (1 - xi) * (1 - eta^2); dN(1,8) = -0.5 * (1 - eta^2);      dN(2,8) = -eta * (1 - xi);
end

function [N, dN] = shape_Q4(xi, eta)
    N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    dN = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta); -(1-xi),  -(1+xi),  (1+xi),   (1-xi)];
end

function [Nb, dNb] = shape_Bubble(xi, eta)
    Nb = (1 - xi^2) * (1 - eta^2);
    dNb = [-2*xi*(1-eta^2); -2*eta*(1-xi^2)];
end

function nodes = get_nodes_by_tag(mesh, target_tags)
    % 检索匹配特定几何物理边界组的网格拓扑节点集合
    if ~isfield(mesh, 'edge_tags') || ~isfield(mesh, 'edge_elems'), nodes = []; return; end
    
    mask = ismember(mesh.edge_tags, target_tags);
    edge_indices = find(mask);
    if isempty(edge_indices), nodes = []; return; end
    
    relevant_elems = mesh.edge_elems(edge_indices, :);
    nodes = unique(relevant_elems(:));
end