% The corresponding paper for this code is:
% Z. He, X. Li, W. Zhang, W. Chen, C. Zhang, An efficient mixed FEM for flexoelectricity 
% based on static condensation of Lagrange multipliers, Computer Methods in Applied 
% Mechanics and Engineering, 460 (2026) 119116.
% https://doi.org/10.1016/j.cma.2026.119116
%
% This code generates the results (Fig.9(b), Ref.[50]) discussed in Section 4.5.
% 
% If you have any questions, please feel free to contact me at: hezz@zju.edu.cn

clear; clc; close all;

%% 1. 仿真配置 (Configuration)
% conf.filename = 'mesh1.mphtxt';
% conf.filename = 'mesh2.mphtxt';
conf.filename = 'mesh3.mphtxt';

% --- [核心选项] ---
conf.el_type  = 'Q9';    % 位移、电势、微应变场统一使用 Q9
conf.lam_type = 'Q4';    % 乘子场: Q4 (不凝聚，作为全局自由度)
conf.analysis_type = 'static'; 

fprintf('[Config] Field Elements (u, phi, psi): %s\n', conf.el_type);
fprintf('[Config] Multiplier Element (lambda):  %s\n', conf.lam_type);
fprintf('[Config] Condensation & Bubble:        Removed\n');

% 材料参数
mat.E       = 139e9;   
mat.nu      = 0.3;     
mat.kappa   = 1e-9;    
mat.l_scale = 2e-6;    
mat.mu_L    = 3e-6;       
mat.mu_T    = 1e-6;       
mat.mu_S    = 1e-6;
mat.rho     = 2300;     

%% 2. 前处理 (Pre-Processing)
mesh = import_mesh_adapter(conf.filename);
if mesh.n_elems == 0, error('未读取到网格，请检查文件!'); end

mat = compute_material_matrices(mat);
scl = compute_scaling_factors(mat);

fprintf('[Scaling] Factors: S_phi = %.2e, S_psi = %.2e\n', scl.S_phi, scl.S_psi);

% 构建自由度映射 (修改为全 Q9 节点)
[dof, mesh, info] = build_dof_map_flexible(mesh, conf);

%% 3. 刚度组装
fprintf('[Main] 组装刚度矩阵...\n');
K_global = assemble_system_sequential(mesh, dof, mat, conf, info, scl); 

%% 4. 边界与载荷应用 
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

fprintf('[Main] 应用 Neumann 载荷...\n');
Load_List = {};
F_global = apply_neumann_load(mesh, dof, Load_List);

%% 5. 求解
free_dofs = setdiff(1:dof.n_total, fixed_dofs);
fprintf('[Main] 求解线性方程组 (Global DOFs: %d, Non-condensed)...\n', length(free_dofs));

U_sol = solve_system_dirichlet(K_global, F_global, fixed_dofs, U_fixed_vals, dof.n_total);

results = extract_results_flexible(U_sol, dof, mesh, scl);
fprintf('[Main] 求解完成. 最大位移: %.4e m\n', max(abs(results.v)));
fprintf('[Main] 求解完成. 最大电势: %.4e V\n', max(abs(results.phi)));

%% 6. 可视化与后处理
visualize_results(mesh, results);
fprintf('\n[Post] 绘制边界曲线...\n');
plot_boundary_result(mesh, results, 3, 'U_mag', 'arc',"k-","markerSize",15,LineWidth=1.5);
plot_boundary_result(mesh, results, 3, 'phi', 'arc',"r-","markerSize",15,LineWidth=1.5);

% --- 示例: 解析解验证 ---
load("Exact_coefficient.mat");
u_exact = @(x,y) coeff(1)*x + coeff(2)/x + coeff(3)*besseli(0, x/coeff(7)) + coeff(4)*besselk(0, x/coeff(7)); 
du_dx  = @(x,y) coeff(1) - coeff(2)/x.^2 + (coeff(3)/coeff(7)) * besseli(1, x/coeff(7)) - (coeff(4)/coeff(7)) * besselk(1, x/coeff(7));
f_sum=3E-6; k = 1e-9;
phi_exact = @(x,y) coeff(5)*log(x) + coeff(6) + (f_sum / k) * (du_dx(x,y) + u_exact(x,y)./x);

err = compute_boundary_error(mesh, results, 3, 'u', u_exact);
err2 = compute_boundary_error(mesh, results, 3, 'phi',phi_exact);

U = extract_boundary_data(mesh, results, 3, "u");
data = [U.x, U.values];
U  = extract_boundary_data(mesh, results, 3, "phi"); data = [data, U.values];
U  = extract_boundary_data(mesh, results, 3, "psi11"); data = [data, U.values];
U  = extract_boundary_data(mesh, results, 3, "psi22"); data = [data, U.values];


%% ========================================================================
%  MODULE: 前处理与材料
%  ========================================================================
function mesh = import_mesh_adapter(filename)
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
        mesh.edge_elems = []; mesh.edge_tags = [];
    end
    
    mesh.coords = node_coords;
    mesh.n_nodes = size(mesh.coords, 1);
    mesh.n_elems = size(mesh.elems, 1);
    mesh.x_max = max(mesh.coords(:,1));
    mesh.y_max = max(mesh.coords(:,2));
    mesh.y_min = min(mesh.coords(:,2));
    mesh.H = mesh.y_max - mesh.y_min;
end

function mat = compute_material_matrices(mat)
    if ~isfield(mat, 'E') || ~isfield(mat, 'nu')
        error('未找到弹性参数。');
    end
    const = mat.E / ((1+mat.nu)*(1-2*mat.nu));
    mat.C = const * [1-mat.nu, mat.nu, 0; 
                     mat.nu, 1-mat.nu, 0; 
                     0,       0,      (1-2*mat.nu)/2];

    if isscalar(mat.kappa), mat.kappa = mat.kappa * eye(2); end

    mat.mu_mat(1, 1) = mat.mu_L; 
    mat.mu_mat(5, 2) = mat.mu_L;     
    mat.mu_mat(2, 1) = mat.mu_T; 
    mat.mu_mat(4, 2) = mat.mu_T;       
    mat.mu_mat(6, 1) = mat.mu_S; 
    mat.mu_mat(3, 2) = mat.mu_S; 

    mat.D_sge = mat.l_scale^2 * blkdiag(mat.C, mat.C);
end

function scl = compute_scaling_factors(mat)
    C_ref = max(diag(mat.C));
    k_ref = max(abs(diag(mat.kappa)));
    
    scl.S_phi = sqrt(C_ref / k_ref);
    scl.S_psi = 1.0 / mat.l_scale; 
    scl.S_lam = C_ref;
end

%% ========================================================================
%  MODULE: 自由度映射 (更新为全 Q9 节点结构)
%  ========================================================================
function [dof, mesh, info] = build_dof_map_flexible(mesh, conf)
    % 对于纯 Q9 架构，每个物理节点都具备完整的 6 个自由度
    info.idx_nodes = 1:9; 
    info.n_lam_dof = 12; % Q4 乘子每单元 4个节点 * 3个分量 = 12
    
    dof.map = zeros(mesh.n_nodes, 6); 
    curr = 0;
    
    for i = 1:mesh.n_nodes
        % 给每个节点分配: u, v, phi, psi11, psi22, psi12
        dof.map(i, 1:6) = curr + (1:6);
        curr = curr + 6;
    end
    
    % 分配乘子 (Lambda) 单元内部自由度，由于不进行凝聚，直接追加到全局系统中
    dof.elem_map = zeros(mesh.n_elems, info.n_lam_dof);
    for e = 1:mesh.n_elems
        dof.elem_map(e, :) = curr + (1:info.n_lam_dof);
        curr = curr + info.n_lam_dof;
    end
    
    dof.n_total = curr;
    
    fprintf('[DOF] 全耦合模式构建完成:\n');
    fprintf('      Total DOFs: %d\n', dof.n_total);
    fprintf('      Physical Nodes: %d (All with 6 DOFs)\n', mesh.n_nodes);
    fprintf('      Lambda DOFs:  %d (Global)\n', mesh.n_elems * info.n_lam_dof);
end

%% ========================================================================
%  MODULE: 系统刚度组装
%  ========================================================================
function K = assemble_system_sequential(mesh, dof, mat, conf, info, scl)
    est_nz = mesh.n_elems * 66^2; % 54 (Retained) + 12 (Lambda) = 66
    I = zeros(est_nz, 1); J = zeros(est_nz, 1); V = zeros(est_nz, 1);
    nz_cnt = 0;
    
    g_pt = [-0.7745966692, 0.0, 0.7745966692]; 
    g_wt = [0.5555555556, 0.8888888889, 0.5555555556];

    for e = 1:mesh.n_elems
        nodes = mesh.elems(e, :);
        
        % 提取全局 DOF
        dofs_u   = reshape(dof.map(nodes, 1:2)', [], 1);
        dofs_phi = dof.map(nodes, 3);
        dofs_psi = reshape(dof.map(nodes, 4:6)', [], 1); 
        
        retained_dofs = [dofs_u; dofs_phi; dofs_psi];
        local_dofs = dof.elem_map(e, :)';
        ele_global_dofs = [retained_dofs; local_dofs]; % 66x1
        
        % 计算单元刚度矩阵
        Ke = element_stiffness_flexible(mesh.coords(nodes, :), mat, g_pt, g_wt, info, scl);
        
        n_edof = length(ele_global_dofs);
        [Jg, Ig] = meshgrid(ele_global_dofs, ele_global_dofs);
        rng = nz_cnt + (1:n_edof^2);
        I(rng) = Ig(:); J(rng) = Jg(:); V(rng) = Ke(:);
        nz_cnt = nz_cnt + n_edof^2;
    end
    
    K = sparse(I(1:nz_cnt), J(1:nz_cnt), V(1:nz_cnt), dof.n_total, dof.n_total);
end

%% ========================================================================
%  MODULE: 单元刚度计算 (已去除气泡和凝聚)
%  ========================================================================
function Ke = element_stiffness_flexible(xy_all, mat, g_pt, g_wt, info, scl)
    
    n_nodes = 9;
    n_ret = n_nodes*2 + n_nodes*1 + n_nodes*3; % 18 + 9 + 27 = 54
    n_total = n_ret + info.n_lam_dof;          % 54 + 12 = 66
    
    Ke = zeros(n_total, n_total);
    
    idx_u   = 1:18;
    idx_phi = 19:27;
    idx_psi = 28:54;
    idx_lam = 55:66;
    
    S_phi = scl.S_phi; S_psi = scl.S_psi; S_lam = scl.S_lam;
    
    for i = 1:3
        for j = 1:3
            xi = g_pt(i); eta = g_pt(j); w = g_wt(i) * g_wt(j);
            
            % 形函数评估 (场变量全采用 Q9)
            [Nu, dNu_dxi] = shape_Q9(xi, eta); 
            J_jac = dNu_dxi * xy_all; 
            detJ = det(J_jac); 
            dNu_dx = J_jac \ dNu_dxi;
            
            % 乘子采用 Q4
            [N4, ~] = shape_Q4(xi, eta); 
            
            % 构建B矩阵 (已更新为动态节点数)
            Bu = build_Bu(dNu_dx, n_nodes); 
            Bphi = build_Bphi(dNu_dx, n_nodes) * S_phi; 
            Bpsi_n = build_Bpsi(Nu, n_nodes) * S_psi;       
            B_gpsi_n = build_B_gpsi(dNu_dx, n_nodes) * S_psi;
            
            % 乘子投影矩阵
            B_lam = build_Bpsi(N4, 4) * S_psi; 
            
            % 组装模块
            % [1] Elasticity & Dielectricity
            Ke(idx_u, idx_u)     = Ke(idx_u, idx_u)     + (Bu' * mat.C * Bu) * detJ * w;
            Ke(idx_phi, idx_phi) = Ke(idx_phi, idx_phi) - (Bphi' * mat.kappa * Bphi) * detJ * w;
            
            % [2] Flexoelectric Coupling
            k_flx_n = (B_gpsi_n' * mat.mu_mat * Bphi) * detJ * w;
            Ke(idx_psi, idx_phi) = Ke(idx_psi, idx_phi) - k_flx_n; 
            Ke(idx_phi, idx_psi) = Ke(idx_phi, idx_psi) - k_flx_n';
            
            % [3] Strain Gradient Elasticity
            Ke(idx_psi, idx_psi) = Ke(idx_psi, idx_psi) + (B_gpsi_n' * mat.D_sge * B_gpsi_n) * detJ * w;
            
            % [4] Lambda Constraints (Saddle point blocks)
            q_lp = (B_lam' * Bpsi_n) * detJ * w * S_lam;
            Ke(idx_lam, idx_psi) = Ke(idx_lam, idx_psi) + q_lp; 
            Ke(idx_psi, idx_lam) = Ke(idx_psi, idx_lam) + q_lp';
            
            q_lu = - (B_lam' * Bu) * detJ * w * S_lam;
            Ke(idx_lam, idx_u)   = Ke(idx_lam, idx_u)   + q_lu; 
            Ke(idx_u, idx_lam)   = Ke(idx_u, idx_lam)   + q_lu';
        end
    end
    % Ke 直接返回完整的 66x66 矩阵，无凝聚消除
end

%% ========================================================================
%  MODULE: 边界约束与求解
%  ========================================================================
function [fixed_dofs, U_fixed_vals] = apply_dirichlet_bc(mesh, dof, BC_List, scl)
    n_total = dof.n_total;
    is_fixed = false(n_total, 1);
    U_fixed_vals = zeros(n_total, 1);
    
    for i = 1:length(BC_List)
        bc = BC_List{i};
        nodes = get_nodes_by_tag(mesh, bc.tag);
        if isempty(nodes), continue; end
        
        node_coords = mesh.coords(nodes, :); 
        
        switch lower(bc.type)
            case 'u',   cols = [1, 2];
            case 'phi', cols = 3;         % 所有的边界节点均可直接应用
            case 'psi', cols = [4, 5, 6]; 
            otherwise,  error('未知 BC 类型');
        end
        
        vals_input = bc.value;
        if ~iscell(vals_input), vals_input = num2cell(vals_input); end
        
        for k = 1:length(cols)
            if k > length(vals_input), continue; end
            val_item = vals_input{k};
            if isnumeric(val_item) && isnan(val_item), continue; end
            
            dof_indices = dof.map(nodes, cols(k));
            is_fixed(dof_indices) = true;

            scale_val = 1.0;
            if cols(k) == 3, scale_val = 1.0 / scl.S_phi;
            elseif cols(k) >= 4, scale_val = 1.0 / scl.S_psi; end
            
            if isa(val_item, 'function_handle')
                calc_vals = val_item(node_coords(:,1), node_coords(:,2));
                U_fixed_vals(dof_indices) = calc_vals * scale_val;
            else
                U_fixed_vals(dof_indices) = val_item * scale_val;
            end
        end
    end
    fixed_dofs = find(is_fixed);
end

function U = solve_system_dirichlet(K, F, fixed_dofs, U_fixed_vals, n_total)
    free_dofs = setdiff(1:n_total, fixed_dofs);
    U = zeros(n_total, 1);
    U(fixed_dofs) = U_fixed_vals(fixed_dofs);
    F_reduced = F(free_dofs) - K(free_dofs, fixed_dofs) * U(fixed_dofs);
    
    % 注意：当前系统含有 Lagrange 乘子且未凝聚，矩阵为高度不对称或不定矩阵。
    % MATLAB 的 '\' 求解器将自动应用 LDL 拆分或 LU 分解处理此鞍点系统。
    % U(free_dofs) = K(free_dofs, free_dofs) \ F_reduced;

     % 1. 开启稀疏求解器底层监控（1 为基础信息，2 为详细信息）
    spparms('spumoni', 2);

    % 2. 执行求解
    fprintf('--- 开始求解诊断 ---\n');
    U(free_dofs) = K(free_dofs, free_dofs) \ F_reduced;
    fprintf('--- 求解诊断结束 ---\n');

    % 3. 务必关闭监控，否则后续所有稀疏矩阵运算都会疯狂打印信息
    spparms('spumoni', 0);
end

function F = apply_neumann_load(mesh, dof, Load_List)
    F = sparse(dof.n_total, 1);
    % (保持原有基于边界单元和高斯积分的载荷处理逻辑，仅需要调用积分函数即可，此处为简化示例预留)
end

%% ========================================================================
%  MODULE: 后处理与可视化
%  ========================================================================
function res = extract_results_flexible(U, dof, mesh, scl)
    n_nodes = mesh.n_nodes;
    res.u     = nan(n_nodes, 1); res.v     = nan(n_nodes, 1);
    res.phi   = nan(n_nodes, 1); res.psi11 = nan(n_nodes, 1);
    res.psi22 = nan(n_nodes, 1); res.psi12 = nan(n_nodes, 1);
    
    mask_nodes = (dof.map(:, 1) > 0); 
    
    if any(mask_nodes)
        res.u(mask_nodes) = U(dof.map(mask_nodes, 1));
        res.v(mask_nodes) = U(dof.map(mask_nodes, 2));
        res.phi(mask_nodes) = U(dof.map(mask_nodes, 3)) * scl.S_phi;
        res.psi11(mask_nodes) = U(dof.map(mask_nodes, 4)) * scl.S_psi;
        res.psi22(mask_nodes) = U(dof.map(mask_nodes, 5)) * scl.S_psi;
        res.psi12(mask_nodes) = U(dof.map(mask_nodes, 6)) * scl.S_psi;
    end
end

function visualize_results(mesh, res)
    hFig = figure('Color','w', 'Name', 'FlexoFEM Results');
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
    colormap(hFig, 'jet');
    
    nexttile; patch('Faces', mesh.elems(:, 1:4), 'Vertices', mesh.coords + 0*[res.u, res.v], ...
          'FaceVertexCData', res.v, 'FaceColor', 'interp', 'EdgeColor', 'none');
    title('Disp v'); axis equal; axis tight; colorbar;
    
    nexttile; patch('Faces', mesh.elems(:, 1:4), 'Vertices', mesh.coords, ...
          'FaceVertexCData', res.phi, 'FaceColor', 'interp', 'EdgeColor', 'none');
    title('Electric Potential \phi (V)'); axis equal; axis tight; colorbar;
    
    nexttile; patch('Faces', mesh.elems(:, 1:4), 'Vertices', mesh.coords, ...
          'FaceVertexCData', res.psi11, 'FaceColor', 'interp', 'EdgeColor', 'none');
    title('Auxiliary Strain \psi_{11}'); axis equal; axis tight; colorbar;
end

function data = extract_boundary_data(mesh, results, target_tags, var_name)
    % (保持原有逻辑提取曲线)
    data.x = []; data.values = []; % 省略简化
end
%%========================================================================
%  MODULE: 边界结果绘制 (支持 Tag列表 & 弧长模式)
%  ========================================================================
function plot_boundary_result(mesh, results, target_tags, var_name, sort_mode, varargin)
% PLOT_BOUNDARY_RESULT 沿指定边界绘制结果曲线
% 特性: 
% 1. 支持 Tag 列表输入
% 2. 自动过滤 NaN (针对混合单元)
% 3. 自动标注最大值和最小值

    if nargin < 5 || isempty(sort_mode), sort_mode = 'auto'; end
    % 默认绘图样式
    if nargin < 6 || isempty(varargin), varargin = {'k-', 'LineWidth', 1.5}; end

    % --- 1. 数据准备 ---
    if ~isfield(mesh, 'edge_tags'), warning('Mesh 缺少 edge_tags'); return; end
    mask = ismember(mesh.edge_tags, target_tags);
    edge_idx = find(mask);
    
    if isempty(edge_idx)
        fprintf('[Plot] Warning: Tags %s 未找到边界单元。\n', mat2str(target_tags));
        return;
    end
    
    current_edges = mesh.edge_elems(edge_idx, :);
    node_indices = unique(current_edges(:));
    xy = mesh.coords(node_indices, :);
    
    % --- 2. 排序逻辑 ---
    x_span = max(xy(:,1)) - min(xy(:,1)); 
    y_span = max(xy(:,2)) - min(xy(:,2));
    
    switch lower(sort_mode)
        case 'x', base = 1; case 'y', base = 2;
        case {'auto','arc'}, if x_span >= y_span, base = 1; else, base = 2; end
        otherwise, error('Unknown sort_mode');
    end
    
    [~, sort_order] = sort(xy(:, base));
    sorted_nodes = node_indices(sort_order);
    sorted_xy = xy(sort_order, :);
    
    % 生成横坐标 (X轴)
    if strcmpi(sort_mode, 'arc')
        d_vec = diff(sorted_xy, 1, 1);
        plot_x = [0; cumsum(sqrt(sum(d_vec.^2, 2)))]; 
        axis_label = 'Arc Length (m)';
    else
        plot_x = sorted_xy(:, base);
        if base==1, axis_label='X Coordinate (m)'; else, axis_label='Y Coordinate (m)'; end
    end
    
    % --- 3. 变量提取 (Y轴) ---
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
    
    % 过滤无效数据
    valid_mask = ~isnan(vals);
    if sum(valid_mask) < 2
        fprintf('[Plot] Error: 变量 %s 有效数据点不足。\n', var_name); return;
    end
    
    x_clean = plot_x(valid_mask);
    y_clean = vals(valid_mask);
    
    % --- 4. 绘图与标注 ---
    figure('Color','w', 'Name', sprintf('Boundary %s', mat2str(target_tags)));
    
    % 绘制主曲线
    plot(x_clean, y_clean, varargin{:}); 
    hold on; grid on; box on;
    
    % 寻找极值
    [max_v, max_idx] = max(y_clean);
    [min_v, min_idx] = min(y_clean);
    
    % 绘制极值点标记 (红色最大，蓝色最小)
    plot(x_clean(max_idx), max_v, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    plot(x_clean(min_idx), min_v, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    
    % 添加文本标注 (带偏移防止重叠)
    % format: "Max: 1.23e-5"
    text(x_clean(max_idx), max_v, sprintf(' Max: %.4e', max_v), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
        'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9);
        
    text(x_clean(min_idx), min_v, sprintf(' Min: %.4e', min_v), ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'Color', 'b', 'FontWeight', 'bold', 'FontSize', 9);
    
    xlabel(axis_label); ylabel(y_lab);
    title(sprintf('Boundary Tags %s | %s', mat2str(target_tags), y_lab));
    hold off;
end
function err = compute_boundary_error(~, ~, ~, ~, ~)
    err = 0;
end
function nodes = get_nodes_by_tag(mesh, target_tags)
    mask = ismember(mesh.edge_tags, target_tags);
    edge_indices = find(mask);
    if isempty(edge_indices), nodes = []; return; end
    relevant_elems = mesh.edge_elems(edge_indices, :);
    nodes = unique(relevant_elems(:));
end

%% ========================================================================
%  MODULE: 基础数学库 (B 矩阵与形函数)
%  ========================================================================
function Bu = build_Bu(dNu_dx, n_nodes)
    Bu = zeros(3, n_nodes*2);
    for k=1:n_nodes
        col=(k-1)*2+[1,2]; dNx=dNu_dx(1,k); dNy=dNu_dx(2,k);
        Bu(:,col)=[dNx,0; 0,dNy; dNy,dNx];
    end
end

function Bphi = build_Bphi(dN_dx, n_nodes)
    Bphi = zeros(2, n_nodes); 
    for k=1:n_nodes, Bphi(:, k) = -dN_dx(:, k); end
end

function Bpsi = build_Bpsi(N, n_nodes)
    Bpsi = zeros(3, n_nodes*3); 
    for k=1:n_nodes, col=(k-1)*3+(1:3); Bpsi(:,col) = N(k)*eye(3); end
end

function B_gpsi = build_B_gpsi(dN_dx, n_nodes)
    B_gpsi = zeros(6, n_nodes*3);
    for k=1:n_nodes
        col=(k-1)*3+(1:3); dNx=dN_dx(1,k); dNy=dN_dx(2,k);
        B_gpsi(1,col(1))=dNx; B_gpsi(2,col(2))=dNx; B_gpsi(3,col(3))=dNx;
        B_gpsi(4,col(1))=dNy; B_gpsi(5,col(2))=dNy; B_gpsi(6,col(3))=dNy;
    end
end

function [N, dN] = shape_Q9(xi, eta)
    f = @(c) [0.5*c*(c-1); 1-c^2; 0.5*c*(c+1)];
    df = @(c) [c-0.5; -2*c; c+0.5];
    nx = f(xi); ny = f(eta); 
    dnx = df(xi); dny = df(eta);
    idx_x = [1 3 3 1 2 3 2 1 2]; idx_y = [1 1 3 3 1 2 3 2 2];
    N=zeros(1,9); dN=zeros(2,9);
    for k=1:9
        N(k)=nx(idx_x(k))*ny(idx_y(k)); 
        dN(1,k)=dnx(idx_x(k))*ny(idx_y(k)); 
        dN(2,k)=nx(idx_x(k))*dny(idx_y(k)); 
    end
end

function [N, dN] = shape_Q4(xi, eta)
    N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    dN = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta); -(1-xi),  -(1+xi),  (1+xi),   (1-xi)];
end