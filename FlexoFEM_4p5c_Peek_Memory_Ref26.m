% The corresponding paper for this code is:
% Z. He, X. Li, W. Zhang, W. Chen, C. Zhang, An efficient mixed FEM for flexoelectricity 
% based on static condensation of Lagrange multipliers, Computer Methods in Applied 
% Mechanics and Engineering, 460 (2026) 119116.
% https://doi.org/10.1016/j.cma.2026.119116
%
% This code generates the results (Fig.9(b), Ref.[26]) discussed in Section 4.5.
% 
% If you have any questions, please feel free to contact me at: hezz@zju.edu.cn

clear; clc; close all;

%% 1. 前处理：读取与修正网格
% conf.filename = 'mesh1.mphtxt';
% conf.filename = 'mesh2.mphtxt';
conf.filename = 'mesh3.mphtxt';

% 材料参数
mat.E       = 139e9;   
mat.nu      = 0.3;     
mat.kappa   = 1e-9;    
mat.l_scale = 2e-6;    

mat.mu_L = 3e-6;       
mat.mu_T = 1e-6;       
mat.mu_S = 1e-6;

mat.rho = 2300; 

mesh = import_mesh_adapter(conf.filename);
mat = compute_material_matrices(mat);
scl = compute_scaling_factors(mat);

%% 2. 自由度映射 (封装调用)
[dof, mesh] = build_dof_map(mesh);

%% 3. 刚度矩阵组装 (封装调用)
[K_global, F_global] = assemble_system(mesh, dof, mat, scl);

%% 4. 施加边界条件与载荷
BC_List = {}; 

% 示例边界条件 (根据实际情况调整)
BC_List{end+1} = struct('tag', [4, 6], 'type', 'u', 'value', {{@(x,y) x./(x.^2+y.^2).^0.5*0.045E-6, @(x,y) y./(x.^2+y.^2).^0.5*0.045E-6}});
BC_List{end+1} = struct('tag', [5, 7], 'type', 'u', 'value', {{@(x,y) x./(x.^2+y.^2).^0.5*0.05E-6, @(x,y) y./(x.^2+y.^2).^0.5*0.05E-6}});
BC_List{end+1} = struct('tag', 3, 'type', 'u', 'value', {{NaN, 0}});
BC_List{end+1} = struct('tag', 1, 'type', 'u', 'value', {{0, NaN}});
BC_List{end+1} = struct('tag', [4, 6], 'type', 'phi', 'value', 0);
BC_List{end+1} = struct('tag', [5, 7], 'type', 'phi', 'value', 1);
BC_List{end+1} = struct('tag', [1, 3], 'type', 'psi',   'value', {{NaN, 0, 0,NaN}}); 

[fixed_dofs, U_fixed_vals] = apply_dirichlet_bc(mesh, dof, BC_List, scl);

%% 5. 求解
free_dofs = setdiff(1:dof.n_total, fixed_dofs);
fprintf('[Main] 求解线性方程组 (DOF: %d)...\n', length(free_dofs));

% 如果有体载荷，F_global 在 assemble_system 中会被填充，这里仅需处理边界
U_sol = solve_system_dirichlet(K_global, F_global, fixed_dofs, U_fixed_vals, dof.n_total);

fprintf('求解完成. 最大位移: %.2e m\n', max(abs(U_sol(dof.map(:,2)))));
fprintf('求解完成. 最大电势: %.2e V\n', max(abs(U_sol(dof.map(dof.map(:,3)>0,3)) * scl.S_phi)));

%% 6. 可视化
results = extract_results_flexible(U_sol, dof, mesh, conf, scl);
visualize_results(mesh, results);

plot_boundary_result(mesh, results, 3, 'U_mag', 'arc',"k-","markerSize",15,LineWidth=1.5);
plot_boundary_result(mesh, results, 3, 'phi', 'arc',"r-","markerSize",15,LineWidth=1.5);

% =========================================================================
%                               核心函数封装
% =========================================================================

function [dof, mesh] = build_dof_map(mesh)
% BUILD_DOF_MAP 构建混合单元的自由度映射
% 输入:
%   mesh: 网格结构体
% 输出:
%   dof: 包含 map (n_nodes x 11) 和 n_total
%   mesh: 更新后的 mesh (增加了 is_corner 标记)

    fprintf('[DOF] 构建自由度映射...\n');
    dof_map = zeros(mesh.n_nodes, 11); 
    current_dof = 0;

    % 识别角点 (Q9 单元的前 4 个节点为角点)
    is_corner = false(mesh.n_nodes, 1);
    unique_corners = unique(mesh.elems(:, 1:4));
    is_corner(unique_corners) = true;
    
    % 将 is_corner 存回 mesh，供后续 BC 或后处理使用
    mesh.is_corner = is_corner;

    for i = 1:mesh.n_nodes
        % 所有节点都有位移 u, v (Lagrange Q9)
        dof_map(i, 1:2) = current_dof + [1, 2]; 
        current_dof = current_dof + 2;
        
        % 仅角点有混合变量 (Lagrange Q4)
        if is_corner(i)
            dof_map(i, 3) = current_dof + 1;       % Phi (电势)
            dof_map(i, 4:7) = current_dof + (2:5); % Psi (辅助应变梯度)
            dof_map(i, 8:11) = current_dof + (6:9);% Lambda (拉格朗日乘子)
            current_dof = current_dof + 9; 
        end
    end
    
    dof.map = dof_map;
    dof.n_total = current_dof;
    fprintf('[DOF] 总自由度数: %d\n', dof.n_total);
end

function [K_global, F_global] = assemble_system(mesh, dof, mat, scl)
% ASSEMBLE_SYSTEM 组装全局刚度矩阵和载荷向量
% 输入:
%   mesh: 网格信息
%   dof:  自由度映射
%   mat:  材料属性矩阵
% 输出:
%   K_global: 稀疏刚度矩阵
%   F_global: 稀疏载荷向量 (目前仅初始化为0，可在此添加体载荷)

    fprintf('[Assembly] 开始组装矩阵 (Using full C_mat for SGE)...\n');
    
    % 预分配稀疏矩阵的三元组
    % 估算非零元: 每个单元约 54*54 个非零元
    est_nz = mesh.n_elems * 54^2; 
    I_idx = zeros(est_nz, 1); 
    J_idx = zeros(est_nz, 1); 
    V_val = zeros(est_nz, 1);
    nz_count = 0;
    
    % 全局载荷向量
    F_global = sparse(dof.n_total, 1);

    % 高斯积分点 (3x3 Gauss Rule)
    g_pt = [-0.7745966692, 0.0, 0.7745966692];
    g_wt = [0.5555555556, 0.8888888889, 0.5555555556];

    % --- 缩放因子提取 ---
    S_phi = scl.S_phi;
    S_psi = scl.S_psi;
    S_lam = scl.S_lam;




    for e = 1:mesh.n_elems
        nodes = mesh.elems(e, :);
        corners = nodes(1:4); % Q4 角点用于混合变量
        
        % 获取当前单元所有相关 DOF 索引
        dofs_u   = reshape(dof.map(nodes, 1:2)', [], 1);
        dofs_phi = dof.map(corners, 3);
        dofs_psi = reshape(dof.map(corners, 4:7)', [], 1);
        dofs_lam = reshape(dof.map(corners, 8:11)', [], 1);
        
        ele_dof_indices = [dofs_u; dofs_phi; dofs_psi; dofs_lam];
        
        % 局部矩阵索引范围
        idx_u=1:18; idx_phi=19:22; idx_psi=23:38; idx_lam=39:54;
        
        Ke = zeros(54, 54);
        xy = mesh.coords(nodes, :); % 单元节点坐标
        
        for i = 1:3
            for j = 1:3
                xi = g_pt(i); eta = g_pt(j);
                w = g_wt(i) * g_wt(j);
                
                % 计算形函数及其导数
                [Nu, dNu_dxi] = shape_Q9(xi, eta);
                J_jac = dNu_dxi * xy; 
                detJ = det(J_jac);
                dNu_dx = J_jac \ dNu_dxi;
                
                [N4, dN4_dxi] = shape_Q4(xi, eta);
                dN4_dx = J_jac \ dN4_dxi; 
                
                % 1. Elasticity K_uu (18x18)
                Bu = zeros(3, 18);
                for k=1:9
                    col=(k-1)*2+[1,2]; 
                    Bu(:,col)=[dNu_dx(1,k),0; 0,dNu_dx(2,k); dNu_dx(2,k),dNu_dx(1,k)]; 
                end
                Ke(idx_u, idx_u) = Ke(idx_u, idx_u) + (Bu' * mat.C * Bu) * detJ * w;
                
                % 2. Dielectric K_pp (4x4)
                Bphi = zeros(2, 4);
                for k=1:4, Bphi(:, k) = -dN4_dx(:, k); end
                Bphi = Bphi*S_phi;
                Ke(idx_phi, idx_phi) = Ke(idx_phi, idx_phi) - (Bphi' * mat.kappa * Bphi) * detJ * w;
                
                % 3. Flexoelectricity & SGE (Coupling & Gradient)
                % Bpsi: 用于约束项 (N4)
                Bpsi = zeros(4, 16);
                for k=1:4, col=(k-1)*4+(1:4); Bpsi(:,col) = N4(k)*eye(4); end
                Bpsi = Bpsi*S_psi;
                % B_gpsi: 用于挠曲电和 SGE 的梯度算子
                B_gpsi = zeros(6, 16);
                for k=1:4
                    col=(k-1)*4+(1:4); dNx=dN4_dx(1,k); dNy=dN4_dx(2,k);
                    % 排列对应 mat.mu_mat 的行序
                    B_gpsi(1,col(1))=dNx; B_gpsi(2,col(4))=dNx; 
                    B_gpsi(3,col(2))=dNx; B_gpsi(3,col(3))=dNx; 
                    B_gpsi(4,col(1))=dNy; B_gpsi(5,col(4))=dNy; 
                    B_gpsi(6,col(2))=dNy; B_gpsi(6,col(3))=dNy; 
                end
                B_gpsi = B_gpsi*S_psi;
                % K_flexo (Couple Psi and Phi)
                k_flexo = (B_gpsi' * mat.mu_mat * Bphi) * detJ * w;
                Ke(idx_psi, idx_phi) = Ke(idx_psi, idx_phi) - k_flexo;
                Ke(idx_phi, idx_psi) = Ke(idx_phi, idx_psi) - k_flexo';
                
                % K_sge (Stabilization/Gradient Elasticity)
                Ke(idx_psi, idx_psi) = Ke(idx_psi, idx_psi) + (B_gpsi' * mat.D_sge * B_gpsi) * detJ * w;
                
                % 4. Constraints (Lagrange Multipliers)
                % B_lam = Bpsi; % Lambda 使用与 Psi 相同的形函数
                B_lam = zeros(4, 16);
                for k=1:4, col=(k-1)*4+(1:4); B_lam(:,col) = N4(k)*eye(4); end

                Bg = zeros(4, 18); % u 的梯度，用于约束
                for k=1:9
                    col=(k-1)*2+[1,2]; 
                    Bg(:,col)=[dNu_dx(1,k),0; dNu_dx(2,k),0; 0,dNu_dx(1,k); 0,dNu_dx(2,k)]; 
                end
                
                % Psi - Lambda 耦合
                q_lp = (B_lam' * eye(4) * Bpsi) * detJ * w*S_lam;
                % u - Lambda 耦合 (弱形式约束: Psi = grad(u))
                q_lu = - (B_lam' * eye(4) * Bg) * detJ * w*S_lam;
                
                Ke(idx_lam, idx_psi) = Ke(idx_lam, idx_psi) + q_lp;
                Ke(idx_psi, idx_lam) = Ke(idx_psi, idx_lam) + q_lp';
                Ke(idx_lam, idx_u)   = Ke(idx_lam, idx_u) + q_lu;
                Ke(idx_u, idx_lam)   = Ke(idx_u, idx_lam) + q_lu';
            end
        end
        
        % 填入三元组
        [Jg, Ig] = meshgrid(ele_dof_indices, ele_dof_indices);
        len = 54^2; 
        rng = nz_count+(1:len);
        I_idx(rng) = Ig(:); 
        J_idx(rng) = Jg(:); 
        V_val(rng) = Ke(:); 
        nz_count = nz_count + len;
    end
    
    % 创建稀疏矩阵
    K_global = sparse(I_idx(1:nz_count), J_idx(1:nz_count), V_val(1:nz_count), dof.n_total, dof.n_total);
end

% =========================================================================
%                               辅助函数 (保持不变)
% =========================================================================

function [coords, elems] = read_comsol_corrected(filename)
    fid = fopen(filename, 'r'); coords = []; elems_raw = [];
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, '# number of mesh vertices')
            n = sscanf(line, '%d');
            while ~contains(line, 'coordinates'), line = fgetl(fid); end
            coords = fscanf(fid, '%f', [2, n])';
        end
        if contains(line, 'quad2') 
            while ~contains(line, '# number of elements'), line = fgetl(fid); end
            n = sscanf(line, '%d');
            while ~contains(line, '# Elements'), line = fgetl(fid); end
            elems_raw = fscanf(fid, '%d', [9, n])' + 1;
        end
    end
    fclose(fid);
    p = [1, 2, 4, 3, 5, 8, 9, 6, 7]; elems = elems_raw(:, p);
end

function [N, dN] = shape_Q9(xi, eta)
    f = @(c) [0.5*c*(c-1); 1-c^2; 0.5*c*(c+1)]; df = @(c) [c-0.5; -2*c; c+0.5];
    nx = f(xi); ny = f(eta); dnx = df(xi); dny = df(eta);
    idx_x = [1 3 3 1 2 3 2 1 2]; idx_y = [1 1 3 3 1 2 3 2 2];
    N=zeros(1,9); dN=zeros(2,9);
    for k=1:9, N(k)=nx(idx_x(k))*ny(idx_y(k)); dN(1,k)=dnx(idx_x(k))*ny(idx_y(k)); dN(2,k)=nx(idx_x(k))*dny(idx_y(k)); end
end

function [N, dN] = shape_Q4(xi, eta)
    N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    dN = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta); -(1-xi),  -(1+xi),  (1+xi),   (1-xi)];
end

function mesh = import_mesh_adapter(filename)
    [node_coords, elem_struct, boundaries] = import_mesh_mixed_v3(filename);
    if isempty(elem_struct.quad2), error('网格文件中未找到 quad2 (Q9) 单元！'); end
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
    fprintf('[Mesh Adapter] 已适配: %d 节点, %d 单元\n', mesh.n_nodes, mesh.n_elems);
end

function mat = compute_material_matrices(mat)
    if isfield(mat, 'C') && ~isempty(mat.C)
        if ~all(size(mat.C) == [3, 3]), error('C must be 3x3'); end
    else
        const = mat.E / ((1+mat.nu)*(1-2*mat.nu));
        mat.C = const * [1-mat.nu, mat.nu, 0; mat.nu, 1-mat.nu, 0; 0, 0, (1-2*mat.nu)/2];
    end
    if isscalar(mat.kappa), mat.kappa = mat.kappa * eye(2); end
    if ~isfield(mat, 'mu_mat') || isempty(mat.mu_mat)
        mat.mu_mat = zeros(6,2);
        mat.mu_mat(1, 1) = mat.mu_L; mat.mu_mat(5, 2) = mat.mu_L; 
        mat.mu_mat(2, 1) = mat.mu_T; mat.mu_mat(4, 2) = mat.mu_T; 
        mat.mu_mat(6, 1) = mat.mu_S; mat.mu_mat(3, 2) = mat.mu_S; 
    end
    mat.D_sge = mat.l_scale^2 * blkdiag(mat.C, mat.C);
end

function scl = compute_scaling_factors(mat)
    % 目标: 让介电项 (kappa) 的量级提升到与弹性项 (E) 相当
    % 取材料属性的特征值
    C_ref = max(diag(mat.C));
    
    % 提取介电常数的特征值 (取范数或最大对角元)
    if isscalar(mat.kappa)
        k_ref = mat.kappa;
    else
        k_ref = max(abs(diag(mat.kappa)));
    end
    
    % 1. 电势缩放因子 S_phi
    % 使得: S_phi^2 * kappa ~ E
    scl.S_phi = sqrt(C_ref / k_ref);
    
    % 2. 辅助应变缩放因子 S_psi (可选)
    % 使得: S_psi^2 * D_sge ~ E  => S_psi^2 * l^2 * E ~ E => S_psi ~ 1/l
    % 这里我们暂时设为 1.0 (如果收敛仍不好，可改为 1/mat.l_scale)
    scl.S_psi = 1.0 / mat.l_scale; 

    % 3. 乘子缩放因子 S_lam
    scl.S_lam = C_ref;
    
    fprintf('[Scaling] S_phi = %.2e (Voltage Scaling)\n', scl.S_phi);
    fprintf('[Scaling] S_psi = %.2e (Gradient Scaling)\n', scl.S_psi);
end

function [fixed_dofs, U_fixed_vals] = apply_dirichlet_bc(mesh, dof, BC_List, scl)
    n_total = dof.n_total; is_fixed = false(n_total, 1); U_fixed_vals = zeros(n_total, 1);
    for i = 1:length(BC_List)
        bc = BC_List{i};
        nodes = get_nodes_by_tag(mesh, bc.tag);
        if isempty(nodes), continue; end
        node_coords = mesh.coords(nodes, :); 
        switch lower(bc.type)
            case 'u',   target_nodes = nodes; cols = [1, 2];
            case 'phi', target_nodes = nodes(mesh.is_corner(nodes)); cols = 3;
            case 'psi', target_nodes = nodes(mesh.is_corner(nodes)); cols = [4, 5, 6, 7];
        end
        vals_input = bc.value; if ~iscell(vals_input), vals_input = num2cell(vals_input); end
        for k = 1:length(cols)
            if k > length(vals_input), continue; end
            val_item = vals_input{k};
            if isnumeric(val_item) && isnan(val_item), continue; end
            if length(target_nodes) ~= length(nodes), current_coords = mesh.coords(target_nodes, :);
            else, current_coords = node_coords; end
            dof_indices = dof.map(target_nodes, cols(k));
            valid_mask = dof_indices > 0;
            dof_indices = dof_indices(valid_mask);
            current_coords = current_coords(valid_mask, :);
            is_fixed(dof_indices) = true;
             scale_val = 1.0;
             if cols(k) == 3, scale_val = 1.0 / scl.S_phi;
             elseif cols(k) >= 4, scale_val = 1.0 / scl.S_psi; end
            if isa(val_item, 'function_handle')
                U_fixed_vals(dof_indices) = val_item(current_coords(:,1), current_coords(:,2)) * scale_val;
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

function nodes = get_nodes_by_tag(mesh, target_tags)
    if ~isfield(mesh, 'edge_tags'), nodes = []; return; end
    mask = ismember(mesh.edge_tags, target_tags);
    edge_indices = find(mask);
    if isempty(edge_indices), nodes = []; return; end
    nodes = unique(mesh.edge_elems(edge_indices, :));
end

function res = extract_results_flexible(U, dof, mesh, conf, scl)
    n_nodes = mesh.n_nodes;
    res.u = nan(n_nodes, 1); res.v = nan(n_nodes, 1);
    res.phi = nan(n_nodes, 1); res.psi11 = nan(n_nodes, 1);
    res.psi22 = nan(n_nodes, 1); res.psi12 = nan(n_nodes, 1);
    mask_u = (dof.map(:, 1) > 0); 
    if any(mask_u), res.u(mask_u)=U(dof.map(mask_u,1)); res.v(mask_u)=U(dof.map(mask_u,2)); end
    mask_phi = (dof.map(:, 3) > 0); 
    if any(mask_phi)
        res.phi(mask_phi) = U(dof.map(mask_phi, 3))* scl.S_phi;
        res.psi11(mask_phi) = U(dof.map(mask_phi, 4))* scl.S_psi;
        res.psi22(mask_phi) = U(dof.map(mask_phi, 5))* scl.S_psi;
        res.psi12(mask_phi) = U(dof.map(mask_phi, 6))* scl.S_psi;
    end
end

function visualize_results(mesh, res)
    hFig = figure('Color','w', 'Name', 'FlexoFEM Results');
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
    title(t, 'Flexoelectric Simulation Results (\psi_{ij} formulation)', 'FontSize', 14);
    colormap(hFig, 'jet');
    nexttile;
    max_disp = max(abs(res.v));
    scale = 1.0; if max_disp > 1e-12, scale = 0.1 * mesh.x_max / max_disp; end
    deformed_coords = mesh.coords + scale * [res.u, res.v];
    patch('Faces', mesh.elems(:, 1:4), 'Vertices', deformed_coords, ...
          'FaceVertexCData', res.v, 'FaceColor', 'interp', 'EdgeColor', 'none');
    title(sprintf('Disp v (Scale: %.1f)', scale)); axis equal; axis tight; colorbar;
    nexttile; plot_field_patch(mesh, res.phi, 'Electric Potential \phi (V)');
    nexttile; plot_field_patch(mesh, res.psi11, 'Auxiliary Strain \psi_{11}');
    nexttile; plot_field_patch(mesh, res.psi12, 'Auxiliary Strain \psi_{12}');
end

function plot_field_patch(mesh, data, title_str)
    if all(isnan(data)), text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center'); axis off;
    else
        patch('Faces', mesh.elems(:, 1:4), 'Vertices', mesh.coords, ...
              'FaceVertexCData', data, 'FaceColor', 'interp', 'EdgeColor', 'none');
        title(title_str); axis equal; axis tight; colorbar;
    end
end

function plot_boundary_result(mesh, results, target_tags, var_name, sort_mode, varargin)
    if nargin < 5 || isempty(sort_mode), sort_mode = 'auto'; end
    if nargin < 6 || isempty(varargin), varargin = {'k-', 'LineWidth', 1.5}; end
    if ~isfield(mesh, 'edge_tags'), warning('Mesh 缺少 edge_tags'); return; end
    mask = ismember(mesh.edge_tags, target_tags);
    edge_idx = find(mask);
    if isempty(edge_idx), fprintf('[Plot] Warning: Tags %s 未找到边界单元。\n', mat2str(target_tags)); return; end
    current_edges = mesh.edge_elems(edge_idx, :);
    node_indices = unique(current_edges(:));
    xy = mesh.coords(node_indices, :);
    x_span = max(xy(:,1)) - min(xy(:,1)); y_span = max(xy(:,2)) - min(xy(:,2));
    switch lower(sort_mode)
        case 'x', base = 1; case 'y', base = 2;
        case {'auto','arc'}, if x_span >= y_span, base = 1; else, base = 2; end
    end
    [~, sort_order] = sort(xy(:, base));
    sorted_nodes = node_indices(sort_order); sorted_xy = xy(sort_order, :);
    if strcmpi(sort_mode, 'arc'), d_vec = diff(sorted_xy, 1, 1); plot_x = [0; cumsum(sqrt(sum(d_vec.^2, 2)))]; axis_label = 'Arc Length (m)';
    else, plot_x = sorted_xy(:, base); if base==1, axis_label='X (m)'; else, axis_label='Y (m)'; end, end
    switch lower(var_name)
        case 'u',      vals = results.u(sorted_nodes); y_lab = 'u (m)';
        case 'v',      vals = results.v(sorted_nodes); y_lab = 'v (m)';
        case 'phi',    vals = results.phi(sorted_nodes); y_lab = '\phi (V)';
        case 'u_mag',  vals = sqrt(results.u(sorted_nodes).^2 + results.v(sorted_nodes).^2); y_lab = '|U| (m)';
        case {'psi11'}, vals = results.psi11(sorted_nodes); y_lab = '\psi_{11}';
        case {'psi22'}, vals = results.psi22(sorted_nodes); y_lab = '\psi_{22}';
        case {'psi12'}, vals = results.psi12(sorted_nodes); y_lab = '\psi_{12}';
    end
    valid_mask = ~isnan(vals); if sum(valid_mask)<2, return; end
    x_clean = plot_x(valid_mask); y_clean = vals(valid_mask);
    figure('Color','w', 'Name', sprintf('Boundary %s', mat2str(target_tags)));
    plot(x_clean, y_clean, varargin{:}); hold on; grid on; box on;
    [max_v, max_idx] = max(y_clean); [min_v, min_idx] = min(y_clean);
    plot(x_clean(max_idx), max_v, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    plot(x_clean(min_idx), min_v, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    text(x_clean(max_idx), max_v, sprintf(' Max: %.4e', max_v), 'VerticalAlignment', 'bottom', 'Color', 'r');
    text(x_clean(min_idx), min_v, sprintf(' Min: %.4e', min_v), 'VerticalAlignment', 'top', 'Color', 'b');
    xlabel(axis_label); ylabel(y_lab); title(sprintf('Boundary Tags %s | %s', mat2str(target_tags), y_lab)); hold off;
end