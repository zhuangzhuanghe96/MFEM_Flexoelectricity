function plot_mesh_check_v3(node_coords, elems, boundaries)
% PLOT_MESH_CHECK 网格与域Tag检查
% 功能:
%   1. 根据 Domain Tag (最后一列) 对单元进行彩色填充区分。
%   2. 标注域编号 (D1, D2...)。
%   3. 绘制边界和二次单元节点。

    figure('Name', 'Mesh Inspection (Quadratic Support)', 'Color', 'w', 'NumberTitle', 'off');
    ax = axes('NextPlot', 'add', 'Box', 'on');
    axis equal; %grid on;
    xlabel('X'); ylabel('Y');
    title('Domain Tag Inspection');
    nTri = 0; nQuad = 0; nTri2 = 0; nQuad2 = 0;
    % --- 1. 收集所有存在的 Tag ---
    all_tags = [];
    if isfield(elems, 'tri') && ~isempty(elems.tri), all_tags = [all_tags; elems.tri(:,end)]; end %#ok<*AGROW>
    if isfield(elems, 'tri2') && ~isempty(elems.tri2), all_tags = [all_tags; elems.tri2(:,end)]; end
    if isfield(elems, 'quad') && ~isempty(elems.quad), all_tags = [all_tags; elems.quad(:,end)]; end
    if isfield(elems, 'quad2') && ~isempty(elems.quad2), all_tags = [all_tags; elems.quad2(:,end)]; end
    
    unique_tags = unique(all_tags);
    if isempty(unique_tags)
        warning('未检测到任何单元或Tag，请检查 import 是否成功。');
        return;
    end
    
    % 生成颜色表 (为每个域分配不同颜色)
    nDomains = length(unique_tags);
    dom_colors = parula(nDomains); 
    legend_entries = [];
    legend_labels = {};
    
    fprintf('正在绘制 %d 个子域...\n', nDomains);

    % --- 2. 循环绘制每个域 ---
    for i = 1:nDomains
        dom_id = unique_tags(i);
        color = dom_colors(i, :);
        
        % 收集该域的所有中心点用于标注
        dom_centers = [];
        
        % (A) 绘制 Quad / Quad2
        % ------------------------------------------------
        % 提取属于当前 Domain 的单元
        q_idx = []; q2_idx = [];
        if isfield(elems, 'quad') && ~isempty(elems.quad)
            q_idx = find(elems.quad(:, end) == dom_id);
        end
        if isfield(elems, 'quad2') && ~isempty(elems.quad2)
            q2_idx = find(elems.quad2(:, end) == dom_id);
        end
        
        % 统一绘制 (取前4个角点)
        faces_to_plot = [];
        if ~isempty(q_idx),  faces_to_plot = [faces_to_plot; elems.quad(q_idx, 1:4)]; end
        if ~isempty(q2_idx), faces_to_plot = [faces_to_plot; elems.quad2(q2_idx, 1:4)]; end
        
        if ~isempty(faces_to_plot)
            h = patch('Faces', faces_to_plot, ...
                'Vertices', node_coords, ...
                'FaceColor', color, ...
                'EdgeColor', 'k', ...
                'EdgeAlpha', 0.2, ...
                'FaceAlpha', 0.6);
            
            % 计算中心点用于后续标注
            % 简单取第一个单元的中心作为域标签位置 (或者取平均)
            vtx_indices = faces_to_plot; 
            pts_x = node_coords(vtx_indices, 1);
            pts_y = node_coords(vtx_indices, 2);
            dom_centers = [dom_centers; mean(mean(reshape(pts_x, size(vtx_indices)))), mean(mean(reshape(pts_y, size(vtx_indices))))];
            
            if isempty(legend_entries) || i > length(legend_entries)
                legend_entries(end+1) = h; 
                legend_labels{end+1} = sprintf('Domain %d', dom_id);
            end
        end
        
        % (B) 绘制 Tri / Tri2
        % ------------------------------------------------
        t_idx = []; t2_idx = [];
        if isfield(elems, 'tri') && ~isempty(elems.tri)
            t_idx = find(elems.tri(:, end) == dom_id);
        end
        if isfield(elems, 'tri2') && ~isempty(elems.tri2)
            t2_idx = find(elems.tri2(:, end) == dom_id);
        end
        
        faces_tri = [];
        if ~isempty(t_idx),  faces_tri = [faces_tri; elems.tri(t_idx, 1:3)]; end
        if ~isempty(t2_idx), faces_tri = [faces_tri; elems.tri2(t2_idx, 1:3)]; end
        
        if ~isempty(faces_tri)
            h = patch('Faces', faces_tri, 'Vertices', node_coords, ...
                  'FaceColor', color, 'EdgeColor', 'k', 'FaceAlpha', 0.6);
              
            vtx_indices = faces_tri;
            pts_x = node_coords(vtx_indices, 1);
            pts_y = node_coords(vtx_indices, 2);
            dom_centers = [dom_centers; mean(mean(reshape(pts_x, size(vtx_indices)))), mean(mean(reshape(pts_y, size(vtx_indices))))];
            
            if length(legend_entries) < i
                legend_entries(end+1) = h;
                legend_labels{end+1} = sprintf('Domain %d', dom_id);
            end
        end
        
        % (C) 在域中心标注 "D1", "D2"
        if ~isempty(dom_centers)
            % 计算整个域的质心
            center_of_mass = mean(dom_centers, 1);
            text(center_of_mass(1), center_of_mass(2), sprintf('D%d', dom_id), ...
                'FontSize', 14, 'FontWeight', 'bold', 'Color', 'w', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', 'k', 'Margin', 1);
        end
    end    
     %% 3. 绘制边界 (Edges)
    % ---------------------------------------------------------
    if nargin > 2 && ~isempty(boundaries)
        % 自动识别列结构：最后的一列是 Tag
        num_cols = size(boundaries, 2);
        tags = boundaries(:, end); 
        unique_tags = unique(tags);
        colors = lines(length(unique_tags)); 
        
        % legend_entries = [];
        % legend_labels = {};
        
        for i = 1:length(unique_tags)
            tag = unique_tags(i);
            % 找到属于该 Tag 的所有边
            bnd_indices = find(tags == tag);
            subset = boundaries(bnd_indices, :);
            
            % 提取绘图坐标 (使用 NaN 分隔技巧加速绘图)
            X_plot = []; Y_plot = [];
            X_mid = [];  Y_mid = []; % 用于绘制二次中间节点
            
            for k = 1:length(bnd_indices)
                n1 = subset(k, 1); 
                n2 = subset(k, 2);
                
                X_plot = [X_plot, node_coords(n1,1), node_coords(n2,1), NaN]; %#ok<AGROW>
                Y_plot = [Y_plot, node_coords(n1,2), node_coords(n2,2), NaN]; %#ok<AGROW>
                
                % 如果是二次边界 (cols >= 4: n1, n2, mid, tag)
                if num_cols >= 4
                    nm = subset(k, 3);
                    X_mid = [X_mid, node_coords(nm,1)]; %#ok<AGROW>
                    Y_mid = [Y_mid, node_coords(nm,2)]; %#ok<AGROW>
                end
            end
            
            % 绘制连线 (端点)
            h = plot(X_plot, Y_plot, 'Color', colors(i,:), 'LineWidth', 2);
            legend_entries(end+1) = h; %#ok<AGROW>
            legend_labels{end+1} = sprintf('Bnd Tag %d', tag); %#ok<AGROW>
            
            % % 绘制中间节点 (如果存在)
            % if ~isempty(X_mid)
            %     plot(X_mid, Y_mid, '.', 'Color', colors(i,:)*0.8, 'MarkerSize', 8, 'HandleVisibility', 'off');
            % end
            
            % 添加 Tag 文字标签 (取中间的一条边)
            mid_k = floor(length(bnd_indices)/2) + 1;
            sample_edge = subset(mid_k, :);
            % 计算文字位置
            txt_x = mean(node_coords(sample_edge(1:end-1), 1));
            txt_y = mean(node_coords(sample_edge(1:end-1), 2));
            
            text(txt_x, txt_y, sprintf('B%d', tag), ...
                 'Color', colors(i,:),'FontSize', 11, 'FontWeight', 'bold', ...
                 'BackgroundColor', 'w', 'Margin', 1);
        end
        
        if ~isempty(legend_entries)
           legend(legend_entries, legend_labels, 'Location', 'bestoutside');
        end
    end

    %% 4. 调试编号与标题
    % ---------------------------------------------------------
    nNode = size(node_coords, 1);
    
    % 如果节点数少，显示节点编号方便调试
    if nNode < 200
        text(node_coords(:,1), node_coords(:,2), string(1:nNode), ...
             'Color', [0.3 0.3 0.3], 'FontSize', 8, 'VerticalAlignment', 'bottom');
    end
    
    title_str = sprintf('Mesh: %d Nodes', nNode);

    if isfield(elems, 'quad ') && ~isempty(elems.quad ), nQuad = size(elems.quad, 1); end
    if isfield(elems, 'quad2') && ~isempty(elems.quad2), nQuad2 = size(elems.quad2, 1); end
    if isfield(elems, 'tri') && ~isempty(elems.tri), nTri = size(elems.tri, 1); end
    if isfield(elems, 'tri2') && ~isempty(elems.tri2), nTri2 = size(elems.tri2, 1); end

    if nQuad > 0, title_str = [title_str, sprintf(', %d Q1', nQuad)]; end
    if nQuad2 > 0, title_str = [title_str, sprintf(', %d Q2', nQuad2)]; end
    if nTri > 0, title_str = [title_str, sprintf(', %d T1', nTri)]; end
    if nTri2 > 0, title_str = [title_str, sprintf(', %d T2', nTri2)]; end
    
    title(title_str);
    hold off;
end