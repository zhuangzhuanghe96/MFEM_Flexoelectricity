function [node_coords, elems, boundaries] = import_mesh_mixed_v3(filename)
% IMPORT_MESH_MIXED_V2 读取 COMSOL mphtxt (支持二次单元)

    if ~isfile(filename), error('Mesh file not found!'); end
    lines = readlines(filename);

    %% 1. 解析节点 (Vertices)
    num_vtx_idx = find(contains(lines, 'number of mesh vertices'), 1);
    nNode = str2double(extractBefore(lines(num_vtx_idx), '#'));
    
    start_coord = find(contains(lines, 'Mesh vertex coordinates'), 1);
    node_coords = read_data_block(lines, start_coord, nNode, 2);

    %% 2. 解析单元 (Elements)
    elems.tri = [];   elems.tri2 = [];
    elems.quad = [];  elems.quad2 = [];
    boundaries = [];

    % 遍历所有类型块
    type_indices = find(contains(lines, '# Type #'));
    
    for i = 1:length(type_indices)
        start_idx = type_indices(i);
        header_line = lines(start_idx + 2); % 读取类型名称行
        
        % --- 处理边界 (Edges) ---
        if contains(header_line, 'edg')
            is_quadratic = contains(header_line, 'edg2');
            nodes_per_elem = ifthen(is_quadratic, 3, 2);
            [conn, tags] = parse_element_block(lines, start_idx, nodes_per_elem);
            % 存储格式: [节点列..., Tag]
            boundaries = [boundaries; conn, tags+1];%边界tag+1%与COMSOL界面显示对齐

        % --- 处理三角形 (Triangles) ---
        elseif contains(header_line, 'tri')
            if contains(header_line, 'tri2')
                [conn, tags] = parse_element_block(lines, start_idx, 6);
                elems.tri2 = [elems.tri2; conn, tags];
            else
                [conn, tags] = parse_element_block(lines, start_idx, 3);
                elems.tri = [elems.tri; conn, tags];
            end

        % --- 处理四边形 (Quadrilaterals) ---
        elseif contains(header_line, 'quad')
            if contains(header_line, 'quad2')
                % COMSOL quad2 默认 9 节点
                % raw_quad2 = parse_element_block(lines, start_idx, 9);
                [conn, tags] = parse_element_block(lines, start_idx, 9);
                % 二次四边形 (Lagrange): 
                % COMSOL 原序通常为张量积顺序: 
                % 1(0,0), 2(1,0), 3(0,1), 4(1,1), 5(0.5,0), 6(0,0.5), 7(0.5,0.5), 8(1,0.5), 9(0.5,1)
                %
                % 通用 FEM 顺序 (CCW 逆时针):
                % 角点: 1(BL)->2(BR)->3(TR)->4(TL)
                % 中点: 5(B)->6(R)->7(T)->8(L)
                % 中心: 9
                %
                % 对应的列重排映射:
                % 原列: 1 2 3 4 5 6 7 8 9
                % 目标: 1 2 4 3 5 8 9 6 7
                conn = conn(:, [1, 2, 4, 3, 5, 8, 9, 6, 7]);
                elems.quad2 = [elems.quad2; conn, tags];
                % elems.quad2 = raw_quad2; 
            else
                [conn, tags] = parse_element_block(lines, start_idx, 4);
                % 修正 0-based 并调整 CCW 顺序: 1,2,4,3
                conn = conn(:, [1, 2, 4, 3]);
                elems.quad = [elems.quad; conn, tags];
            end
        end
    end
%% 4. 灵活输出导入信息
    fprintf('导入完成: ');
    fprintf('%d Nodes,', nNode);
    if ~isempty(elems.tri),  fprintf(' %d Tri,', size(elems.tri,1)); end
    if ~isempty(elems.tri2), fprintf(' %d Tri2,', size(elems.tri2,1)); end
    if ~isempty(elems.quad), fprintf(' %d Quad,', size(elems.quad,1)); end
    if ~isempty(elems.quad2),fprintf(' %d Quad2,', size(elems.quad2,1)); end
    if ~isempty(boundaries), fprintf(' %d Edges\n', size(boundaries,1)); end
    % fprintf('导入完成: %d Nodes, %d Quad2, %d Edges.\n', ...
        % size(node_coords,1), size(elems.quad2,1), size(boundaries,1));
end

%% --- 内部辅助函数 ---

function [conn, tags] = parse_element_block(lines, start_idx, nodes_per_elem)
    % 提取子块
    subset = lines(start_idx:end);
    nElem_line = subset(find(contains(subset, 'number of elements'), 1));
    nElem = str2double(extractBefore(nElem_line, '#'));
    
    % 读取拓扑连接 (Elements)
    elem_start = find(contains(subset, '# Elements'), 1) + start_idx - 1;
    conn = read_data_block(lines, elem_start, nElem, nodes_per_elem);
    conn = conn + 1; % 转换为 1-based 索引
    
    % 读取几何实体标记 (Tags)
    tag_marker = find(contains(subset, 'geometric entity indices'), 1);
    if ~isempty(tag_marker)
        tag_start = tag_marker + start_idx - 1;
        tags = read_data_block(lines, tag_start, nElem, 1);
        % tags = tags + 1;%域单元无需+1转换
    else
        tags = ones(nElem, 1);
    end
end

function data = read_data_block(lines, start_line, count, cols)
    % 高效读取数据块，处理多行换行情况
    data = zeros(count, cols);
    curr = start_line + 1;
    c = 0;
    while c < count
        row_str = strtrim(lines(curr));
        if isempty(row_str) || startsWith(row_str, '#') || startsWith(row_str, '[')
            curr = curr + 1; continue;
        end
        vals = sscanf(row_str, '%f')';
        num_vals = length(vals);
        
        % 逻辑：填充直到填满当前行的列数，或者处理跨行数据
        % 这里简化处理：假设 COMSOL 的单元节点可能跨行排列
        temp_data = vals;
        while length(temp_data) < cols && c < count
             curr = curr + 1;
             next_vals = sscanf(strtrim(lines(curr)), '%f')';
             temp_data = [temp_data, next_vals];
        end
        
        % 分配数据
        num_to_fill = min(floor(length(temp_data)/cols), count - c);
        if num_to_fill > 0
            reshaped = reshape(temp_data(1:num_to_fill*cols), cols, num_to_fill)';
            data(c+1 : c+num_to_fill, :) = reshaped;
            c = c + num_to_fill;
        end
        curr = curr + 1;
    end
end

function val = ifthen(cond, v1, v2)
    if cond, val = v1; else, val = v2; end
end