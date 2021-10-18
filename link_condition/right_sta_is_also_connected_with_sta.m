function [A, b] = right_sta_is_also_connected_with_sta(table, n, m)
%%
% EQUALITY CONSTRAINTS

% Right station is also connected with station.

tableVarName = table.Properties.VariableNames;
A = [];
RowNames = {};

for k = 1 : n
    for q = 1 : m
        mas = zeros(1, width(table));
        var_z = {};
        var_x = ['x', num2str(k), num2str(q)];
        [~, row_x, ~] = intersect(tableVarName, var_x);
        mas(1, row_x) = 1;
        for i = k + 1 : n + 1
            if i ~= k
                if i == n + 1 
                    var_z = [var_z, ...
                        ['z', num2str(i), '_', ... 
                        num2str(m+1), '_', num2str(k), ...
                        '_', num2str(q)]];
                else
                    for j = 1 : m
                        if j ~= q
                            var_z = [var_z, ...
                            ['z', num2str(i), '_', ... 
                            num2str(j), '_', num2str(k), ...
                            '_', num2str(q)]];
                        end
                    end  
                end
            end 
        end

        [~, row_z, ~] = intersect(tableVarName, var_z);
        mas(1, row_z) = -1;
        if k == 0 || k == n + 1
            RowNames = [RowNames, ['x', num2str(k), num2str(m+1), ...
                '-SUMzijkq']];
        else
            RowNames = [RowNames, ['x', num2str(k), num2str(q), ...
                '-SUMzijkq']];
        end
        A = [A; mas];       
        
    end
        
end

A = array2table(A,'VariableNames', tableVarName);

b = zeros(height(A),1);
end

