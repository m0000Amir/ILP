function [A, b] = sta_is_connected_with_left_sta(table, n, m)
%%
% EQUALITY CONSTRAINTS
% 
% Station must be connected with any left station or left gateway.
% 
tableVarName = table.Properties.VariableNames;
A = [];
RowNames = {};

for i = 1 : n
    for j = 1 : m
        mas = zeros(1, width(table));
        var_z = {};
        var_x = ['x', num2str(i), '_', num2str(j)];

        [~, row_x, ~] = intersect(tableVarName, var_x);
        mas(1, row_x) = 1;
        for k = 0 : i - 1
            if (k == 0)
                var_z = [var_z, ...
                            ['z', num2str(i), '_', ... 
                            num2str(j), '_', num2str(k), ...
                            '_', num2str(m+1)]];   
            else
                for q = 1 : m
                    if q ~= j
                        var_z = [var_z, ...
                        ['z', num2str(i), '_', ... 
                        num2str(j), '_', num2str(k), ...
                        '_', num2str(q)]];
                    end 
                end
                
            end
        end
        
        [~, row_z, ~] = intersect(tableVarName, var_z);
        mas(1, row_z) = -1;
        if i == 0 || i == n + 1
            RowNames = [RowNames, ['x', num2str(i), '_', num2str(m+1), ...
                '-SUMzijkq']];
        else
            RowNames = [RowNames, ['x', num2str(i), '_', num2str(j), ...
                '-SUMzijkq']];
        end
        A = [A; mas];       
        
    end
        
end

A = array2table(A,'VariableNames', tableVarName);
% A_29.Properties.RowNames = RowNames;
b = zeros(height(A),1);
end
