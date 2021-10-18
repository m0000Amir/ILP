function [A, b] = link_to_the_left_sta(table, place, place_end, ...
    R, n, m)
%%
% INEQUALITY CONSTRAINTS

% 
% Zijkq(Rjq 0 (ai - ak) >= 0 k = i - 1, ..., 0.
% 

tableVarName = table.Properties.VariableNames;
RowNames = {};
A = [];
place = [0, place, place_end];

for i = 1 : n
    for j = 1 : m
        for k = i - 1 : -1 : 0
            if k == 0
                var_z = ['z', num2str(i), '_', num2str(j), '_', ...
                num2str(k), '_', num2str(m+1)];

                RowNames = [RowNames, ['-z', num2str(i), num2str(j), ...
                num2str(k),num2str(m+1) '(R', num2str(j), '_',...
                num2str(m+1), '-(a', num2str(i), ...
                '-a', num2str(k), '))']];

                [~, row_z, ~] = intersect(tableVarName, var_z);

                mas = zeros(1, width(table));
                mas(1, row_z) = -(R(j,m+1) - (place(i+1) - place(k+1)));
                A = [A; mas];
                
            else
                for q = 1 : m
                    if j ~= q
                        var_z = ['z', num2str(i), '_', num2str(j), '_', ...
                        num2str(k), '_', num2str(q)];
                    
                        RowNames = [RowNames, ['-z', num2str(i), num2str(j), ...
                        num2str(k),num2str(q) '(R', num2str(j), '_',...
                        num2str(q), '-(a', num2str(i), ...
                        '-a', num2str(k), '))']];
                    
                        [~, row_z, ~] = intersect(tableVarName, var_z);
                    
                        mas = zeros(1, width(table));
                        mas(1, row_z) = -(R(j,q) - (place(i+1) - place(k+1)));
                        A = [A; mas];
                        
                    end
                end
            end
        end
    end
end

i = n + 1;
j = m + 1;

for k = i - 1 : -1 : 0
    if k == 0
        var_z = ['z', num2str(i), '_', num2str(j), '_', ...
        num2str(k), '_', num2str(m+1)];

        RowNames = [RowNames, ['-z', num2str(i), num2str(j), ...
        num2str(k),num2str(m+1) '(R', num2str(j), '_',...
        num2str(m+1), '-(a', num2str(i), ...
        '-a', num2str(k), '))']];

        [~, row_z, ~] = intersect(tableVarName, var_z);

        mas = zeros(1, width(table));
        mas(1, row_z) = -(R(j,m+1) - (place(i+1) - place(k+1)));
        A = [A; mas];

    else
        for q = 1 : m
            if j ~= q
                var_z = ['z', num2str(i), '_', num2str(j), '_', ...
                num2str(k), '_', num2str(q)];

                RowNames = [RowNames, ['-z', num2str(i), num2str(j), ...
                num2str(k),num2str(q) '(R', num2str(j), '_',...
                num2str(q), '-(a', num2str(i), ...
                '-a', num2str(k), '))']];

                [~, row_z, ~] = intersect(tableVarName, var_z);

                mas = zeros(1, width(table));
                mas(1, row_z) = -(R(j,q) - (place(i+1) - place(k+1)));
                A = [A; mas];

            end
        end
    end
end

A = array2table(A,'VariableNames', tableVarName);
A.Properties.RowNames = RowNames;
b = zeros(height(A),1);
end