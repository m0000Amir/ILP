function [A, b] = sta_must_be_placed_in_only_one_point(table, n, m)
%%
% INEQUALITY CONSTRAINTS

% Station must be placed in only on point

% $  \sum\limits_{i=1}^n x_{ij} \leq 1, \quad j = \overline{1,m}$

tableVarName = table.Properties.VariableNames;
RowNames = cell(1,m);
var_x = cell(m, n);
A = zeros(m, width(table));
b = zeros(m, 1);
for j = 1 : m
    for i = 1 : n
        var_x{j,i} = ['x', num2str(i), '_', num2str(j)];
    end
    [~, row_x, ~] = intersect(tableVarName, var_x(j,:));
    A(j, row_x) = 1;
    b(j, 1) = 1;
    RowNames{1,j} = ['SUMxi', num2str(j), '=1'];
end

A = array2table(A,'VariableNames', tableVarName);
A.Properties.RowNames = RowNames;
end
