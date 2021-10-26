function [A, b] = point_is_include_sta(table, n, m)
%%
% EQUALITY CONSTRAINTS

% Variable e_i is eqal to 1, if any station is placed on it.

%   e_i =  \sum\limits_{j=1}^m x_{ij}, \quad i = \overline{1,n}.

tableVarName = table.Properties.VariableNames;
A = zeros(n, width(table));
RowNames = {};

for i = 1 : n
    var_e = ['e', num2str(i)];
    [~, row_e, ~] = intersect(tableVarName, var_e);
    A(i, row_e) = 1;
    for j = 1 : m
        var_x = ['x', num2str(i), '_', num2str(j)];
        [~, row_x, ~] = intersect(tableVarName, var_x);
        A(i, row_x) = -1;
    end
    RowNames = [RowNames, ['e', num2str(i)]];
end

A = array2table(A,'VariableNames', tableVarName);
A.Properties.RowNames = RowNames;
b = zeros(height(A),1);
end

