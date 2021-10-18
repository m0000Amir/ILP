function [A, b] = ... 
    sta_coverage_is_no_more_than_coverage_radius(table, r, symbol, n, m)
%%
% INEQUALITY CONSTRAINTS

% Coverage value of placed station is no more than station's coverage 
% radius.

% symbol == 'plus'
% $ y_i^+ \leq \sum\limits_{j=1}^m x_{ij} \cdot r_j, \quad i = \overline{1,n};$

% symbol == 'minus'
% $ y_i^- \leq \sum\limits_{j=1}^m x_{ij} \cdot r_j, \quad i = \overline{1,n}.$

tableVarName = table.Properties.VariableNames;
A = zeros(n, width(table));
RowNames = {};
for i = 1 : n
    var_y = ['y', num2str(i), symbol];
    [~, row_y, ~] = intersect(tableVarName, var_y);
    A(i, row_y) = 1;
    for j = 1 : m
        var_x = ['x', num2str(i), num2str(j)];
        [~, row_x, ~] = intersect(tableVarName,var_x);
    A(i, row_x) = -1 * r(j);
    end
    RowNames = [RowNames, ['y', num2str(i), symbol]];
end
A = array2table(A,'VariableNames', tableVarName);
A.Properties.RowNames = RowNames;
b = zeros(height(A),1);
end
