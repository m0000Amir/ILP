function [A, b] = gateway_condition(table, n)
%%
% EQUALITY CONSTRAINTS

% Gateway is already placed:
%     - $e_0$ = 1;
%     - $e_{n+1}$ = 1;

% Gateway coverage is equal to 0:
%     - $y^+_0$ = 0;
%     - $y^-_0$ = 0;
%     - $y^+_{n+1}$ = 0;
%     - $y^-_{n+1}$ = 0;

tableVarName = table.Properties.VariableNames;
A = [];
RowNames = {};

mas = zeros(1, width(table));
var_y = ['y', num2str(0), 'plus'];
[~, row_y, ~] = intersect(tableVarName, var_y);
mas(1, row_y) = 1;
A = [A; mas];

mas = zeros(1, width(table));
var_y = ['y', num2str(0), 'minus'];
[~, row_y, ~] = intersect(tableVarName, var_y);
mas(1, row_y) = 1;
A = [A; mas];

mas = zeros(1, width(table));
var_y = ['y', num2str(n + 1), 'plus'];
[~, row_y, ~] = intersect(tableVarName, var_y);
mas(1, row_y) = 1;
A = [A; mas];

mas = zeros(1, width(table));
var_y = ['y', num2str(n + 1), 'minus'];
[~, row_y, ~] = intersect(tableVarName, var_y);
mas(1, row_y) = 1;
A = [A; mas];


mas = zeros(1, width(table));
var_e = ['e', num2str(0)];
[~, row_e, ~] = intersect(tableVarName, var_e);
mas(1, row_e) = 1;
A = [A; mas];

mas = zeros(1, width(table));
var_e = ['e', num2str(n + 1)];
[~, row_e, ~] = intersect(tableVarName, var_e);
mas(1, row_e) = 1;
A = [A; mas];

A = array2table(A,'VariableNames', tableVarName);
A.Properties.RowNames = RowNames;
b = [zeros(4,1); ones(2,1)];
end
