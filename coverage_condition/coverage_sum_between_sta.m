function [A, b] = coverage_sum_between_sta(...
    table, place, place_end, n, key)
%%
% INEQUALITY CONSTRAINTS

% Total coverage between placed station is no more than distance 
% between them.

% key == 'a'

% $ y_i^+ + y_k^- \leq \frac{l_k - l_i}{2} \cdot (e_i + e_k ) + 
% (2 - e_i - e_k ) \cdot L, \quad i = \overline{1,n},  \quad k = 
% \overline{i+1,n+1};$


% key == 'b'

% $   y_i^- + y_k^+  \leq \frac{l_i-l_k}{2} \cdot (e_i + e_k) + 
% (2 - e_i - e_k) \cdot L, \quad i = \overline{1,n}, \quad k = 
% \overline{i-1,0}$

tableVarName = table.Properties.VariableNames;
A = [];
RowNames = {};
if key == 'a'
    symbol1 = 'plus';
    symbol2 = 'minus';
elseif key == 'b'
    symbol1 = 'minus';
    symbol2 = 'plus';
end
place = [0, place, place_end];
for i = 1 : n
    yi = ['y', num2str(i), symbol1];
    [~, row_yi, ~] = intersect(tableVarName, yi);
    
    ei = ['e', num2str(i)];
    [~, row_ei, ~] = intersect(tableVarName, ei);
    
    switch key
        case 'a'        
            for k = i + 1 : n + 1
                yk = ['y', num2str(k), symbol2];
                [~, row_yk, ~] = intersect(tableVarName, yk);
                ek = ['e', num2str(k)];
                [~, row_ek, ~] = intersect(tableVarName, ek);
                
                mas = zeros(1, width(table));
                mas(1, row_yi) = 1;
                mas(1, row_ei) = -1*(0.5 * (place(k+1) - place(i+1)) ...
                    - place_end);
                
                mas(1, row_yk) = 1;
                mas(1, row_ek) = -1*(0.5 * (place(k+1) - place(i+1)) ...
                    - place_end);
                A = [A; mas];
                RowNames = [RowNames, ...
                    ['y', num2str(i), symbol1, '-y', num2str(k)]];
            end
        case 'b'
            for k = i - 1 : -1: 0
                yk = ['y', num2str(k), symbol2];
                [~, row_yk, ~] = intersect(tableVarName, yk);
                ek = ['e', num2str(k)];
                [~, row_ek, ~] = intersect(tableVarName, ek);
                
                mas = zeros(1, width(table));
                mas(1, row_yi) = 1;
                mas(1, row_ei) = -1*(0.5 * (place(i+1) - place(k+1)) ...
                    - place_end);
                
                mas(1, row_yk) = 1;
                mas(1, row_ek) = -1*(0.5 * (place(i+1) - place(k+1)) ...
                    - place_end);
                A = [A; mas];
                RowNames = [RowNames, ...
                    ['y', num2str(i), symbol1, '-y', num2str(k)]];
            end    
    end
end
A = array2table(A,'VariableNames', tableVarName);
A.Properties.RowNames = RowNames;
b = ones(height(A),1) * 2 * place_end;
end