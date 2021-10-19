function [A, b] = sta_must_be_placed_to_link(table, n, m, key)
%%
% INEQUALITY CONSTRAINTS

%Station must be placed in two points $a_i$ и $a_k$ to link them.

% key = 'a'

%  $$ z_{ijkq} \leq e_i , \quad i = \overline{1, n}; \quad j = 
% \overline{1, m}; \quad k = \overline{1,n}, k \neq i; \quad q = 
% \overline{1,m}, q \neq j; $$

% key = 'b'

% $$  z_{ijkq} \leq e_k , \quad k = \overline{1, n};  \quad j =
% \overline{1, m}; \quad i = \overline{1,n}, i \neq k; \quad q = \overline{1,m}, q \neq j.$$


tableVarName = table.Properties.VariableNames;
A = [];
RowNames = {};

%%
for i = 1 : n
    for j = 1 : m
        for k = 0 : n + 1
            if i ~= k
                switch key
                    case 'a'
                        index_e = i;
                    case 'b'
                        index_e = k;
                end
                if k == 0 || k == n + 1 ... gateway S_0 and S_{n+1}
                        
                    var_e = ['e', num2str(index_e)];
                    var_z = ['z', num2str(i), '_', ... 
                            num2str(j), '_', num2str(k), ...
                            '_', num2str(m+1)];
                    RowNames = [RowNames, ...
                        ['z', num2str(i), '_', ... 
                        num2str(j), '_', num2str(k), ...
                        '_', num2str(m+1), '-e', num2str(index_e)];];
                    
                    [~, row_e, ~] = intersect(tableVarName, var_e);
                    [~, row_z, ~] = intersect(tableVarName, var_z);
                    mas = zeros(1, width(table));
                    mas(1, row_e) = -1;
                    mas(1, row_z) = 1;
                    A = [A; mas];
                
                else                
                    for q = 1 : m ... Stations
                        if j~= q
                            var_e = ['e', num2str(index_e)];
                            var_z = ['z', num2str(i), '_', ... 
                                    num2str(j), '_', num2str(k), ...
                                    '_', num2str(q)];
                            RowNames = [RowNames, ...
                                ['z', num2str(i), '_', ... 
                                num2str(j), '_', num2str(k), ...
                                '_', num2str(q), '-e', num2str(index_e)];];
                            [~, row_e, ~] = intersect(tableVarName, var_e);
                            [~, row_z, ~] = intersect(tableVarName, var_z);
                            mas = zeros(1, width(table));
                            mas(1, row_e) = -1;
                            mas(1, row_z) = 1;
                            A = [A; mas];                            
                        end
                    end
                end
              
            end
        end
    end
end
%5

A = array2table(A,'VariableNames', tableVarName);
A.Properties.RowNames = RowNames;
b = zeros(height(A),1);
end
