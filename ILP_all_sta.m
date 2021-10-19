 function [print_solution, solution, Xname, xInt, fInt] = ILP_all_sta(...
     cost_limit, l, l_end, r, R, c)

tic
addpath('./coverage_condition/')
addpath('./link_condition/')


n = length(l);
m = length(r);

%% objective_function
[T, Yname, Xname] = objective_funcion(l, n, m);

% =========================================================================
%% ������� ��� ����������� ��������

%% Condition 7 - EQUALITY

% � ����� ����� ���� ��������� ������ ���� �����
% $e_i = \sum\limits_{j=1}^m x_{ij}$
[A_7, b_7] = point_is_include_sta(T, n, m);
% [A_7, b_7] = add_condition_7(T, n, m);


%% condition 7 hatch SUM(xi1) = 1 -- INEQUALITY

% ������ ������� ������ ���� ��������� ������ � ����� �����.
[A_7h, b_7h] = sta_must_be_placed_in_only_one_point(T, n, m);
% [A_7h, b_7h] = add_condition_7hatch(T, n, m);

%% condition 7�,�  -- INEQUALITY

% �������� �������� �� ��������� ������ �������� �������, ����������� � 
% ����� $a_i$, � ����� 0, ���� � ����� $a_i$  ��� �������
[A_7a, b_7a] = sta_coverage_is_no_more_than_coverage_radius(...
    T, r, 'plus', n, m);
[A_7b, b_7b] = sta_coverage_is_no_more_than_coverage_radius(...
    T, r, 'minus', n, m);
% [A_7a, b_7a] = add_condition_7ab(T, r, 'plus', n, m);
% [A_7b, b_7b] = add_condition_7ab(T, r, 'minus', n, m);

%% Condition 8�,�  -- INEQUALITY


% ����� ������� �������� ����� ������ ����� ������� $a_i$ � $a_k$, ��� 
% ����������� �������, �� ����� ��������� ���������� ����� ����� �������.
[A_8a, b_8a] = coverage_sum_between_sta(T, l, l_end, n, 'a');
[A_8b, b_8b] = coverage_sum_between_sta(T, l, l_end, n, 'b');

% [A_8a, b_8a] = add_condition_8ab(T, l, l_end, n, 'a');
% [A_8b, b_8b] = add_condition_8ab(T, l, l_end, n, 'b');

%% Condition 12 Y0,n+1 = 0; E0,n+1 = 1. -- EQUALITY

% ��� ��������������� ������� ����� ����� ���������,������: 
%     - $y^+_0$ = 0;
%     - $y^-_0$ = 0;
%     - $y^+_{n+1}$ = 0;
%     - $y^-_{n+1}$ = 0;
%     - $e_0$ = 1;
%     - $e_{n+1}$ = 1;
[A_12, b_12] = gateway_condition(T, n);
% [A_12, b_12] = add_condition_12(T, n);

% =========================================================================
%% ������� ��� ����������� ����� ����� ������� ����������� �������

%% Condition (9�, 9�)  -- INEQUALITY

% ������� ������ ���� ��������� � ����� ������ $a_i$ � $a_k$
[A_9a, b_9a] = sta_must_be_placed_to_link(T, n, m, 'a'); 
[A_9b, b_9b] = sta_must_be_placed_to_link(T, n, m, 'b'); 

% [A_9a, b_9a] = add_condition_9ab(T, n, m, 'a'); 
% [A_9b, b_9b] = add_condition_9ab(T, n, m, 'b'); 

%% Condition 20 ����� -- EQUALITY

% ����������, ����� ������� $s_j$ � ����� $a_i$ ���� ������� �  ����� 
% ��������, ������������� � ����� $a_k$, ������ �� $a_i$ ($ k>i $) ��� 
% � ������ ������ $s_{m + 1}$

[A_20, b_20] = sta_is_connected_with_right_sta(T, n, m); % N10a
% [A_20, b_20] = add_condition_20(T, n, m); % N10a

%% ��� �� 25 ����� ��� ������� ������ �� k -- EQUALITY

% ����������, ����� ����� ������� $s_q$ � ����� $a_k$ ������ ��� ������ 
% ���� $s_{m + 1}$ ���� ������� � �� �������� $s_j$ � ����� $a_i$ ($k>i$) 

[A_28, b_28] = right_sta_is_also_connected_with_sta(T, n, m); 
% [A_28, b_28] = add_condition_28(T, n, m);

%% ��� �� 25 ����� ��� ������� ����� �� k -- EQUALITY

% ����� ������� $s_j$ � ����� $a_i$ ������ ���� ������� � ����� ��������, 
% ������������� � ����� $a_k$ ����� �� ����� $a_i$ ($k<i$) ��� � ����� 
% ������ $s_{m + 1}$
[A_29, b_29] = sta_is_connected_with_left_sta(T, n, m);
% [A_29, b_29] = add_condition_29(T, n, m);


%% 25 18 �������  -- EQUALITY
% ����� ����� ������� $s_q$ � ����� $a_k$ ����� ��� ����� ���� $s_{m + 1}$
% ������ ���� ������� �� ��������  $s_j$ � ����� $a_i$ ($k<i$)
[A_25, b_25] = left_sta_is_also_connectd_with_sta(T, n, m); 
% [A_25, b_25] = add_condition_25(T, n, m);

%%
%  ���� ������� $s_j$ � $s_q$ �������, �� ������������� ������ ����� 
%  ����������� ������� ������ ���� �� ������ ���������� ����� ������� 
%  $a_i$ � $a_k$, ��� ����������� ������� $s_i$ � $s_q$

% 23 �����
% ��������

[A_23, b_23] = link_to_the_left_sta(T, l, l_end, R, n, m);

% [A_23, b_23] = add_condition_23(T, l, l_end, R, n, m);

% 24 �����
% ��������
[A_24, b_24] = link_to_the_right_sta(T, l, l_end, R, n, m);
% [A_24, b_24] = add_condition_24(T, l, l_end, R, n, m);


%% Condition 13 Cost limit
% ��������� �����������
[A_13, b_13] = cost_condition(T, n, m, c, cost_limit);

% =========================================================================
%% Matrix preparation
f = table2array(T);
% coverage is not integer
% gateway_coverage + placed_sta_coverage + gateway_coverage) * 2
% and next values is integer
% intcon = (1 + n + 1) * 2 + 1 : length(f);

intcon = 1: length(f);
%% CONSTRAINTS
% linear inequality constraints.
total_A = [A_7h; A_7a; A_7b; A_8a; A_8b; A_9a; A_9b; ...
     A_23; A_24; A_13; ];
A = table2array(total_A);

b = [b_7h; b_7a; b_7b; b_8a; b_8b; b_9a; b_9b; ...
     b_23; b_24; b_13;]; 

% % linear equality constraints.
total_Aeq = [A_7; A_12; A_20;  A_28; A_29; A_25];
Aeq = table2array(total_Aeq);
beq = [b_7; b_12;  b_20;  b_28; b_29; b_25];

% =========================================================================
%% BOUNDS

% bound constraints.
total_AVarName = total_A.Properties.VariableNames;
[~, row_Yname, ~] = intersect(total_AVarName,Yname,'stable');
% bound constraints.
% low bound 
lb = zeros(1, width(total_A));

% upper bound 
ub = ones(1, width(total_A));
ub(row_Yname(1: end-2)) = inf;


% =========================================================================
%% Solution
options = optimoptions(@intlinprog,'OutputFcn',@savem_feasible_solutions);
[x, fval] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub, options);

    function stop = savem_feasible_solutions(x,optimValues,state,varargin)
        %SAVEMILPSOLUTIONS save all the integer solution found by 
        %intlinprog.

        % No reason to stop the solver from this function.
        stop = false;

        switch (state)
          case 'init'
            % Create new (or reset) variables in the base workspace  
            xInt = [];
            fInt = [];
          case 'iter'
            % 'x' is not empty only when there is a new integer solution 
            % found during branch and bound phase.
            if ~isempty(x)
              % Save integer xInt and fval
                fval = optimValues.fval;
                xInt = [xInt, x(:)];
                fInt = [fInt, optimValues.fval];
            end
          case 'done'
            % Nothing to do here.
        end
    end

solution = array2table(x');
solution.Properties.VariableNames = total_AVarName;

t = toc

% =========================================================================
%% Print Solution
Cost_estimate = calculate_constraints(solution, Xname, ...
    A_13);
Placed_stations = get_placed_sta(solution, Xname, n, m);
% Placed_stations
print_solution = ['Placed stations = [', num2str(Placed_stations),']', ...
    ' ; Total coverage = ', num2str(-fval), ' ; Cost = ', ...
    num2str(Cost_estimate)];

end

function [table, Yname, Xname] = objective_funcion(place, n, m)
%% yi
index = sort([0 : n + 1, 0 : n + 1]);
Yname = strings(1, length(index));

for i = 1:2:length(index)
    Yname(i) = ['y', num2str(index(i)), 'plus'];
    Yname(i+1) = ['y', num2str(index(i+1)), 'minus'];
end

f =[zeros(1,2), -1*ones(1,length(place)*2), zeros(1,2)];

%% xij
i = 1 : n;
j = 1 : m;
Xname = strings(1, length(i)*length(j));
index = 1;
for i = 1 : n 
    for j = 1 : m
        Xname(index) = ['x', num2str(i), num2str(j)];
        index = index + 1;
    end
end
f = [f, 0*ones(1, n * m)]; 

%% ei
Ename = strings(1, length(0 : n + 1));
index = 1;
for i = 0 : n + 1
    Ename(index) = ['e', num2str(i)];
    index = index +  1;
end
f = [f, 0 * ones(1,n + 2)];

%% zijk
Zname = {};

for i = 1 : n 
    for j = 1 : n
        Z = ['z0_', ... 
             num2str(m+1), '_', num2str(i), ...
             '_', num2str(j)];
        Zname = [Zname, Z];
    end
end


for i = 1 : n
    for j = 1 : m
        for k = 0 : n + 1
            if i ~= k
                if k == 0 || k == n + 1 % gateway S_0 and S_{n+1}
                    Z = ['z', num2str(i), '_', ... 
                        num2str(j), '_', num2str(k), ...
                        '_', num2str(m+1)];
                    Zname = [Zname, Z];
                
                else                
                    for q = 1 : m  %Stations
                        if (j~= q)
                            Z = ['z', num2str(i), '_', ... 
                                num2str(j), '_', num2str(k), ...
                                '_', num2str(q)];
                            Zname = [Zname, Z];
                        end
                    end
                end
              
            end
        end
    end
end

for i = 1 : n 
    for j = 1 : n
        Z = ['z', num2str(n+1), '_', ... 
             num2str(m+1), '_', num2str(i), ...
             '_', num2str(j)];
        Zname = [Zname, Z];
    end
end

f = [f, 0*ones(1,length(Zname))];

%% Total objective function
VarName = [Yname, Xname, Ename, Zname];

table = array2table(f);

table.Properties.VariableNames = VarName;
table.Properties.RowNames = {'f'};
end

function [A_13, b_13] = cost_condition(table, n, m, c, cost_limit)
tableVarName = table.Properties.VariableNames;
RowNames = {'SUMxij is less than Cost Limit'};
var_x = cell(m, n);
A = zeros(1, width(table));
b_13 = cost_limit;

for j = 1 : m
    for i = 1 : n
        var_x{j,i} = ['x', num2str(i), num2str(j)];
    end
    [~, row_x, ~] = intersect(tableVarName, var_x(j,:));
    A(1, row_x) = c(j);
end

A_13 = array2table(A,'VariableNames', tableVarName);
A_13.Properties.RowNames = RowNames;
end

function Cost = calculate_constraints(solution, xname, ...
    cost_ineq)
    solution_name = solution.Properties.VariableNames;
    [~, row_x, ~] = intersect(solution_name, xname);
    
    solution_array = table2array(solution(1, row_x));
    cost_array = table2array(cost_ineq(1, row_x));
    
    Cost = sum(cost_array .* solution_array);
end

function [Placed] = get_placed_sta(solution, xname, n, m)
    Placed = ones(1, n)*inf;
    
    solution_name = solution.Properties.VariableNames;
    [~, row_x, ~] = intersect(solution_name, xname);
    p = zeros(1, n*m);
    s = zeros(1, n*m);
    solution_x = table2array(solution(1, row_x));
    index = 1;

    for i = 1 : n
        for j = 1 : m
            p(index) = i;
            s(index) = j;
            
            
            if int8(solution_x(index)) == 1
                Placed(i) = j;
            end
            index = index + 1;
        end
    end
end
