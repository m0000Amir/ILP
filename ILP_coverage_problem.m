function ILP_coverage_problem()
tic
%% Input Data
% 2x3
l = [20, 30, 40];
l_end = 50;
r = [20, 5];
R = [40, 20];

n = length(l);
m = length(r);
%% objective function
[T, Xname, Yname] = objective_funcion(l, n, m);

%% ADD Inequality conditions

%% condition 7à,á
[A_7a, b_7a] = add_condition_7ab(T, r, 'plus', n, m);
[A_7b, b_7b] = add_condition_7ab(T, r, 'minus', n, m);

%% Condition 8à,á
[A_8a, b_8a] = add_condition_8ab(T, l, l_end, n, 'a');
[A_8b, b_8b] = add_condition_8ab(T, l, l_end, n, 'b');

%% Condition (9à, 9á)
[A_9a, b_9a] = add_condition_9ab(T, n, m, 'a');
[A_9b, b_9b] = add_condition_9ab(T, n, m, 'b');

%% Condition (10à, 10á)
[A_10a, b_10a] = add_condition_10ab(T, n, m, 'a');
[A_10b, b_10b] = add_condition_10ab(T, n, m, 'b');

%% Condition (10â)SUMZij1 >= e1
[A_10c, b_10c] = add_condition_10c(T, n, m);

%% Condition (11à, 11á)
[A_11a, b_11a] = add_condition_11(T, l, l_end, r, R, n, m, 'a');
[A_11b, b_11b] = add_condition_11(T, l, l_end, r, R, n, m, 'b');

%% ADD Equality conditions
%% Condition 7
[A_7, b_7] = add_condition_7(T, n, m);

%% Equality Condition 7 hatch SUM(xi1) = 1
[A_7h, b_7h] = add_condition_7hatch(T, n, m);

%% Equality Condition 13 Y0,n+1 = 0; E0,n+1 = 1.
[A_12, b_12] = add_condition_12(T, n);

%% Matrix preparation
f = table2array(T);
intcon = 1 : width(T);
% linear inequality constraints.
total_A = [A_7a; A_7b; A_8a; A_8b; A_9a; A_9b; ...
    A_10a; A_10b; A_11a; A_11b; A_10c];
A = table2array(total_A);

b = [b_7a; b_7b; b_8a; b_8b; b_9a; b_9b; ...
    b_10a; b_10b; b_11a; b_11b; b_10c]; 

% linear equality constraints.
total_Aeq = [A_7; A_7h; A_12];
Aeq = table2array(total_Aeq);
beq = [b_7; b_7h; b_12;];

% bound constraints.
total_AVarName = total_A.Properties.VariableNames;
[~, row_Yname, ~] = intersect(total_AVarName,Yname,'stable');
% bound constraints.
% low bound 
lb = zeros(1, width(total_A));

% upper bound 
ub = ones(1, width(total_A));
ub(row_Yname([1, 2])) = inf;
ub(row_Yname([end-1, end])) = inf;
ub(row_Yname(3 : end-2)) = inf;
%% Solution
[x,fval, ~, ~] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
solution = array2table(x');
solution.Properties.VariableNames = total_AVarName;
t = toc

Placed_stations = get_placed_sta(solution, Xname, n, m);
% Placed_stations
print_solution = ['Placed_stations = [', num2str(Placed_stations),']', ...
    ' ; Noncoverage = ', num2str(l_end + fval)];
fval
end

function [table, Xname, Yname] = objective_funcion(place, n, m)
%% yi
Yname = {'y0plus', 'y0minus'};
for i = 1 : n + 1
    Yname = [Yname, ['y', num2str(i), 'plus']];
    Yname = [Yname, ['y', num2str(i), 'minus']];
end
f =[zeros(1,2), -1*ones(1,length(place)*2), zeros(1,2)];

%% xij
Xname = {};
for i = 1 : n 
    for j = 1 : m
        Xname = [Xname, ['x', num2str(i), num2str(j)]];
    end
end
f = [f, 0*ones(1, n * m)]; 

%% ei
Ename = {};
for i = 0 : n + 1
    Ename = [Ename, ['e', num2str(i)]];
end
f = [f, 0 * ones(1,n + 2)];

%% zijk
index_place = [0, 1 : n, n + 1];
Zname = {};
for i = 1 : n
    for j = 1 : m
        for k = 0 : n + 1
            if i ~= k
                Zname = [
                    Zname, ['z', num2str(i), '_', ... 
                    num2str(j), '_', num2str(k)]
                    ];                
            end
        end
    end
end

f = [f, 0*ones(1,length(Zname))];

%% Total objective function
VarName = [Yname, Xname, Ename, Zname];

table = array2table(f);

table.Properties.VariableNames = VarName;
table.Properties.RowNames = {'f'};
end

function [A_7, b_7] = add_condition_7(table, n, m)
tableVarName = table.Properties.VariableNames;
A = zeros(n, width(table));
RowNames = {};

for i = 1 : n
    var_e = ['e', num2str(i)];
    [~, row_e, ~] = intersect(tableVarName, var_e);
    A(i, row_e) = 1;
    for j = 1 : m
        var_x = ['x', num2str(i), num2str(j)];
        [~, row_x, ~] = intersect(tableVarName, var_x);
        A(i, row_x) = -1;
    end
    RowNames = [RowNames, ['e', num2str(i)]];
end

A_7 = array2table(A,'VariableNames', tableVarName);
A_7.Properties.RowNames = RowNames;
b_7 = zeros(height(A_7),1);
end

function [A_7ab, b_7ab] = add_condition_7ab(table, r, symbol, n, m)
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
A_7ab = array2table(A,'VariableNames', tableVarName);
A_7ab.Properties.RowNames = RowNames;
b_7ab = zeros(height(A_7ab),1);
end

function [A_8, b_8] = add_condition_8ab(table, place, place_end, n, key)
tableVarName = table.Properties.VariableNames;
A = [];
b_2 = [];
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
A_8 = array2table(A,'VariableNames', tableVarName);
A_8.Properties.RowNames = RowNames;
b_8 = ones(height(A_8),1) * 2 * place_end;
end

function [A_9, b_9] = add_condition_9ab(table, n, m, key)
tableVarName = table.Properties.VariableNames;
A = [];
RowNames = {};
for i = 1 : n  
    for j = 1 : m
        for k = 0 : n + 1
            if i ~= k
                switch key
                    case 'a'
                        var_e = ['e', num2str(i)];
                        var_z = ['z', num2str(i), '_', num2str(j), '_', ...
                            num2str(k)];
                        RowNames = [RowNames, ...
                            ['z', num2str(i), '_', num2str(j), '_', ...
                            num2str(k), '-e', num2str(i)];];
                    case 'b'
                        var_e = ['e', num2str(k)];
                        var_z = ['z', num2str(i), '_', num2str(j), '_', ...
                            num2str(k)];
                        RowNames = [RowNames, ...
                            ['z', num2str(i), '_', num2str(j), '_', ...
                            num2str(k), '-e', num2str(k)];];
                end       
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

A_9 = array2table(A,'VariableNames', tableVarName);
A_9.Properties.RowNames = RowNames;
b_9 = zeros(height(A_9),1);
end

function [A_10, b_10] = add_condition_10ab(table, n, m, key)
tableVarName = table.Properties.VariableNames;
A = [];
RowNames = {};
for i = 1 : n
    for j = 1 : m
        mas = zeros(1, width(table));
        var_z = {};
        var_x = ['x', num2str(i), num2str(j)];
        [~, row_x, ~] = intersect(tableVarName, var_x);
        mas(1, row_x) = 1;
        switch key
            case 'a'
                for k = i + 1 : n + 1
                    var_z = [var_z, ...
                        ['z', num2str(i), '_', num2str(j), '_', ...
                        num2str(k)]];       
                end
            case 'b'
                for k = i - 1 : -1 : 0
                    var_z = [var_z, ...
                        ['z', num2str(i), '_', num2str(j), '_', ...
                        num2str(k)]];       
                end
        end
        [~, row_z, ~] = intersect(tableVarName, var_z);
        mas(1, row_z) = -1;
        RowNames = [RowNames, ['x', num2str(i), num2str(j), ...
            '-SUMz', num2str(i), num2str(j), num2str(k)]];
        A = [A; mas];
    end
end
A_10 = array2table(A,'VariableNames', tableVarName);
A_10.Properties.RowNames = RowNames;
b_10 = zeros(height(A_10),1);
end

function [A_11, b_11] = add_condition_11(table, place, place_end, ...
    r, R, n, m, key)
tableVarName = table.Properties.VariableNames;
RowNames = {};
A = [];
place = [0, place, place_end];

for i = 1 : n
    switch key
        case 'a'
            for k = i - 1 : -1 : 0
                for j = 1 : m
                    var_z = ['z', num2str(i), '_', num2str(j), '_', ...
                        num2str(k)];
                    RowNames = [RowNames, ['-z', num2str(i), num2str(j), ...
                        num2str(k), '(R', num2str(j), '-(a', ...
                        num2str(i), '-a', num2str(k), '))']];
                    [~, row_z, ~] = intersect(tableVarName, var_z);
                    
                    mas = zeros(1, width(table));
                    mas(1, row_z) = -(R(j) - (place(i+1) - place(k+1)));
                    A = [A; mas];
                end        
            end
        case 'b'
            for k = i + 1 : n + 1
                for j = 1 : m
                    var_z = ['z', num2str(i), '_', num2str(j), '_', ...
                        num2str(k)];
                    RowNames = [RowNames, ...
                        ['-z', num2str(i), '_', num2str(j),'_', ...
                        num2str(k), '(R', num2str(j), '-(a', ...
                        num2str(k), '-a', num2str(i), '))']];
                    [~, row_z, ~] = intersect(tableVarName, var_z);
                    
                    mas = zeros(1, width(table));
                    mas(1, row_z) = -(R(j) - (place(k+1) - place(i+1)));
                    A = [A; mas];
                end        
            end
    end
    
end
A_11 = array2table(A,'VariableNames', tableVarName);
A_11.Properties.RowNames = RowNames;
b_11 = zeros(height(A_11),1);
end

function [A_7h, b_7h] = add_condition_7hatch(table, n, m)
tableVarName = table.Properties.VariableNames;
RowNames = cell(1,m);
var_x = cell(m, n);
A = zeros(m, width(table));
b_7h = zeros(m, 1);
for j = 1 : m
    for i = 1 : n
        var_x{j,i} = ['x', num2str(i), num2str(j)];
    end
    [~, row_x, ~] = intersect(tableVarName, var_x(j,:));
    A(j, row_x) = 1;
    b_7h(j, 1) = 1;
    RowNames{1,j} = ['SUMxi', num2str(j), '=1'];
end

A_7h = array2table(A,'VariableNames', tableVarName);
A_7h.Properties.RowNames = RowNames;
end

function [A_12, b_12] = add_condition_12(table, n)
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

A_12 = array2table(A,'VariableNames', tableVarName);
A_12.Properties.RowNames = RowNames;
b_12 = [zeros(4,1); ones(2,1)];
end

function [A_10c, b_10c] = add_condition_10c(table, n, m)
tableVarName = table.Properties.VariableNames;
A = [];
RowNames = {};
for k = 0 : n + 1
    mas = zeros(1, width(table));
    var_z = {};
    var_e = ['e', num2str(k)];
    [~, row_e, ~] = intersect(tableVarName, var_e);
    mas(1, row_e) = 1;
    for i = 1 : n
        if i ~= k
            for j = 1 : m
                var_z = [var_z, ...
                        ['z', num2str(i), '_', num2str(j), '_', ...
                        num2str(k)]];
                [~, row_z, ~] = intersect(tableVarName, var_z);
                mas(1, row_z) = -1;
            end
        end
    end
    RowNames = [RowNames, ['e', num2str(k), ...
            '-SUMz', num2str(i), num2str(j), num2str(k)]];
    A = [A; mas];
end             

A_10c = array2table(A,'VariableNames', tableVarName);
A_10c.Properties.RowNames = RowNames;
b_10c = zeros(height(A_10c),1);
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
