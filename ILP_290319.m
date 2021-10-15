function example()
tic
l = [20, 30, 40];
l_end = 50;
r = [20, 5];
R = [40, 20];
%% 5x7
%%
% l = [80, 120, 190, 250, 310, 370, 460];
% l_end = 500;
% R = [100, 120, 120, 120, 100];
% %
% r = [90, 60,  100, 90, 100];
%%
% r = [60  60 100  80 100];
%%
% r = [60  60  60  80 100];
%%
% r = [60  70  60  60 100 ];
%%
% r = [60  70  70  60 100];
%%
% r = [60 70 70 60 90];
%%
% r = [60 60 70 60 90];
%%
% r = [60 80 80 60 90];
%%
% r = [60 80 80 90 90];
%%
% r = [90 80 80 90 90 ];
%%

%% 5x9
% l = [40 110 135 200 260 310 350 390 455];
% l_end = 500;
% R = [100, 120, 100, 100, 120];
%%
% r = [50  50  40  50 100];
%%
% r = [40  60  40  80 100];
%%
% r = [40  60  40  80 100];
%%
% r = [60  60  40  80 100];
%%
% r = [100  60  50  80 100];
%%
% r = [60  60  70  80 100];
%%
% r = [80 90 70 80 50];
%%
% r = [80  80  70  80 100];
%%
% r = [80 80 80 80 70 ];
%%
% r = [80 80 90 80 70];
%%
% r = [100  80  90  80  70];
%%

%% 5x10
% l = [90 150 190 220 280 310 340 390 410 450];
% l_end = 500;
% R = [100 120 100 100 120];
%%
% r = [100  80  50  50  70];
%%
% r = [70 80 50 50 70];
%%
% r = [70 60 50 60 70];
%%
% r = [90 60 50 60 60];
%%
% r = [90 60 50 60 90];
%%
% r = [90 50 50 60 90];
%%
% r = [60 50 50 60 60];
%%
% r = [60 50 50 80 80];
%%
% r = [60  50  40  80 100];
%%
% r = [50  50  40  50 100 ];
%%

% %% 6x12
% l = [90 150 190 220 280 310 340 390 410 450 500 560];
% l_end = 600;
% R = [100 120 100 100 120 100];
%%
% r = [80  60  70 100  60  50];
%%
% r = [80  60  80 100  60  50];
%%
% r = [80 60 80 50 60 50];
%%
% r = [80 60 80 50 70 50];
%%
% r = [80 60 80 50 70 60];
%%
% r = [80 80 80 50 70 60];
%%
% r = [80 80 80 50 70 80];
%%
% r = [100  80  80  50  70  80];
%%
% r = [80 80 50 50 70 80];
%%
% r = [100  80  50  50  70  80];
%%

%% 6x13
% l = [90 150 190 220 280 310 340 390 410 450 500 560 600];
% l_end = 650;
% R = [100 120 100 100 120 100];
%%
% r = [90 50 80 70 70 90];
%%
% r = [90 80 80 70 70 90];
%%
% r = [90 80 80 90 90 90];
%%
% r = [90 60 80 90 90 90];
%%
% r = [100  60  80 100  90  90];
%%
% r = [100  60  80 100  90  80];
%%
% r = [100  60  70 100  90  80];
%%
% r = [80  60  70 100  90  80];
%%
% r = [80  60  70 100  60  80];
%%
% r = [80  60  70 100  60  50];
%%

n = length(l);
m = length(r);
%% objective function
[T, Yname] = objective_funcion(l, n, m);

%% Equality Condition 7
[A_4, b_4] = add_condition_7(T, n, m);

%% condition 7à,á
[A_7a, b_7a] = add_condition_7ab(T, r, 'plus', n, m);
[A_7b, b_7b] = add_condition_7ab(T, r, 'minus', n, m);

%% Condition 8à,á
[A_8a, b_8a] = add_condition_8(T, l, l_end, n, 'a');
[A_8b, b_8b] = add_condition_8(T, l, l_end, n, 'b');

%% Condition (9à, 9á)
[A_9a, b_9a] = add_condition_9(T, n, m, 'a');
[A_9b, b_9b] = add_condition_9(T, n, m, 'b');

%% Condition (10à, 10á)
[A_10a, b_10a] = add_condition_10_1(T, n, m, 'a');
[A_10b, b_10b] = add_condition_10_1(T, n, m, 'b');

%% Condition (11à, 11á)
[A_11a, b_11a] = add_condition_11(T, l, l_end, r, R, n, m, 'a');
[A_11b, b_11b] = add_condition_11(T, l, l_end, r, R, n, m, 'b');

%% ADD conditions
%% Equality Condition 12 SUM(xi1) = 1
[A_12, b_12] = add_condition_12(T, n, m);

%% Equality Condition 13 Y0,n+1 = 0; E0,n+1 = 1.
[A_13, b_13] = add_condition_13(T, n, m);

%% SUMZij1 = e1
[A_14, b_14] >= add_condition_14(T, n, m);

%% Matrix preparation
f = table2array(T);
intcon = 1 : width(T);
% linear inequality constraints.
total_A = [A_7a; A_7b; A_8a; A_8b; A_9a; A_9b; A_10a; A_10b; A_11a; A_11b; A_14];
A = table2array(total_A);
b = [b_7a; b_7b; b_8a; b_8b; b_9a; b_9b; b_10a; b_10b; b_11a; b_11b; b_14]; 

% linear equality constraints.
total_Aeq = [ A_4;  A_12; A_13; ];
Aeq = table2array(total_Aeq);
beq = [ b_4;  b_12; b_13;];

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
[x,fval, exitflag, output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
solution = array2table(x');
solution.Properties.VariableNames = total_AVarName;
t = toc
fval
end

function [table, Yname] = objective_funcion(place, n, m)
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
b_7 = [zeros(height(A_7),1)];
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

function [A_8, b_8] = add_condition_8(table, place, place_end, n, key)
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

function [A_9, b_9] = add_condition_9(table, n, m, key)
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
                        var_z = ['z', num2str(i), '_', num2str(j), '_', num2str(k)];
                        RowNames = [RowNames, ['z', num2str(i), '_', ...
                            num2str(j), '_', num2str(k), '-e', num2str(i)];];
                    case 'b'
                        var_e = ['e', num2str(k)];
                        var_z = ['z', num2str(i), '_', num2str(j), '_', num2str(k)];
                        RowNames = [RowNames, ['z', num2str(i), '_', ...
                            num2str(j), '_', num2str(k), '-e', num2str(k)];];
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

function [A_10, b_10] = add_condition_10(table, n, m, key)
tableVarName = table.Properties.VariableNames;
A = zeros(n, width(table));
RowNames = {};
for i = 1 : n
    var_z = {};
    var_e = ['e', num2str(i)];
    [~, row_e, ~] = intersect(tableVarName, var_e);
    A(i, row_e) = 1;
    switch key
        case 'a'
            for k = i + 1 : n + 1
                for j = 1 : m
                    var_z = [var_z, ...
                        ['z', num2str(i), '_', num2str(j), '_', num2str(k)]];
                end        
            end
        case 'b'
            for k = i - 1 : -1 : 0
                for j = 1 : m
                    var_z = [var_z, ...
                        ['z', num2str(i), '_', num2str(j), '_', num2str(k)]];
                end        
            end
    end
    RowNames = [RowNames, ['e', num2str(i),'-SUMz', num2str(i), ...
        num2str(j), num2str(k)]];
    [~, row_z, ~] = intersect(tableVarName, var_z);
    A(i, row_z) = -1;                 
end
A_10 = array2table(A,'VariableNames', tableVarName);
A_10.Properties.RowNames = RowNames;
b_10 = zeros(height(A_10),1);
end

function [A_10, b_10] = add_condition_10_1(table, n, m, key)
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
                        ['z', num2str(i), '_', num2str(j), '_', num2str(k)]];       
                end
            case 'b'
                for k = i - 1 : -1 : 0
                    var_z = [var_z, ...
                        ['z', num2str(i), '_', num2str(j), '_', num2str(k)]];       
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

function [A_11, b_11] = add_condition_11(table, place, place_end, r, R, n, m, key)
tableVarName = table.Properties.VariableNames;
RowNames = {};
A = [];
place = [0, place, place_end];

for i = 1 : n
    switch key
        case 'a'
            for k = i - 1 : -1 : 0
                for j = 1 : m
                    var_z = ['z', num2str(i), '_', num2str(j), '_', num2str(k)];
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
                    var_z = ['z', num2str(i), '_', num2str(j), '_', num2str(k)];
                    RowNames = [RowNames, ['-z', num2str(i), '_', num2str(j),...
                        '_', num2str(k), '(R', num2str(j), '-(a', ...
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

function [A_12, b_12] = add_condition_12(table, n, m)
tableVarName = table.Properties.VariableNames;
RowNames = cell(1,m);
var_x = cell(m, n);
A = zeros(m, width(table));
b_12 = zeros(m, 1);
for j = 1 : m
    for i = 1 : n
        var_x{j,i} = ['x', num2str(i), num2str(j)];
    end
    [~, row_x, ~] = intersect(tableVarName, var_x(j,:));
    A(j, row_x) = 1;
    b_12(j, 1) = 1;
    RowNames{1,j} = ['SUMxi', num2str(j), '=1'];
end

A_12 = array2table(A,'VariableNames', tableVarName);
A_12.Properties.RowNames = RowNames;
end

function [A_13, b_13] = add_condition_13(table, n, m)
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

A_13 = array2table(A,'VariableNames', tableVarName);
A_13.Properties.RowNames = RowNames;
b_13 = [zeros(4,1); ones(2,1)];
end

function [A_14, b_14] = add_condition_14(table, n, m)
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
                        ['z', num2str(i), '_', num2str(j), '_', num2str(k)]];
                [~, row_z, ~] = intersect(tableVarName, var_z);
                mas(1, row_z) = -1;
            end
        end
    end
    RowNames = [RowNames, ['e', num2str(k), ...
            '-SUMz', num2str(i), num2str(j), num2str(k)]];
    A = [A; mas];
end             

A_14 = array2table(A,'VariableNames', tableVarName);
A_14.Properties.RowNames = RowNames;
b_14 = zeros(height(A_14),1);
end
