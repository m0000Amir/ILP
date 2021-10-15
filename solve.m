%% LOAD JSON data
fname = 'data.json'; 
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

output_MILP = {};
%% GET SOLUTIONS
% for i = 1 : size(struct2table(val), 2)
% TODO: change for loop
for i = 4
    data = val.(append('x',int2str(i)));

    l = data.placement';
    l_end = data.gateway_placement(2);
    
    cost_limit = data.cost_limit;
    delay_limit = data.delay_limit;
    arrival_rate = data.arrival_rate;
    
    r = zeros(1, size(struct2table(val.(append('x',int2str(i))).sta), 2));
    R = zeros(1, size(struct2table(val.(append('x',int2str(i))).sta), 2));
    c = zeros(1, size(struct2table(val.(append('x',int2str(i))).sta), 2));
    mu = zeros(1, size(struct2table(val.(append('x',int2str(i))).sta), 2));
    
    for j = 1 : size(struct2table(val.(append('x',int2str(i))).sta), 2)
        r(j) = val.(append('x',int2str(i))).sta.(append('x',int2str(j))).r;
        R(j) = val.(append('x',int2str(i))).sta.(append('x',int2str(j))).R;
        c(j) = val.(append('x',int2str(i))).sta.(append('x',int2str(j))).c;
        mu(j) = val.(append('x',int2str(i))).sta.(append('x',int2str(j))).mu;
    end
    
    print_solution = ILP_coverage_problem_cost(cost_limit, ...
        arrival_rate, delay_limit, l, l_end, r, R, c, mu);
    
    output_MILP = [output_MILP; print_solution];
%     csvwrite('output.csv',print_solution)
          
end
% csvwrite('output.csv',{output_data{1}})

%  Convert cell to a table and use first row as variable names
T = cell2table(output_MILP)
 
% Write the table to a CSV file
writetable(T,'output_MILP.csv', 'Delimiter',';')

% fid = fopen('output.csv','w')
%     for i = 1 : length(output_data)
%         fprintf(fid,'%s, %s, %s, %s\n',c{i})
%     end
%   fclose(fid)