function [Placed] = print_stations_placement(solution, n, m)
    Placed = ones(1, n)*inf;
    
    for i = 1 : n
        for j = 1 : m            
            var_x = ['x', num2str(i), '_', num2str(j)];
            if int8(solution{1, var_x}) == 1
                Placed(i) = j;
            end
        end
    end
end