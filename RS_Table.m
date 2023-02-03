function T = RS_Table(L, S, N, heading, title)
    
   table = cell(length(N_values), length(dt_values));
table(:) = {''};

% row and column headings
row_headings = {'N = 3', 'N = 7', 'N = 15', 'N = 31'};
col_headings = {'dt = 1/64', 'dt = 1/128', 'dt = 1/256', 'dt = 1/512', 'dt = 1/1024', 'dt = 1/2048', 'dt = 1/4096'};

% Fill in the table with 'unstable' for the corresponding elements
for i = 1:length(N_values)
    for j = 1:length(dt_values)
        %if abs(T_inner(i,j,4)) > 1
            table{i,j} = 'unstable';
        end
    end
end

% Create the table
T = array2table(table, 'RowNames', row_headings, 'VariableNames', col_headings);

% Display the table
disp(T);