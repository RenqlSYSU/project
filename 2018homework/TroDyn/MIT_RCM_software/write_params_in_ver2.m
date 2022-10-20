function write_params_in(p)

% This does not produce a fully-commented params_ver2.in file, but it is
% readable by the model. It would be better to produce a fully-descriptive
% file.

% After parameter skiplines(:,1), add skiplines(:,2) of whitespace
% The whitespace is an unfortunate byproduct of taking what used to be a
% human-readable params_in file and making it procedurally-generated. 
% See params.txt for the human-readable version.
skiplines = [5,2;
    20,2;
    28,2;
    30,2;
    34,2;
    40,2];



if length(p) ~= 42
    error('Inconsistent parameter input.');
end

b = fopen('params_ver2.in','w');
% This won't be read, but we need the whitespace.
fprintf(b,'generated by write_params_in\n\n\n');

s = 1; % where we are in skiplines.

intlist = [12,13,14,33,34]; %These parameters need to be listed as integers.
for i = 1:length(p)
    if isstr(p{i})
        fprintf(b,p{i});     
    else
        if i==20 % Albedo needs two decimal places.
            fprintf(b,'%-4.2f',p{i});
        elseif ~isempty(intersect(i,intlist)); %A couple of things need to be integers with no decimal.
            fprintf(b,'%-4.0f',p{i});
        else %Everything needs one decimal place.
            fprintf(b,'%-7.1f',p{i});
        end
    end
    fprintf(b,'\n');
    if i == skiplines(s,1)
        fprintf(b,repmat('\n',1,skiplines(s,2)));
        if s+1 <= length(skiplines) % if this condition isn't true, then we've reached the end of skiplines.
            s = s+1;
        end
    end 
    
end
fprintf(b,'\n');
fclose(b);