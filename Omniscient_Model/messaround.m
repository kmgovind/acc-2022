values = 1:1:5;
for i = values
    fprintf('%d', i);
    if i == values(end)
        fprintf('We''re done!');
    end
end