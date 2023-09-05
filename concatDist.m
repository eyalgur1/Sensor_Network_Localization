data1=data;
%%
nr=50;
fields=fieldnames(data.realizations{1});

for l=1:nr
    for j=2:length(fields)
        data1.realizations{l}.(fields{j})=data.realizations{l}.(fields{j});
    end
end

for j=2:length(fields)
    data1.sum_norm2_diff.(fields{j})=data.sum_norm2_diff.(fields{j});
end
%%
data=data1;