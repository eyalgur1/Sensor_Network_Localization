%% Check Sequential Time for Benchmark Networks

fields=fieldnames(data.realizations{end});
fields=fields(2:end); % methods in the data structure
num_methods=length(fields);

Xreal=data.original_network.net.Matrices.X_real;
nr=50;
N=data.original_network.net.K-data.original_network.net.anchors;
T=table;

for j=1:num_methods
    atst=0;
    for l=1:nr
        if contains(fields{j},'EML')
            atst=atst+data.realizations{l}.(fields{j}).time/nr;
        elseif contains(fields{j},'AMFD')
            atst=atst+data.realizations{l}.(fields{j}).total_estimated_parallel_time/nr;
        elseif contains(fields{j},'AMCC_100AGl')
            atst=atst+(sum(data.realizations{l}.(fields{j}).estimated_parallel_time(1:100))+sum(data.realizations{l}.(fields{j}).time(101:end)))/nr;
        else
            atst=atst+data.realizations{l}.(fields{j}).total_time/nr;
        end
    end
    T=[T;table({fields{j}},atst)];
end

T.Properties.VariableNames={'method_rando3000','avg_tot_seq_time'};
T