%% Check Bias for Benchmark Networks

fields=fieldnames(data.realizations{end});
fields=fields(2:end); % methods in the data structure
num_methods=length(fields);

Xreal=data.original_network.net.Matrices.X_real;
nr=50;
N=data.original_network.net.K-data.original_network.net.anchors;
T=table;

for j=1:num_methods
    bias=zeros(size(Xreal(:,1:N)));
    for l=1:nr
        Xoutput=data.realizations{l}.(fields{j}).location_estimation;
        bias=bias+(Xoutput(:,1:N)-Xreal(:,1:N))/nr;
    end
    bias=reshape(bias,2*N,1);
    T=[T;table({fields{j}},norm(bias),norm(bias,'Inf'))];
end

T.Properties.VariableNames={'method_rando3000','norm2_bias','normInf_bias'};
T