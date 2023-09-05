function [PP,dd]=truncateBench500(m,N,max_deg)
%% INFORMATION

% DESCRIPTION: truncates the Benchmark K=500 network by deleting sensors 
%              and edges. The code checks it the truncated network is
%              connected, if not it continues to generate networks.
% INPUTS:
% * m - total number of required anchor sensors (must not exceed 10)
% * N - total number of required non-anchor sensors (must not exceed 490)
% * max_deg - required maximum degree for non-anchors (if possible)

% OUTPUTS:
% * PP - real positions where last m columns are anchors
% * dd - sesnor adjacency matrix of the network with real distances

%% Generate a connected netwoek

load('net500bench.mat','net') % load the network

bins=2; % for checking connected components
while max(bins)>1 || min(sum(L))==1
    anchor_indexes=find(net.is_anchor);
    new_anchor_indeces=zeros(m,1);
    for i=1:m
        index=randi(length(anchor_indexes));
        new_anchor_indeces(i)=anchor_indexes(index);
        anchor_indexes(index)=[];
    end
    L=net.Matrices.Laplacian_net;
    diag_L=diag(L);
    L=L-diag(diag_L);
    Left_indeces=(1:490)';
    Network_indeces=new_anchor_indeces;
    for i=1:N-max_deg
        indexes_need_neigbors=find(sum(-L(Network_indeces,Network_indeces),2)<max_deg);
        index_set=find(((sum(-L(Left_indeces,Network_indeces(indexes_need_neigbors)),2)>=1)+(isempty(indexes_need_neigbors))).*(diag_L(Left_indeces)>=max_deg));
        index=randi(length(index_set));
        new_node=Left_indeces(index_set(index));
        Network_indeces=[new_node;Network_indeces];
        Left_indeces=setdiff(Left_indeces,new_node);
    end
    for i=1:max_deg
        indexes_need_neigbors=find(sum(-L(Network_indeces,Network_indeces),2)<max_deg);
        degree=max_deg;
        index_set=[];
        while isempty(index_set)
            index_set=find(((sum(-L(Left_indeces,Network_indeces(indexes_need_neigbors)),2)>=1)+isempty(indexes_need_neigbors)).*(sum(-L(Left_indeces,Network_indeces),2)>=degree));
            degree=degree-1;
        end
        index=randi(length(index_set));
        new_node=Left_indeces(index_set(index));
        Network_indeces=[new_node;Network_indeces];
        Left_indeces=setdiff(Left_indeces,new_node);
    end
    L=L(Network_indeces,Network_indeces);
    diag_L=(-sum(L,2));
    connections_num=diag_L;
    too_many_connect=find(connections_num>max_deg);
    s=[];
    t=[];
    for i=1:length(Network_indeces)
        new_indexes=i+1:N+m;
        find_connect=find(L(i,i+1:end));
        s=[s,repmat(i,1,length(find_connect))];
        t=[t,new_indexes(find_connect)];
    end
    G = graph(s,t);
    plot(G)
    bins = conncomp(G);
end

indeces=find(ismember(s,too_many_connect).* ismember(t,too_many_connect));
while ~isempty(indeces) &&  max(connections_num)>max_deg
    if connections_num(s(indeces(1)))>max_deg && connections_num(t(indeces(1)))>max_deg
        s_temp=s;
        s_temp(indeces(1))=[];
        t_temp=t;
        t_temp(indeces(1))=[];
        G_temp=graph(s_temp,t_temp);
        bins = conncomp(G_temp);
        if max(bins)==1
            connections_num(s(indeces(1)))=connections_num(s(indeces(1)))-1;
            connections_num(t(indeces(1)))=connections_num(t(indeces(1)))-1;
            s=s_temp;
            t=t_temp;
            indeces(indeces>indeces(1))=indeces(indeces>indeces(1))-1; % update indices after remove
        end
    end
    indeces(1)=[];
end
Edges=[Network_indeces(s),Network_indeces(t)];
Nodes=Network_indeces;
deleted_nodes=unique(setdiff(1:500,Nodes));

%% Generate the output matrices PP and dd

load('test10-500.mat','PP','dd')

PP(:,deleted_nodes)=[];

rel_edges=zeros(500);
for i=1:size(Edges,1)
    rel_edges(Edges(i,1),Edges(i,2))=1;
    rel_edges(Edges(i,2),Edges(i,1))=1;
end
dd=dd.*rel_edges;
dd(deleted_nodes,:)=[];
dd(:,deleted_nodes)=[];

end