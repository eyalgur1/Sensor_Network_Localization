function cluster_stc=MakeClusters(Network,method)
%% INFORMATION
% DESCRIPTION: generates clusters according to a specified method.

% INPUTS:
% * Network - the true network structure
% * method - {} - diagonal coloring is used.
%            {q,'r'} - q>=1 random clusters is used.
%            {q,'g'} - q>=1 geographical clusters is used.

%% INITIALIZATION

[Qtilde,~,Atilde,~,Btilde,~,~,~,~,~,...
    ~,~,~,K,~,num_non_anchor,~,~,~,~,...
    ~,~,~,~,~,~,~]=general_init(Network);

anchor_indices=num_non_anchor+1:K; % we assum anchors are always the last nodes

%% CREATING CLUSTERS
cluster_stc=struct;
if length(method)==0 % generates clusters according to diagonal coloring
    cluster_stc.method='diagonal_coloring';
    coloring=cell(K,1);
    Delta = max(diag(Network.Matrices.Laplacian_net));
    for i=1:K
        coloring{i}=[i;i];
    end
    for i=(Delta+2):K
        neigh_colors=[];
        for j=Network.node{i}.neighbors
            neigh_colors=[neigh_colors,coloring{j}(2)];
        end
        coloring{i}=[i; min(setdiff(1:Delta+1,neigh_colors))];
    end
    ordered_colors=[];
    for i=1:K
        ordered_colors=[ordered_colors,coloring{i}(2)];
    end
    q=max(ordered_colors);
    clusters=cell(q,1);
    for j=1:q
        cluster_j=[];
        for i=1:K
            if coloring{i}(2)==j
                cluster_j=[cluster_j,coloring{i}(1)];
            end
        end
        clusters{j}.nodes=setdiff(cluster_j,anchor_indices);
    end
elseif contains(method{2},'r') % generates q random clusters
    cluster_stc.method='random_clusters';
    q=method{1};
    non_anchor_num=K-numel(anchor_indices);
    q=min(q,non_anchor_num); % no empty clusters
    clusters=cell(q,1);
    perm_nodes=randperm(non_anchor_num);
    size_cluster=floor(non_anchor_num/q);
    for i=1:q-1
        clusters{i}.nodes=setdiff(perm_nodes(size_cluster*(i-1)+1:size_cluster*i),anchor_indices);
    end
    clusters{q}.nodes=setdiff(perm_nodes(size_cluster*(q-1)+1:end),anchor_indices);
else
    cluster_stc.method='geographical_clusters';
    K=Network.K;
    m=numel(anchor_indices);
    q=method{1};
 
    
    A=(Network.Matrices.Q_tilde'*Network.Matrices.Q_tilde);
    A=-(A-diag(diag(A)));
    G = graph(A);
    H = distances(G);
    
    B=(Network.Matrices.QD_tilde'*Network.Matrices.QD_tilde);
    B=-(B-diag(diag(B)));
    B=sqrt(B);
    G2 = graph(B);
    D = distances(G2);
    
    wanted_cluster_size=floor((K-m)/q);
    clusters=cell(q,1);
    Nodes=cell(q,1);
    cluster_size=zeros(q,1);
    clusters_not_balances=true;
    best_balanced=10;
    iter=0;
    MAX_ITER=100;
    while clusters_not_balances && iter<MAX_ITER
        permute=randperm(K-m);
        Cluster_Heads=permute(1:q);
        iter=iter+1;
        for inner_iter=1:2
            hops_from_cluster_heads=H(:,Cluster_Heads);
            dist_from_cluster_heads=D(:,Cluster_Heads);
            for i=1:q
                 Nodes{i}=[Cluster_Heads(i)];
            end
            for l=setdiff([1:K-m],Cluster_Heads)
                min_value=min(hops_from_cluster_heads(l,:));
                best_clusters=find(hops_from_cluster_heads(l,:)==min_value);
                [~,best_index]=min(dist_from_cluster_heads(l,best_clusters));
                i=best_clusters(best_index);
                Nodes{i}=[Nodes{i},l];
            end            
            for i=1:q               
                cluster_size(i)=length(Nodes{i});
            end
            for i=1:q
                not_in_cluster=setdiff(1:K-m,Nodes{i});
                hop_from_non_cluster=min(H(Nodes{i},not_in_cluster),[],2);
                dist_from_non_cluster=min(D(Nodes{i},not_in_cluster),[],2);
                max_hops=max(hop_from_non_cluster);
                indexes=find(hop_from_non_cluster==max_hops);
                [~,max_index]= max(dist_from_non_cluster(indexes));              
                Cluster_Heads(i)=Nodes{i}(indexes(max_index));
            end
            %Nodes
            [max_size,max_index]=max(cluster_size);
            [min_size,min_index]=min(cluster_size);
            
            clusters_not_balances=(max_size/wanted_cluster_size>1.2) | (min_size/wanted_cluster_size<0.8);
            best_balanced_in=(max_size-min_size)/wanted_cluster_size;
%             if (max_size-min_size)/wanted_cluster_size<best_balanced_in
%                 BestNodes_in=Nodes;
%                 best_balanced_in=(max_size-min_size)/wanted_cluster_size;
%                 updated_best_nodes=true;
%                 best_balanced_in
%             end
        end
        if best_balanced_in<best_balanced
           BestNodes=Nodes;
           best_balanced=best_balanced_in;
        end
        %Clusters_not_balances=false
    end
    for i=1:q
        clusters{i}.nodes=BestNodes{i};
    end
end

for i=1:q % setting data for each cluster
    nodes=clusters{i}.nodes; % relative indices of nodes in nonanchor_indices
    edges_with_nodes=find(sum(Atilde(:,nodes),2)==1);
    clusters{i}.anchor_edges=edges_with_nodes; % finding edges with anchors not in cluster
    edges_with_nodes=find(sum(Qtilde(:,nodes)~=0,2)); % finding edges with nonanchors not in cluster
    clusters{i}.nonanchor_edges=edges_with_nodes;
    if length(clusters{i}.anchor_edges)>1
        connected_anchors=find(sum(Btilde(clusters{i}.anchor_edges,:))>0);
    else
        connected_anchors=find(Btilde(clusters{i}.anchor_edges,:)>0);
    end
    clusters{i}.anchor_neigbhors=connected_anchors;
    relavent_nonanchors=setdiff(1:num_non_anchor,nodes);
    connected_nonanchors=(sum(abs(Qtilde(clusters{i}.nonanchor_edges,relavent_nonanchors)))>0);
    clusters{i}.nonanchor_neigbhors=relavent_nonanchors(connected_nonanchors);
end
cluster_stc.clusters=clusters;
cluster_stc.number_of_clusters=q;

end