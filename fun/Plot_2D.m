function Plot_2D(Network,plot_real_net,varargin)
%% INFORMATION
% DESCRIPTION: plots the 2D real sensor network and the output 2D sensor
%              network vs. the real 2D sensor network.

% INPUTS:
% * Network - the true network structure
% * plot_real_net - a logical argument. If 1 plots the real sensor network.
% * varargin - an optional argument for plotting the output network
%              obtained by an algorithm, or the clusters. Must be of the
%              form {varargin1, varargin2} where
%              - varargin1 - an output structure of an algorithm.
%              - varargin2 - a cluster structure.
%              note that each can be {}.

%% INITIALIZATION

if length(varargin)==1
    if length(varargin{1}{1})==1 % checks whether to plot the output network
        output = varargin{1}{1};
        X=output.location_estimation;
    end
    if length(varargin{1}{2})==1
        clusters=varargin{1}{2}.Clusters;
        clustering_method=varargin{1}{2}.method;
    end
end

[~,~,Atilde,~,~,~,~,~,Xreal,~,...
    IMnet,~,~,~,~,~,~,~,~,~,...
    ~,~,~,~,~,~,~]=general_init(Network);

%% PLOTTING THE REAL NETWORK

if plot_real_net
    origin_r=[];
    dest_r=[];
    
    for i=1:size(IMnet,1) % generates the real edges
        k=1;
        for j=1:size(IMnet,2)
            if IMnet(i,j)==1
                k=j;
            end
            if IMnet(i,j)==-1
                origin_r=[origin_r,[Xreal(1,k);Xreal(1,j)]];
                dest_r=[dest_r,[Xreal(2,k);Xreal(2,j)]];
            end
        end
    end
    
    figure(1) % plot
    hold on   
    axis equal
    axis([-0.5 0.5 -0.5 0.5])
    plot(origin_r,dest_r,'Color',[0.85,0.85,0.85],'HandleVisibility','off')
    markers = {'+','o','*','x','s','d','^','v','.'};
    marker_sizes={50,15,50,50,25,25,15,15,100};
    marker_face_colors={'none','filled','none','none','filled','filled','filled','filled','none'};
    if ~isempty(varargin)
        if ~isempty(varargin{1}{2})
            Labels=cell(length(clusters)+1,1);
            cluster_num=0;
            newDefaultColors = lines(length(clusters));
            for i=1:length(clusters)
                cluster_i_real=Xreal(:,clusters{i}.nodes);
                color=newDefaultColors(i,:);
                marker_ind=mod(i-1,8)+1;
                marker_type=markers{marker_ind};
                marker_size=marker_sizes{marker_ind};
                face_color=marker_face_colors{marker_ind};
                if strcmp(face_color,'filled')
                    face_color=color;
                end
                if strcmp(clustering_method,'geographical_clusters')
                    scatter(cluster_i_real(1,2:end),cluster_i_real(2,2:end),20,color,'filled')
                    scatter(cluster_i_real(1,1),cluster_i_real(2,1),50,'d','filled','MarkerFaceColor',color,'MarkerEdgeColor','black')
                else
                    if size(cluster_i_real,2)>0                        
                        cluster_num=cluster_num+1;
                        scatter(cluster_i_real(1,1:end),cluster_i_real(2,1:end),marker_size,color,'Marker',marker_type,'MarkerFaceColor',face_color,'LineWidth',1)
                        Labels{i}=['Cluster ',num2str(cluster_num)];                    
                    end
                end
            end
            scatter(Xreal(1,size(Atilde,2)+1:end),Xreal(2,size(Atilde,2)+1:end),100,'p','black','filled')
            Labels{cluster_num+1}='Anchors';
             if ~strcmp(clustering_method,'geographical_clusters') 
                legend(Labels{1:cluster_num+1},'Location','east outside','NumColumns',1,'FontSize',9);%,'position',[0.124,0.024,0.78,0.108])
             end
        else
            scatter(Xreal(1,size(Atilde,2)+1:end),Xreal(2,size(Atilde,2)+1:end),100,'p','black','filled')
            scatter(Xreal(1,1:size(Atilde,2)),Xreal(2,1:size(Atilde,2)),'blue','filled')
        end
    else
        scatter(Xreal(1,size(Atilde,2)+1:end),Xreal(2,size(Atilde,2)+1:end),100,'p','black','filled')
        scatter(Xreal(1,1:size(Atilde,2)),Xreal(2,1:size(Atilde,2)),'blue','filled')
    end  
end
hold off

%% PLOTTING THE OUTPUT VS. THE REAL

if ~isempty(varargin)
    if ~isempty(varargin{1}{1})
        origin_e = [Xreal(1,:); X(1,:)]; % generates the error lines between the output network and the real network
        dest_e=[Xreal(2,:);X(2,:)];
        
        figure(2) % plots the output network vs. the real network
        hold on
        axis equal
        plot(origin_e,dest_e,'black')
        scatter(Xreal(1,size(Atilde,2)+1:end),Xreal(2,size(Atilde,2)+1:end),'.r')
        scatter(Xreal(1,1:size(Atilde,2)),Xreal(2,1:size(Atilde,2)),'.b')
        scatter(X(1:2:end),X(2:2:end),'.g')
        hold off
    end
end

%% OUTPUT

if ~plot_real_net && isempty(varargin)
    disp('No output to plot!')
end