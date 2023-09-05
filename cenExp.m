% DESCRIPTION:
% This script runs simulations for centralized methods over the Benchmark
% K=500 network from https://web.stanford.edu/~yyye/Col.html.
% The script truncates the network by choosing a random subset of sensors
% and anchors, and by setting a required radio range. Each network is
% tested over several noise realizations. Comparison graphs are generated
% for PRMSE and function values. Available centralized methods are:
% 1. EML -       YALMIP implementation with SeDuMi according to Problem
%                (12) in the paper https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6725647&tag=1
% 2. EML_CVX -   CVX implemtation of 1
% 3. EMLNR_CVX - YALMIP implementation with SeDuMi accorsing to Problem (6)
%                in the paper https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6725647&tag=1
% 4. AMFC -      fully centralized version of AMU
% 5. ESDPD -     distribued version of the convex relaxation fotmulated in
%                the paper https://stanford.edu/~boyd/papers/pdf/semidef_relax_sensor_net_loc.pdf
% 6. EMLD -      distributed YALMIP version of EML (currently not working)

%% Initialization

clear
clc
warning off

%-----Hyper-parameters-----%
n=2; % dimension of the space (must be set to 2)
m=5; % total number of anchors in the truncated network (must not exceed 10)
max_deg=100; % maximum degree of neighbors for non-anchors (set it high enough so it has no effect, since maximum degree is dictated by the radius)
sigma=0.1; % noise level of the measurements across all networks and their realizations
netsize=[8,16,32,64,128]; % total number of non-anchors in the truncated network (must not exceed 490)
Ye=1; % if Ye=1 then the networks are generated accoridng to Ye 2008, otherwise according to our method.
radius=[0.385,0.49,0.595,0.7;... % each row corrensponds to an entry in netsize. A network is generated for each pais of (netsize(i),radius(i,:))
    0.36,0.47,0.586,0.7;...
    0.235,0.323,0.411,0.5;...
    0.225,0.25,0.275,0.3;...
    0.145,0.18,0.21,0.25];
methods={'EML','AMFC','ESDPD'}; % setting the required methods
nr=50; % number of realizations for each network
AMFC_iter=1000; % maximum number of AMFC iterations
%--------------------------%

%% Set the output structure

dt=datetime;
data.date=datestr(dt); % date and time of commencing the script
data.m=m; % number of anchors
data.sigma=sigma; % noise level
data.number_of_realizations=nr; % number of realizations
for i=1:length(netsize)
    N=netsize(i);
    for r=1:length(radius(i,:))
        R=radius(i,r);
        strR=['R',strrep(num2str(R),'.','')]; % radius of the network
        data.(['net',num2str(N)]).(strR).funv_real=zeros(nr,1); % initailize the average function value in the real location with zeros
        for j=1:length(methods) % setting the required metrics (function values and PRMSE) as 'funv' and 'norm2_diff'    
            data.(['net',num2str(N)]).(strR).(methods{j}).funv=zeros(nr,1); % initialize the function values of each method and realization with zeros
            data.(['net',num2str(N)]).(strR).(methods{j}).norm2_diff=zeros(nr,1); % initialize the MSE of each method and realization with zeros 
            data.(['net',num2str(N)]).(strR).(methods{j}).sum_norm2_diff=0; % initailize the MSE with zeros
        end
    end
end

%% Run algorithms and generate output

for i=1:length(netsize)
    rng(i) % fixes the seed
    N=netsize(i);
    x0=-0.01+0.02*rand(2*N,1); % the same starting point (AMFC) for each network across all radii and realizations
    
    for r=1:length(radius(i,:))
        R=radius(i,r);
        strR=['R',strrep(num2str(R),'.','')]; % string with used radius
        
        if Ye % truncates the network by setting number of non-anchors and radius
            [PP,dd]=generateD_revised(m,N,R,max_deg); % minimum radius for connectivitiy: (8,0.385), (16,0.36), (32,0.235), (64,0.225), (128,0.145)
        else
            [PP,dd]=truncateBench500(m,N,max_deg); % truncates the netwwork into a connected one by setting the number of non-anchors and anchors
        end
        
        for l=1:nr
            net=create_realization(m,PP,dd,sigma,R); % creates a realization of the noise and the output is a network of the required format
            u0=zeros(size(net.Matrices.Q_tilde,1)*2+size(net.Matrices.A_tilde,1)*2,1); % setting u0 as the zero vector (AMFC)
            
            Xreal=net.Matrices.X_real; % real location, last m columns are anchors
            [edges_nn,edges_na]=create_edges(net);
            edges_na(:,1)=edges_na(:,1)+length(edges_nn);
            edges=[edges_nn;edges_na]; % edges with real distances
            
            % run the algorithms and calculate the following: norm2_diff, sum_norm2_diff and the function value of the realization
            for j=1:length(methods)
                methods_alg=str2func(['alg_',methods{j}]);
                if contains('AMFC',methods{j})
                    alg_output=methods_alg(net,x0,u0,AMFC_iter); % run AMFC
                else
                    alg_output=methods_alg(net); % run EML, EML_CVX, EMLNR_CVX and ESDPD
                end
                
                X_output=alg_output.location_estimation; % output location, last m columns are anchors
                data.(['net',num2str(N)]).(strR).(methods{j}).X{l}=X_output; % save the algorithm's output location
                
                data.(['net',num2str(N)]).(strR).(methods{j}).norm2_diff(l)=sum(sum((X_output(:,1:N)-Xreal(:,1:N)).^2)); % norm2_diff for MSE
                data.(['net',num2str(N)]).(strR).(methods{j}).sum_norm2_diff=data.(['net',num2str(N)]).(strR).(methods{j}).sum_norm2_diff+data.(['net',num2str(N)]).(strR).(methods{j}).norm2_diff(l); % sum_norm2_diff is sum over the realizations
                for e=1:size(edges,1) % function value
                    edge=edges(e,:);
                    data.(['net',num2str(N)]).(strR).(methods{j}).funv(l)=data.(['net',num2str(N)]).(strR).(methods{j}).funv(l)+(norm(X_output(:,edge(2))-X_output(:,edge(3)))-edge(4))^2;
                end
            end
            
            for e=1:size(edges,1) % calculate the real function value of the realization
                edge=edges(e,:);
                data.(['net',num2str(N)]).(strR).funv_real(l)=data.(['net',num2str(N)]).(strR).funv_real(l)+(norm(Xreal(:,edge(2))-Xreal(:,edge(3)))-edge(4))^2;
            end
            
            data.(['net',num2str(N)]).(strR).nets{l}=net; % save the network structure of the realization
            
        end
        
        for j=1:length(methods) % calculate PRMSE, average function value and function value std over the realizations
            sum_norm2_diff=data.(['net',num2str(N)]).(strR).(methods{j}).sum_norm2_diff;
            data.(['net',num2str(N)]).(strR).(methods{j}).PRMSE=sqrt(sum_norm2_diff/nr)/N; % PRMSE of each method for a specfiic network
            data.(['net',num2str(N)]).(strR).(methods{j}).funv_avg=sum(data.(['net',num2str(N)]).(strR).(methods{j}).funv)/nr; % average function value of each method for a specfiic network
            data.(['net',num2str(N)]).(strR).(methods{j}).funv_std=std(data.(['net',num2str(N)]).(strR).(methods{j}).funv); % std of the function values of each method for a specfiic network
        end
        
        % calculate the average and std of the real function value over the realizations
        data.(['net',num2str(N)]).(strR).funv_real_avg=sum(data.(['net',num2str(N)]).(strR).funv_real)/nr; % average function value in the real location for a specfiic network
        data.(['net',num2str(N)]).(strR).funv_real_std=std(data.(['net',num2str(N)]).(strR).funv_real); % std of the function value in the real location for a specfiic network
    end
end

mkdir output
save(['output\cenExp_',datestr(now,'yyyy_mm_dd_HH_MM'),'.mat'],'data') % save the data file with a unique filename

%% Plot PRMSE

colors=get(gca,'ColorOrder'); % set colors

for i=1:length(netsize) % plot a figure for each network size
    N=netsize(i);
    for j=1:length(methods) % each plot contains a grpah for each method
        for r=1:length(radius(i,:)) % R values are the x axis
            R=radius(i,r);
            strR=['R',strrep(num2str(R),'.','')];
            PRMSEplot(r)=data.(['net',num2str(N)]).(strR).(methods{j}).PRMSE; % values to plot for each method
        end
        semilogy(PRMSEplot,'color',colors(j,:))
        hold on
        xticks(1:length(radius(i,:))) % number of x axis ticks
        xticklabels(radius(i,:)) % x axis labels
        set(gca,'YGrid','on','XGrid','off') % y axis grid lines
        title(['$\frac{\textrm{PRMSE}}{N}$ for $K=\ $',num2str(N+m),', $N=\ $',num2str(N),', $m=\ $',num2str(m),', $\sigma=\ $',num2str(sigma),', $L=\ $',num2str(nr)],'Interpreter','latex','fontsize',14)
        legend(methods,'fontsize',10)
    end
    hold off
    if i<length(netsize) % do not open a new figure after the last graph
        figure()
    end
end

%% Plot function values with std

colors=get(gca,'ColorOrder'); % set colors

figure()
for i=1:length(netsize) % plot a figure for each network size
    N=netsize(i);
    for r=1:length(radius(i,:)) % R values are compose the x axis
        R=radius(i,r);
        strR=['R',strrep(num2str(R),'.','')];
        funvrealplot(r)=data.(['net',num2str(N)]).(strR).funv_real_avg; % avg of function value for real location
    end
    semilogy(funvrealplot,'color','black','LineWidth',2) % plot function value of real location
    hold on
    for j=1:length(methods) % each plot contains a grpah for each method
        for r=1:length(radius(i,:)) % R values are compose the x axis
            R=radius(i,r);
            strR=['R',strrep(num2str(R),'.','')];
            funvplot(r)=data.(['net',num2str(N)]).(strR).(methods{j}).funv_avg; % values to plot for each method
        end
        semilogy(funvplot,'color',colors(j,:))
        xticks(1:length(radius(i,:))) % number of x axis ticks
        xticklabels(radius(i,:)) % x axis labels
        set(gca,'YGrid','on','XGrid','off') % y axis grid lines
        title(['Final function value with std for $K=\ $',num2str(N+m),', $N=\ $',num2str(N),', $m=\ $',num2str(m),', $\sigma=\ $',num2str(sigma),', $L=\ $',num2str(nr)],'Interpreter','latex','fontsize',14)
    end
    for j=1:length(methods) % plot std
        for r=1:length(radius(i,:)) % R values are compose the x axis
            R=radius(i,r);
            strR=['R',strrep(num2str(R),'.','')];
            max_vec(r)=max(data.(['net',num2str(N)]).(strR).(methods{j}).funv);
            min_vec(r)=min(data.(['net',num2str(N)]).(strR).(methods{j}).funv);
        end
        patch([1:length(radius(i,:)),length(radius(i,:)):-1:1,1],[max_vec,fliplr(min_vec),max_vec(1)],colors(j,:),'FaceAlpha',0.1,'EdgeColor','none')
    end
    legend(['Real location',methods],'fontsize',10)
    hold off
    if i<length(netsize) % do not open a new figure after the last graph
        figure()
    end
end