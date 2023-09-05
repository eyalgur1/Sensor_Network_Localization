%% Initialization

clear
clc
n=2;
m=5; % total number of anchors in the network (must not exceed 10)
max_deg=100;
sigma=0.1;
netsize=[8]; % total number of non-anchors in the network (must not exceed 490)
R=[0.5:0.05:1]; % if Ye=1 then the networks are generated accoridng to Ye 2008, otherwise according to our method. %R=[0.385,0.385,0.235,0.225,0.145];
methods={'EML','AMFC'};%,'ESDPD'};
metrics={'norm2_diff'};
L=20;
AMFC_iter=100;
epsilon=10^(-5);

for r=1:length(R)
    Net=['R' num2str(R(r)*100)];
    data.(Net).m=m;
    data.(Net).sigma=sigma;
    data.(Net).max_deg=max_deg;
    data.(Net).L=L;
    for i=1:length(netsize)
        N=netsize(i);
        for j=1:length(methods)
            data.(Net).(['net',num2str(N)]).(methods{j}).sum_norm2_diff=0;
            data.(Net).(['net',num2str(N)]).(methods{j}).bias_real=zeros(n,N);
            data.(Net).(['net',num2str(N)]).(methods{j}).bias_eps(1:n*N,1)={zeros(n,N)};
            for k=1:length(metrics)
                data.(Net).(['net',num2str(N)]).(methods{j}).(metrics{k})=zeros(L,1);
            end
        end
        data.(Net).(['net',num2str(N)]).R=R(r);
    end
end

%% Run algorithms and generate output
maxN=netsize(end);
maxx0=-0.01+0.02*rand(n*maxN,1);
fields=fieldnames(data);
for k=1:length(fields)
    field=fields{k};
    for i=1:length(netsize)
        rng(i)
        N=netsize(i);
        x0=maxx0(1:n*N);
        R=data.(field).(['net',num2str(N)]).R;
        [PP,dd]=generateD_revised(m,N,R,max_deg); % truncates the network (Ye's code) with eps difference in one coordinate
        for l=1:L
            rng(1000+l)
            net=create_realization(m,PP,dd,sigma); % creates a realization of the noise and the output is a network of the required format
            Xreal=net.Matrices.X_real;
            
            data.(field).(['net',num2str(N)]).nets{l}=net; % save the network data of the realization
            for j=1:length(methods)
                methods_alg=str2func(['alg_',methods{j}]);
                if contains('AMFC',methods{j})
                    u0=zeros(size(net.Matrices.Q_tilde,1)*2+size(net.Matrices.A_tilde,1)*2,1);
                    alg_output=methods_alg(net,x0,u0,AMFC_iter); % run an algorithm
                else
                    alg_output=methods_alg(net);
                end
                Xoutput=alg_output.location_estimation;
                
                data.(field).(['net',num2str(N)]).(methods{j}).X{l}=Xoutput; % save the algorithm's output location
                data.(field).(['net',num2str(N)]).(methods{j}).norm2_diff(l)=sum(sum((Xoutput(:,1:N)-Xreal(:,1:N)).^2)); % norm2_diff
                data.(field).(['net',num2str(N)]).(methods{j}).bias_real=data.(field).(['net',num2str(N)]).(methods{j}).bias_real+(Xoutput(:,1:N)-Xreal(:,1:N))/L;
            end
        end
        
        data.(field).(['net',num2str(N)]).(methods{j}).PDE=zeros(n*N,n*N);
        for c=1:n*N
            [PP_eps,dd_eps]=generateD_revised_eps(m,N,R(i),max_deg,c,epsilon);   %%shimrit something may be wrong here do you change the net structure?
            
            for l=1:L
                dd_eps_noise=dd_eps+data.(field).(['net',num2str(N)]).nets{l}.Matrices.dd_noise;
                net_eps=create_realization(m,PP_eps,dd_eps_noise,0);
                data.(field).(['net',num2str(N)]).nets_eps{c}.net{l}=net_eps;
                
                Xreal_eps=net_eps.Matrices.X_real;
                
                [edges_nn,edges_na]=create_edges(data.(field).(['net',num2str(N)]).nets{l});
                edges_na(:,1)=edges_na(:,1)+length(edges_nn);
                edges=[edges_nn;edges_na];
                
                % run the algorithms and alculate the following: norm2_diff, sum_norm2_diff and the function value of the realization
                for j=1:length(methods)
                    
                    methods_alg=str2func(['alg_',methods{j}]);
                    if contains('AMFC',methods{j})
                        u0=zeros(size(net.Matrices.Q_tilde,1)*2+size(net.Matrices.A_tilde,1)*2,1);
                        alg_output_eps=methods_alg(net_eps,x0,u0,AMFC_iter);
                    else
                        alg_output_eps=methods_alg(net_eps);
                    end
                    Xoutput_eps=alg_output_eps.location_estimation;
                    
                    data.(field).(['net',num2str(N)]).(methods{j}).X_eps{c}.X{l}=Xoutput_eps;
                    data.(field).(['net',num2str(N)]).(methods{j}).sum_norm2_diff=data.(field).(['net',num2str(N)]).(methods{j}).sum_norm2_diff+data.(field).(['net',num2str(netsize(i))]).(methods{j}).norm2_diff(l); % sum_norm2_diff
                    data.(field).(['net',num2str(N)]).(methods{j}).bias_eps{c}=data.(field).(['net',num2str(N)]).(methods{j}).bias_eps{c}+(Xoutput_eps(:,1:N)-Xreal_eps(:,1:N))/L;
                end
            end
        end
        
        
        
        %%
        for j=1:length(methods) % calculate PRMSE, average function value and function value std over the realizations
            for c=1:n*N
                data.(field).(['net',num2str(N)]).(methods{j}).PDE(:,c)=reshape(data.(field).(['net',num2str(N)]).(methods{j}).bias_eps{c}(:,1:N)-data.(field).(['net',num2str(N)]).(methods{j}).bias_real(:,1:N),n*N,1)/epsilon;
            end
            % calculate sqrt(CRLB)/N (% independent of the realization, but needs the network in our format (therefore it is outside the loop at its end))
            bias_real=reshape(data.(field).(['net',num2str(N)]).(methods{j}).bias_real,n*N,1);
            J=data.(field).(['net',num2str(N)]).(methods{j}).PDE;
            %J=data.(['net',num2str(N)]).AMFC.PDE;
            data.(field).(['net',num2str(N)]).(methods{j}).CRLB=sqrt(CRLB_eps(net,sigma,J,bias_real))/N;
            
            sum_norm2_diff=data.(field).(['net',num2str(netsize(i))]).(methods{j}).sum_norm2_diff;
            data.(field).(['net',num2str(N)]).(methods{j}).PRMSE=sqrt(sum_norm2_diff/L)/N;
        end
    end
end