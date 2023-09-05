%% INFORMATION
% DESCRIPTION: runs a set of specified method on networks that are located
%              in directories with titles such as cd\WSNL\N=500\net1 and
%              saves the outputs in this same directory. The algorithms are
%              to be located in cd.

% INPUTS: manually entered according to the specifications in the
%         INITIALIZATION section.

%% INITIALIZATION
% Note: all entries in this section are to be entered manually accoring to
%       the properties of each sensor network.

% Specifications of the applicable methods:
% * ADMM - ADMM-Hybrid (convex/non-convex)
% * SF - ball relaxation method
% * PG - projected gradient method
% * EML - centralized SDP with convex relaxation
% * AMU: AMFC, AMCC, AMFD, AMRC, AMGC
% * for PG, AMCC, AMFD, AMRCa nd AMGC, warm start initializations are availale. E.g., 'AMRC,SF' or 'AMCC,AG'

% All available methods:
% MethodsToRun={'ADMM','SF','PG','PG,SF','PG,AG','AMFC',...
%     'AMCC','AMCC,SF','AMCC,AG','AMFD','AMFD,SF','AMFD,AG'...
%     'AMRC','AMRC,SF','AMRC,AG','AMGC','AMGC,SF','AMGC,AG'};

clear
clc

% Setting the input parameters:
MethodsToRun={'AMAG','AMFC','AMCC,AG'};

% iterations
ADMM_iter=3; SF_iter=3; PG_iter=3; AMFC_iter=2000; AMAG_iter=2000;
AMCC_iter=2000; AMFD_iter=3; AMRC_iter=3; AMGC_iter=3;

warm_start_iter=[100]; L='l'; % warm start (for PG, AMCC, AMFD, AMGC and AMRC)
q_RC=[2,10]; q_GC=[2,10]; % number of clusters (for AMRC and AMGC, respectively)
Ks=[1000]; m=[20]; % netwrok size and number of anchors
nr=10; sigma=0.061*0.07; %number of realizations and noise level
eps=[0.004]; zeta=[0.2]; tau=[0.015]; % ADMM-H parameters (for K=500: 0.004,0.2,0.015 and for K=1000: 0.003,0.05,0.002, respectively)
s_AMAG=99; r_AMAG=10;
IsBench=1; IsFull=0; % benchmark networks
x0_bound=0.01; existing_x0=1; % sets the box [-x0_bound,x0+bound]^2. Use exisitng x0 for the same starting point across all realizations

warning off

%% RUNNING THE METHODS

for i=1:length(Ks)
    %rng(i)
    N=Ks(i)-m(i);
    
    if IsBench                                                                          % sets folder name for loading the network according to benchmark or not
        Network_original=load(['net',num2str(Ks(i)),'bench.mat']);
    else
        Network_original=load(['net',num2str(Ks(i)),'.mat']);
    end
    
    % Network sizes
    n=size(Network_original.net.Matrices.X_real,1);                                          % dimension
    K=Network_original.net.K;                                                                % total number of sensors
    m=Network_original.net.anchors;                                                          % number of anchors
    Ne=size(Network_original.net.Matrices.Q_tilde,1)+size(Network_original.net.Matrices.A_tilde,1);   % number of non-anchor edges
    
    X_real=Network_original.net.Matrices.X_real;
    true_distances=Network_original.net.Matrices.true_distances;
    R=Network_original.net.GI.R;
    
    % Starting point initilization, where other types are acceptable manually
    if existing_x0
        sp=load(['x0',num2str(Ks(i)),'.mat']);
        x0=sp.x0;
    else
        x0=-x0_bound+2*x0_bound*rand(n*K,1);
    end
    
    u0=zeros(Ne*n,1);
    
    data.original_network=Network_original;
    for l=1:nr
        data.realizations{l}=struct;
    end
    
    for l=1:nr
        rng(1000+l)
        fprintf('K=%2d    m=%2d   Realization=%2d\n', Ks(i), m(i), l)
        
        Network=create_realization(m(i),X_real,true_distances,sigma,R);
        
        data.realizations{l}.net.dd_noise=Network.Matrices.dd_noise;
        data.realizations{l}.net.noised_distances=Network.Matrices.noised_distances;
        
        for k=1:numel(MethodsToRun)                                                     % executes all available methods to run
            
            method=MethodsToRun{k};                                                     % the k-th method
            StrCell=strsplit(method,',','CollapseDelimiters',false);                    % parsing the method's name
            
            if numel(StrCell)>=2 && ~isempty(StrCell{2})                                % warm start method if available
                WarmStart=true;
                switch StrCell{2}
                    case 'SF'
                        InitMethod='SF';
                    case 'AG'
                        InitMethod='AG';
                    otherwise
                        error('No such initialization method')
                end
            else
                WarmStart=false;
            end
            
            switch StrCell{1}
                
                %ADMM---------------------------------------------------------------------%
                case 'ADMM'
                    ADMM=alg_ADMM(Network,ADMM_iter,eps(i),zeta(i),tau(i),'out.mat',x0);
                    data.realizations{l}.ADMM=ADMM;
                    
                    %SF-----------------------------------------------------------------------%
                case 'SF'
                    SF=alg_SF(Network,x0,SF_iter);
                    data.realizations{l}.SF=SF;
                    
                    %SF-----------------------------------------------------------------------%
                case 'EML'
                    EML=alg_EML(Network);
                    data.realizations{l}.SF=EML;
                    
                    %PG-----------------------------------------------------------------------%
                case 'PG'
                    if WarmStart
                        for w=warm_start_iter
                            if contains('AG',InitMethod)
                                PG=alg_PG(Network,x0,PG_iter,InitMethod,w,L,u0);
                                data.realizations{l}.(['PG_',num2str(w),StrCell{2},L])=PG;
                            else
                                PG=alg_PG(Network,x0,PG_iter,InitMethod,w,0,u0);
                                data.realizations{l}.(['PG_',num2str(w),StrCell{2}])=PG;
                            end
                        end
                    else
                        PG=alg_PG(Network,x0,PG_iter);
                        data.realizations{l}.PG=PG;
                    end
                    
                    %AMFC---------------------------------------------------------------------%
                case 'AMFC'
                    AMFC=alg_AMFC(Network,x0,u0,AMFC_iter);
                    data.realizations{l}.AMFC=AMFC;
                    
                    %AMFC---------------------------------------------------------------------%
                case 'AMAG'
                    AMAG=alg_AMAG(Network,x0,u0,AMAG_iter,L,r_AMAG,s_AMAG);
                    data.realizations{l}.AMAG=AMAG;
                    
                    %AMCC---------------------------------------------------------------------%
                case 'AMCC'
                    if WarmStart
                        for w=warm_start_iter
                            if contains('AG',InitMethod)
                                AMCC=alg_AMCC(Network,x0,u0,AMCC_iter,InitMethod,w,L);
                                data.realizations{l}.(['AMCC_',num2str(w),StrCell{2},L])=AMCC;
                            else
                                AMCC=alg_AMCC(Network,x0,u0,AMCC_iter,InitMethod,w,0);
                                data.realizations{l}.(['AMCC_',num2str(w),StrCell{2}])=AMCC;
                            end
                        end
                    else
                        AMCC=alg_AMCC(Network,x0,u0,AMCC_iter);
                        data.realizations{l}.AMCC=AMCC;
                    end
                    
                    %AMFD---------------------------------------------------------------------%
                case 'AMFD'
                    if WarmStart
                        for w=warm_start_iter
                            if contains('AG',InitMethod)
                                AMFD=alg_AMFD(Network,x0,u0,AMFD_iter,InitMethod,w,L);
                                data.realizations{l}.(['AMFD_',num2str(w),StrCell{2},L])=AMFD;
                            else
                                AMFD=alg_AMFD(Network,x0,u0,AMFD_iter,InitMethod,w,0);
                                data.realizations{l}.(['AMFD_',num2str(w),StrCell{2}])=AMFD;
                            end
                        end
                    else
                        AMFD=alg_AMFD(Network,x0,u0,AMFD_iter);
                        data.realizations{l}.AMFD=AMFD;
                    end
                    
                    %AMRC---------------------------------------------------------------------%
                case 'AMRC'
                    
                    for q=1:length(q_RC)
                        if WarmStart
                            for w=warm_start_iter
                                if contains('AG',InitMethod)
                                    AMRC=alg_AMRC(Network,x0,u0,AMRC_iter,q_RC(q),InitMethod,w,L);
                                    data.realizations{l}.(['AMRC_',num2str(q_RC(q)),'C_',num2str(w),StrCell{2},L])=AMRC;
                                else
                                    AMRC=alg_AMRC(Network,x0,u0,AMRC_iter,q_RC(q),InitMethod,w);
                                    data.realizations{l}.(['AMRC_',num2str(q_RC(q)),'C_',num2str(w),StrCell{2}])=AMRC;
                                end
                            end
                        else
                            AMRC=alg_AMRC(Network,x0,u0,AMRC_iter,q_RC(q));
                            data.realizations{l}.(['AMRC_',num2str(q_RC(q)),'C'])=AMRC;
                        end
                    end
                    
                    %AMGC---------------------------------------------------------------------%
                case 'AMGC'
                    for q=1:length(q_GC)
                        if WarmStart
                            for w=warm_start_iter
                                if contains('AG',InitMethod)
                                    AMGC=alg_AMGC(Network,x0,u0,AMGC_iter,q_GC(q),InitMethod,w,L);
                                    data.realizations{l}.(['AMGC_',num2str(q_GC(q)),'C_',num2str(w),StrCell{2},L])=AMGC;
                                else
                                    AMGC=alg_AMGC(Network,x0,u0,AMGC_iter,q_GC(q),InitMethod,w);
                                    data.realizations{l}.(['AMGC_',num2str(q_GC(q)),'C_',num2str(w),StrCell{2}])=AMGC;
                                end
                            end
                        else
                            AMGC=alg_AMGC(Network,x0,u0,AMGC_iter,q_GC(q));
                            data.realizations{l}.(['AMGC_',num2str(q_GC(q)),'C'])=AMGC;
                        end
                    end
            end
        end
    end
end

% calculate sum_norm2_diff for PRMSE for each method, for iteration, over the realizations
fields=fieldnames(data.realizations{end});
fields=fields(2:end); % methods in the data structure
num_methods=length(fields); % number of methods in the data strcuture
for j=1:num_methods
    method_iter=length(data.realizations{end}.(fields{j}).norm2_diff); % number of iteartions used in the method
    sum_norm2_diff=zeros(method_iter,1);
    for i=1:method_iter
        for l=1:nr
            sum_norm2_diff(i)=sum_norm2_diff(i)+data.realizations{l}.(fields{j}).norm2_diff(i);
        end
    end
    data.sum_norm2_diff.(fields{j})=sum_norm2_diff;
end

mkdir output
save(['output\RunExpRea_',datestr(now,'yyyy_mm_dd_HH_MM'),'.mat'],'data')

disp('Run Completed!')