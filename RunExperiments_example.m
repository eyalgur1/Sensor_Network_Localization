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
% * AMU: AMFC, AMCC, AMFD, AMRC, AMGC
% * for PG, AMCC, AMFD, AMRCa nd AMGC, warm start initializations are availale. E.g., 'AMRC,SF' or 'AMCC,AG'

% All available methods:
% MethodsToRun={'ADMM','SF','PG','PG,SF','PG,,AG','AMFC',...
%     'AMCC','AMCC,SF','AMCC,AG','AMFD','AMFD,SF','AMFD,AG'...
%     'AMRC','AMRC,SF','AMRC,AG','AMGC','AMGC,SF','AMGC,AG'};

% Setting the input parameters:
% MethodsToRun={'ADMM','SF','PG','PG,SF','PG,,AG','AMFC',...
%     'AMCC','AMCC,SF','AMCC,AG','AMFD','AMFD,SF','AMFD,AG'...
%     'AMRC','AMRC,SF','AMRC,AG','AMGC','AMGC,SF','AMGC,AG'};

MethodsToRun={'AMFD'};


% iterations
ADMM_iter=5; SF_iter=5; PG_iter=5; AMFC_iter=2000; AMAG_iter=2000;
AMCC_iter=5; AMFD_iter=5; AMRC_iter=5; AMGC_iter=5;

warm_start_iter=[1,2]; L='l'; % warm start
q_RC=[2,10]; q_GC=[2,10]; % number of clusters
Ns=[1000]; m=20; % netwrok size
num_nets_in_folder=2;
S=[40];
eps=[0.004,0.003]; zeta=[0.2,0.05]; tau=[0.015,0.002]; % ADMM-H parameters
IsBench=1; IsFull=0; % benchmark networks
x0_bound=0.01; existing_x0=0; % sets the box [-x0_bound,x0+bound]^2 and DirName is to be changes accordingly.                                     

% LoadDir and SaveDir are defined according to the starting point, date, computer, etc.
LoadDir='C:/Users\eyal.gur.STAFF/OneDrive - Technion/06 - Papers/01_WSNL__Alternating_Minimization_Based_First_Order_Method_for_the_Wireless_Sensor_Network_Localization_Problem/00 - Datasets/Network Data and Algorithm Outputs/nets_data/K';
SaveDir='C:/Users/eyal.gur.STAFF/OneDrive - Technion/05 - Ph.D. Candidate/00 - Others/03 - ORSIS 2023 - 25.2.23/out_2023_16_04';
SaveDir=[SaveDir,num2str(x0_bound),'/K']; % Path for saving the outputs. Currently it is defined by x0_bound (random distribution).
warning off

%% RUNNING THE  METHODS

for i=1:length(Ns)                                                                      % for all network types
    
    if IsBench                                                                          % sets folder name for loading the network according to benchmark or not
        num_nets_in_folder=1;                                                           % there is only one folder name for each benchmark
        if IsFull
            LoadFolder=[LoadDir,num2str(Ns(i)),'_benchmark_full'];
        else
            LoadFolder=[LoadDir,num2str(Ns(i)),'_benchmark'];
        end
    else
        LoadFolder=[LoadDir,num2str(Ns(i)),'_m',num2str(m)];
    end
    
    for j=1:num_nets_in_folder
        display(['Network K=',num2str(Ns(i)),', m=',num2str(m),': ',num2str(j),'/',num2str(num_nets_in_folder)])
        if IsBench                                                                      % loads the single benchmark network or a random network, and sets the folder name for saving the outputs
            if IsFull
                Network=load([LoadFolder,'/net',num2str(Ns(i)),'benchfull.mat']);
                SaveFolder=[SaveDir,num2str(Ns(i)),'_benchmark_full'];
            else
                Network=load([LoadFolder,'/net',num2str(Ns(i)),'bench.mat']);
                SaveFolder=[SaveDir,num2str(Ns(i)),'_benchmark'];
            end
        else
            Network=load([LoadFolder,'/net',num2str(j),'/net',num2str(j),'.mat']);
            SaveFolder=[SaveDir,num2str(Ns(i)),'_m',num2str(m),'/net',num2str(j)];
        end
        if ~exist(SaveFolder, 'dir')
            mkdir(SaveFolder)
        end
        
        % Network sizes
        n=size(Network.net.Matrices.X_real,1);                                          % dimension
        K=Network.net.K;                                                                % total number of sensors
        m=Network.net.anchors;                                                          % number of anchors
        Ne=size(Network.net.Matrices.Q_tilde,1)+size(Network.net.Matrices.A_tilde,1);   % number of non-anchor edges
        
        % Starting point initilization, where other types are acceptable manually
        if K==1000
            if existing_x0
                sp=load('x01000.mat');
                x0=sp.x01000;
            else
                x0=-x0_bound+2*x0_bound*rand(n*K,1);
            end
        elseif K==500
            if existing_x0
                sp=load('x0500.mat');
                x0=sp.x0500;
            else
                x0=-x0_bound+2*x0_bound*rand(n*K,1);
            end
        end
        u0=zeros(Ne*n,1);
        
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
                    ADMM=ADMM_alg(Network.net,ADMM_iter,eps(i),zeta(i),tau(i),'out.mat',x0);
                    save([SaveFolder,'/ADMM.mat'],'ADMM')
                    
%SF-----------------------------------------------------------------------%
                case 'SF'
                    SF=SF_alg(Network.net,x0,SF_iter);
                    save([SaveFolder,'/SF.mat'],'SF')
                    
%PG-----------------------------------------------------------------------%
                case 'PG'
                    if WarmStart
                        for w=warm_start_iter                               
                            if contains('AG',InitMethod)
                                PG=PG_alg(Network.net,x0,PG_iter,InitMethod,w,L,u0);
                                save([SaveFolder,'/PG_',num2str(w),StrCell{2},L,'.mat'],'PG')
                            else                                            
                                PG=PG_alg(Network.net,x0,PG_iter,InitMethod,w,0,u0);
                                save([SaveFolder,'/PG_',num2str(w),StrCell{2},'.mat'],'PG')
                            end
                        end
                    else
                        PG=PG_alg(Network.net,x0,PG_iter);
                        save([SaveFolder,'/PG.mat'],'PG')
                    end
                    
%AMFC---------------------------------------------------------------------%
                case 'AMFC'
                    AMFC=alg_AMFC(Network.net,x0,u0,AMFC_iter);
                    save([SaveFolder,'/AMFC.mat'],'AMFC')

%AMAG---------------------------------------------------------------------%
                case 'AMAG'
                    for ss=1:length(S)
                        s = S(ss);
                        AMAG=alg_AMAG(Network.net,x0,u0,AMAG_iter,L,1000,s);
                        save([SaveFolder,'/AMAG_',num2str(s),'_Lip=',L,'.mat'],'AMAG')
                    end
                    
%AMCC---------------------------------------------------------------------%
                case 'AMCC'
                    if WarmStart
                        for w=warm_start_iter                               
                            if contains('AG',InitMethod)
                                AMCC=AMCC_alg(Network.net,x0,u0,AMCC_iter,InitMethod,w,L);
                                save([SaveFolder,'/AMCC_',num2str(w),StrCell{2},L,'.mat'],'AMCC')
                            else                                            
                                AMCC=AMCC_alg(Network.net,x0,u0,AMCC_iter,InitMethod,w,0);
                                save([SaveFolder,'/AMCC_',num2str(w),StrCell{2},'.mat'],'AMCC')
                            end
                        end
                    else
                        AMCC=AMCC_alg(Network.net,x0,u0,AMCC_iter);
                        save([SaveFolder,'/AMCC.mat'],'AMCC')
                    end
                    
%AMFD---------------------------------------------------------------------%
                case 'AMFD'
                    if WarmStart
                        for w=warm_start_iter                               
                            if contains('AG',InitMethod)
                                AMFD=alg_AMFD(Network.net,x0,u0,AMFD_iter,InitMethod,w,L);
                                save([SaveFolder,'/AMFD_',num2str(w),StrCell{2},L,'.mat'],'AMFD')
                            else                                            
                                AMFD=alg_AMFD(Network.net,x0,u0,AMFD_iter,InitMethod,w,0);
                                save([SaveFolder,'/AMFD_',num2str(w),StrCell{2},'.mat'],'AMFD')
                            end
                        end
                    else
                        AMFD=alg_AMFD(Network.net,x0,u0,AMFD_iter);
                        save([SaveFolder,'/AMFD.mat'],'AMFD')
                    end
                    
%AMRC---------------------------------------------------------------------%
                case 'AMRC'
                    
                    for q=1:length(q_RC)                                        
                        if WarmStart                                            
                            for w=warm_start_iter                               
                                if contains('AG',InitMethod)                 
                                    AMRC=AMRC_alg(Network.net,x0,u0,AMRC_iter,q_RC(q),InitMethod,w,L);
                                    save([SaveFolder,'/AMRC_',num2str(q_RC(q)),'C_',num2str(w),StrCell{2},L,'.mat'],'AMRC')
                                else                                            
                                    AMRC=AMRC_alg(Network.net,x0,u0,AMRC_iter,q_RC(q),InitMethod,w);
                                    save([SaveFolder,'/AMRC_',num2str(q_RC(q)),'C_',num2str(w),StrCell{2},'.mat'],'AMRC')
                                end
                            end
                        else                                                    
                            AMRC=AMRC_alg(Network.net,x0,u0,AMRC_iter,q_RC(q));
                            save([SaveFolder,'/AMRC_',num2str(q_RC(q)),'C','.mat'],'AMRC')
                        end
                    end
                    
%AMGC---------------------------------------------------------------------%
                case 'AMGC'                                                       
                    for q=1:length(q_GC)
                        if WarmStart
                            for w=warm_start_iter                              
                                if contains('AG',InitMethod)                
                                    AMGC=AMGC_alg(Network.net,x0,u0,AMGC_iter,q_GC(q),InitMethod,w,L);
                                    save([SaveFolder,'/AMGC_',num2str(q_GC(q)),'C_',num2str(w),StrCell{2},L,'.mat'],'AMGC')
                                else                                            
                                    AMGC=AMGC_alg(Network.net,x0,u0,AMGC_iter,q_GC(q),InitMethod,w);
                                    save([SaveFolder,'/AMGC_',num2str(q_GC(q)),'C_',num2str(w),StrCell{2},'.mat'],'AMGC')
                                end
                            end
                        else 
                            AMGC=AMGC_alg(Network.net,x0,u0,AMGC_iter,q_GC(q));
                            save([SaveFolder,'/AMGC_',num2str(q_GC(q)),'C','.mat'],'AMGC')
                        end
                    end
            end
        end
    end
end

disp('Run Completed!')