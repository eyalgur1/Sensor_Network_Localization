function out = alg_EMLD(net,max_iter,acc)

%--------2020--------%
% NOTE 5/8/20: It is not clear what is the output location of this code.
% Therefore, it is not working (currently, yields no output). Nevertheless,
% the run time of each iteration is ~30 seconds for a N=30 network, which
% is not practicle.
%--------------------%

%Distributed ML- ESDP, as seen in TSP.

%PP: 2xn matrix representing all point on 2D
%n : the number of total points
%m : the number of anchor points; the first m columns of PP
%r : the radio range
%nf: noisy factor
%degree: the limit on the number of edges connected to any sensor point, typically set it to 5-10
%KK: number of iterations
%SNR: a tentative to add noisy links
%acc: the accuracy of the optimization problem

%Output
%PP: 2xn matrix representing all point on 2D, the last m columns are anchors
%and many other things.

%%
%Generate distances:

%--------2020--------%
m=net.anchors;
K=net.K;
r=net.GI.R;
degree=100; % maximum degree of neighbors for non-anchors (set it high enough so it has no effect, since maximum degree is dictated by the radius)
Xreal=net.Matrices.X_real;
amatrix=reshape(net.Matrices.a,2,m);
PP=zeros(size(Xreal));
PP(:,1:m)=Xreal(:,K-m+1:end); % GenerateDistances takes PP with anchors as the first m columns
PP(:,m+1:end)=Xreal(:,1:K-m);
noised_distances=sparse(triu(net.Matrices.noised_distances)); % GenerateDistances takes dd (measured distances) as an upper trinagular matrix
%--------------------%

[PP, dtrue, dnoise, Adjacency] = GenerateDistances(PP,m,K,r,0,degree,noised_distances);
[i1,~]=find(Adjacency(1:end-m,1:end-m));
Edges = size(i1,1);

%%
%Here I solve the centralized problem. First I solved it as in MLESDP
%(so normal way). Then I solved it, as if I were solving the problem
%the ADMM will solve (with the extra equality constraint). In this way
%I have the correct costs, lambda, etcetera.

%SOLVE THE CENTRALIZED PROBLEM, FIRST WAY------------------------------
%[PP,dtrue,m,Xiter] = MLESDP_inside_m1(PP, dtrue, dnoise, Adjacency, K, m);
%----------------------------------------------------------------------

%SOLVE THE CENTRALIZED PROBLEM-----------------------------------------
%[Popt, Costopt, LambdaOpt, vectLam] = MLESDP_inside_v2(PP, dnoise, Adjacency, m,K);
%----------------------------------------------------------------------

%Here I plot the results, just to check everything is fine [Optional]
% for ii = 1:K-m
%     plot(Xiter(1,ii), Xiter(2,ii),'k.'); hold on
%     plot(Popt(end-2*(K-m)+2*ii-1,:), Popt(end-2*(K-m)+2*ii,:), 'rs')
%     axis([-1 1 -1 1])
% end

%%
%Now, I go for the distributed part.

%Some definitions
AdjacencyT = Adjacency + Adjacency'; dnoiseT = dnoise+dnoise';
Yiimat_TR = zeros(K-m,K-m); Yijmat_TR= zeros(Edges,K-m);
Dijmat_TR= zeros(Edges,K-m);dijmat_TR= zeros(Edges,K-m);
Xmat_TR = zeros(2*(K-m), K-m);
YY_TR = [Yiimat_TR; Yijmat_TR; Dijmat_TR; dijmat_TR; Xmat_TR];
Lambdamat = zeros(size(YY_TR,1),2*Edges);
YY_R = YY_TR;
options = sdpsettings('verbose',0);

ivectmat = find(AdjacencyT(1:end-m,1:end-m));
PostprocessingMat_it = zeros(2*m+size(YY_TR,1),max_iter*(K-m));
PostP = zeros(2*m+size(YY_TR,1),max_iter*(K-m));
w = zeros(max_iter); %Iterations

starta=tic;

%Iterations:
for kk=1:size(w,1)
    %tell me which iteration it is:
    kk
    temp1 = tic;
    
    Lambdamat_N = zeros(size(Lambdamat,1),size(Lambdamat,2));
    Pimat_TR = zeros(2*m,K-m);
    
    for ii=1:K-m
        
        %For each sensor:
        %1) construct/recover the local variables:
        ij=find(AdjacencyT(ii,1:end-m));
        ik=find(AdjacencyT(ii,end-m+1:end));
        mi = length(ik); ntoti = length(ij)+1;ni = ntoti+mi;
        loc_Adj = zeros(ni,ni); loc_dnoise = zeros(ni,ni);
        loc_i = sort([ii, ij, ik+K-m]); i_loc = find(loc_i == ii);
        loc_Adj(i_loc, 1:ni) = AdjacencyT(ii, loc_i);
        loc_Adj(1:ni, i_loc) = AdjacencyT(loc_i,ii);
        loc_dnoise(i_loc, 1:ni) = dnoiseT(ii, loc_i);
        loc_dnoise(1:ni, i_loc) = dnoiseT(loc_i,ii);
        zij = zeros(3*(2*ntoti-1),ntoti-1);
        
        vecX = [2*loc_i(1:ntoti)-1, 2*loc_i(1:ntoti)]; vecX = sort(vecX);
        vectori = [loc_i(1:ntoti), K-m+ij, K-m+Edges+ij, K-m+2*Edges+ij, ...
            K-m+3*Edges+vecX];
        vectori_R = [loc_i(1:ntoti), K-m+ii*ones(1,ntoti-1), ...
            K-m+Edges+ii*ones(1,ntoti-1), K-m+2*Edges+ii*ones(1,ntoti-1), ...
            K-m+3*Edges+vecX];
        iv = ij+(K-m)*(ii-1);
        [~,iilambda] = intersect(ivectmat, iv);
        for ee=1:ntoti-1
            ei = ij(ee);
            zij(:,ee) = (YY_TR(vectori,ii)+YY_R(vectori_R,ei))/2;
        end
        
        Yijold = YY_TR(vectori,ii);
        lambdaij = Lambdamat(vectori,iilambda);
        
        %2) Solve for y and p,--------------------------------------------
        %   rho is set to 0.3
        [Yij, pi, lambdaij] = solvemysdp(.3, zij, lambdaij, ...
            Yijold, PP(1:2,loc_i), triu(loc_dnoise), triu(loc_Adj), ni, mi, acc);
        %------------------------------------------------------------------
        
        %3) Save the result in the global variables
        Yiimat_TR(loc_i(1:ntoti),ii) = Yij(1:ntoti);
        Yijmat_TR(ij,ii) = Yij(ntoti+1:2*ntoti-1);
        Dijmat_TR(ij,ii) = Yij(2*ntoti:3*ntoti-2);
        dijmat_TR(ij,ii) = Yij(3*ntoti-1:4*ntoti-3);
        vecX = [2*loc_i(1:ntoti)-1, 2*loc_i(1:ntoti)]; vecX = sort(vecX);
        Xmat_TR(vecX ,ii) = Yij(4*ntoti-2:end);
        Lambdamat_N(vectori,iilambda) = lambdaij;
        if ik>0
            Pimat_TR([ik, ik+m],ii) = pi;
        end
        
    end
    
    %save, clean, postprocess
    Lambdamat = Lambdamat_N;
    clear Lambdamat_N
    
    %Big Y
    YY_TR = [Yiimat_TR; Yijmat_TR; Dijmat_TR; dijmat_TR; Xmat_TR];
    PostprocessingMat_it(:,1+(K-m)*(kk-1):(K-m)*(kk)) = [Pimat_TR; YY_TR];
    
    %Sending procedure --- with possibilities to add SNR (so
    %TRansmitted data is not equal to Received data. (This is
    %Beta-version, not sure is properly working.
    %YY_R = YY_TR+SNR*randn(size(YY_TR,1), size(YY_TR,2));
    
    SNR = 0; %just in case, since I am not sure it is working (not tested properly).
    YY_R = YY_TR+SNR*YY_TR.*(2*rand(size(YY_TR,1), size(YY_TR,2))-1)*.5;
    
    %Get me the time for iteration (so I estimate the total time)
    toc(temp1)
    
end

%%

PostP = PostprocessingMat_it;

%--------2020--------%
enda=toc(starta);
out=struct;
out.location_estimation=[Xiter amatrix];
out.estimated_parallel_time=enda/(K-m);
out.execution_date=date;
%--------------------%

end

function [PP,dtrue,m, Xiter] = MLESDP_inside_m1(PP, dtrue, dnoise, Adjacency, n, m)

%This is exactly the same of the MLESDP, but inside the loop. Some things
%have been deleted therefore. (For more details look at MLESDP function)

%%
%Declarations:

[i1,i2]=find(Adjacency(1:end-m,1:end-m));
Edges = size(i1,1);
[i3,i4]=find(Adjacency(1:end-m,end-m+1:end));
if size(i3,1) == 1
    i3=i3';
end
EdgesA = size(i3,1);


%position variable
X = sdpvar(2,n-m, 'full','real');

%noise
Dij = sdpvar(Edges,1, 'full','real');
dij = sdpvar(Edges,1, 'full','real');
Eik = sdpvar(EdgesA,1, 'full','real');
eik = sdpvar(EdgesA,1, 'full','real');

%Y=X'X as vec(Y) for nonzero element
Yoffdiag = sdpvar(Edges,1,'full','real');
Ydiag = sdpvar(n-m,1,'full','real');
costfunction = 0;

%Constraints and cost:
constraintset = (Dij(:)>=0) + (Eik(:)>=0);
if Edges
    ee=1:Edges;
    ii = i1(ee); jj = i2(ee); ij=ee;
    costfunction = sum(Dij)-2*dij'*diag(dnoise(ii,jj));
    M0 = Ydiag(ii)+Ydiag(jj) - 2*Yoffdiag(ij) - Dij(ij);
    constraintset = constraintset+ (M0(:)==0);
    for ee=1:Edges
        ii = i1(ee); jj = i2(ee); ij=ee;
        Xij = X(1:2, [ii,jj]);
        %Yij = [Ydiag(ii), Yoffdiag(ij); Yoffdiag(ij), Ydiag(jj)];
        constraintset = constraintset+...
            ([eye(2), Xij; Xij', [Ydiag(ii), Yoffdiag(ij); Yoffdiag(ij), Ydiag(jj)]]>=0);
        constraintset = constraintset+...
            ([1, dij(ij); dij(ij), Dij(ij)]>=0);
    end
end
ee=1:EdgesA;
ii = i3(ee); kk = n-m+i4(ee); ik=ee;
costfunction = costfunction+sum(Eik)-2*eik'*diag(dnoise(ii,kk));
M = Ydiag(ii)+diag(PP(1:2,kk)'*PP(1:2,kk)) - diag(2*X(1:2,ii)'*PP(1:2,kk)) - Eik(ik);
constraintset = constraintset+(M(:)==0);
for ee=1:EdgesA
    ik=ee;
    constraintset = constraintset+...
        ([1, eik(ik); eik(ik), Eik(ik)]>=0);
end

%%
%Solving

disp('Solving using SeDuMi...')

solvesdp(constraintset, costfunction);

Xiter = double(X);

end

function [Popt, Costopt, LambdaOpt, vectLam] = MLESDP_inside_v2(PP, dnoise, Adjacency, m,n)

%This is exactly the same of the MLESDP, but inside the loop and with the added
%equality constraint. It is the problem (15) of TSP. (For more details look at MLESDP function)

%%
%Declarations:

[i1,i2]=find(Adjacency(1:end-m,1:end-m));
Edges = size(i1,1);

AdjacencyT = Adjacency + Adjacency'; dnoiseT = dnoise+dnoise';
Yiimat_TR = sdpvar(n-m,n-m, 'full', 'real'); Yijmat_TR= sdpvar(Edges,n-m, 'full', 'real');
Dijmat_TR= sdpvar(Edges,n-m, 'full', 'real');dijmat_TR= sdpvar(Edges,n-m, 'full', 'real');
Xmat_TR = sdpvar(2*(n-m), n-m,'full', 'real');
YY_TR = [Yiimat_TR; Yijmat_TR; Dijmat_TR; dijmat_TR; Xmat_TR];
Zijmat_TR = sdpvar(size(YY_TR,1), Edges, 'full', 'real');
options = sdpsettings('verbose',0);
Eikmat_TR = sdpvar(m,n-m, 'full', 'real');
eikmat_TR = sdpvar(m,n-m, 'full', 'real');
LambdaM = zeros(size(YY_TR,1),2*(Edges));

%Cost/constraints
Cost = 0*sum(Eikmat_TR(:,1));
ConstraintSet = (Dijmat_TR(:)>=0) + (Eikmat_TR(:)>=0);
indexcount=1;
for ii=1:n-m
    
    ij=find(AdjacencyT(ii,1:end-m));
    ik=find(AdjacencyT(ii,end-m+1:end));
    mi = length(ik); ntoti = length(ij)+1;ni = ntoti+mi;
    loc_Adj = zeros(ni,ni); loc_dnoise = zeros(ni,ni);
    loc_i = sort([ii, ij, ik+n-m]); i_loc = find(loc_i == ii);
    loc_Adj(i_loc, 1:ni) = AdjacencyT(ii, loc_i);
    loc_Adj(1:ni, i_loc) = AdjacencyT(loc_i,ii);
    loc_dnoise(i_loc, 1:ni) = dnoiseT(ii, loc_i);
    loc_dnoise(1:ni, i_loc) = dnoiseT(loc_i,ii);
    
    vecX = [2*loc_i(1:ntoti)-1, 2*loc_i(1:ntoti)]; vecX = sort(vecX);
    vectori = [loc_i(1:ntoti), n-m+ij, n-m+Edges+ij, n-m+2*Edges+ij, ...
        n-m+3*Edges+vecX];
    
    
    Adjacencyins = triu(loc_Adj); dnoiseins = triu(loc_dnoise);
    [i1ins,i2ins]=find(Adjacencyins(1:end-mi,1:end-mi));
    Edgesins = size(i1ins,1);
    [i3ins,i4ins]=find(Adjacencyins(1:end-mi,end-mi+1:end));
    if size(i3ins,1) == 1
        i3ins = i3ins';
    end
    EdgesAins = size(i3ins,1);
    
    eeins=1:Edgesins;
    iiins = i1ins(eeins); jjins = i2ins(eeins); ijins=eeins;
    pvec = loc_i(1:ntoti); pvecx = [2*loc_i(1:ntoti)-1, 2*loc_i(1:ntoti)]; pvecx = sort(pvecx);
    
    if Edgesins
        Cost = Cost + 1/2*(sum(Dijmat_TR(ij,ii))-2*dijmat_TR(ij,ii)'*diag(dnoiseins(iiins,jjins))); %1/2 for the separability
        %Cost
        M0 = Yiimat_TR(pvec(iiins),ii)+Yiimat_TR(pvec(jjins),ii) - 2*Yijmat_TR(ij(ijins),ii) - Dijmat_TR(ij(ijins), ii);
        ConstraintSet = ConstraintSet+ (M0(:)==0);
    end
    
    
    for eeins=1:Edgesins
        iiins = i1ins(eeins); jjins = i2ins(eeins); ijins=eeins;
        Xij = [Xmat_TR(pvecx([iiins*2-1, iiins*2]),ii), Xmat_TR(pvecx([jjins*2-1, jjins*2]),ii)];
        vex = [2*[iiins jjins]-1, 2*[iiins jjins]]; vex = sort(vex);
        
        ConstraintSet = ConstraintSet+...
            ([eye(2), Xij; Xij', [Yiimat_TR(pvec(iiins),ii), Yijmat_TR(ij(ijins),ii); Yijmat_TR(ij(ijins),ii), Yiimat_TR(pvec(jjins),ii)]]>=0);
        ConstraintSet = ConstraintSet+...
            ([1, dijmat_TR(ij(ijins),ii); dijmat_TR(ij(ijins),ii), Dijmat_TR(ij(ijins),ii)]>=0);
        veci = [pvec(iiins),pvec(jjins), ni-mi+ijins,ni-mi+Edgesins+ijins,ni-mi+2*Edgesins+ijins,ni-mi+3*Edgesins+vex];
        
        ppp = find(((i1==ii)&(i2==ij(ijins)))|((i1==ij(ijins))&(i2==ii)));
        icount = ppp;
        
        vectZ = [veci(1:2),vectori(veci(3:end))];
        
        vectS = [Yiimat_TR(pvec(iiins),ii)', Yiimat_TR(pvec(jjins),ii)', Yijmat_TR(ij(ijins),ii)', ...
            Dijmat_TR(ij(ijins),ii)', dijmat_TR(ij(ijins),ii)', vec(Xij)']' - Zijmat_TR(vectZ, icount);
        
        ConstraintSet = ConstraintSet+((vectS == 0):'Equality');
        LambdaM(vectZ,ppp*(ii<ij(ijins))+(ppp+Edges)*(ii>ij(ijins))) = indexcount;
        indexcount=indexcount+1;
    end
    
    
    if EdgesAins
        for eeins=1:EdgesAins
            iiins = i3ins(eeins); kkins = ni-mi+i4ins(eeins); ikins=eeins;
            Cost = Cost+sum(Eikmat_TR(ik(ikins),ii))-2*eikmat_TR(ik(ikins),ii)'*diag(dnoiseins(iiins,kkins));
            Xanc = Xmat_TR(pvecx([2*iiins-1, 2*iiins]), ii);
            M = Yiimat_TR(pvec(iiins),ii)+diag(PP(1:2,loc_i(kkins))'*PP(1:2,loc_i(kkins))) - diag(2*Xanc'*PP(1:2,loc_i(kkins))) - Eikmat_TR(ik(ikins),ii);
            ConstraintSet = ConstraintSet+(M(:)==0);
            
            ConstraintSet = ConstraintSet+...
                ([1, eikmat_TR(ik(ikins),ii); eikmat_TR(ik(ikins),ii), Eikmat_TR(ik(ikins),ii)]>=0);
        end
    end
    
    
end

YY_TR = [Yiimat_TR; Yijmat_TR; Dijmat_TR; dijmat_TR; Xmat_TR];

%%
%Solve the problem:
solvesdp(ConstraintSet, Cost, options);

%ConstraintSet
Popt = double([Eikmat_TR; eikmat_TR; YY_TR]);
Costopt = double(Cost);

CC = (ConstraintSet('Equality')); num = length(CC);
LambdaOpt = 0*LambdaM;
for iii = 1:num
    ip=find(LambdaM==iii);
    LambdaOpt(ip) = dual(CC(iii)); %Get the dual variables associated with the equality constraints
end

vectLam = max(LambdaM);

end

function [Yi, pi, lambdaij] = solvemysdp(rho, zij, lambdaij, Yijold, PP, dnoise, Adjacency, n, m, accuracy)
%Function that sets up and solves the local problems. Highly based on
%ML-ESDP with a few tuning due to lambda and z. 

%%
%find some indeces

    [i1,i2]=find(Adjacency(1:end-m,1:end-m));
    Edges = size(i1,1);
    [i3,i4]=find(Adjacency(1:end-m,end-m+1:end));
    if size(i3,1) == 1
        i3 = i3';
    end
    EdgesA = size(i3,1);
    
%Declare some variables    
    
    %position variable
    X = sdpvar(2,n-m, 'full','real');
   
    %noise
    Dij = sdpvar(Edges,1, 'full','real');
    dij = sdpvar(Edges,1, 'full','real');
    Eik = sdpvar(EdgesA,1, 'full','real');
    eik = sdpvar(EdgesA,1, 'full','real');
    
    %Y=X'X as vec(Y) for nonzero element
    Yoffdiag = sdpvar(Edges,1,'full','real');
    Ydiag = sdpvar(n-m,1,'full','real');
    t = sdpvar(Edges,1,'full','real');
    
    %Cost and constraints
    costfunction = 0*sum(sum(Dij)) + 0*sum(sum(Eik));
    constraintset = (Dij(:)>=0) + (Eik(:)>=0);
    
    ee=1:Edges;
    ii = i1(ee); jj = i2(ee); ij=ee;
    if Edges
    costfunction = 1/2*(sum(Dij)-2*dij'*diag(dnoise(ii,jj))); %1/2 for the separability
    M0 = Ydiag(ii)+Ydiag(jj) - 2*Yoffdiag(ij) - Dij(ij);
    constraintset = constraintset+ (M0(:)==0);            
    end
    for ee=1:Edges
        ii = i1(ee); jj = i2(ee); ij=ee;
        Xij = X(1:2, [ii,jj]);
        vex = [2*[ii jj]-1, 2*[ii jj]]; vex = sort(vex);
        %Yij = [Ydiag(ii), Yoffdiag(ij); Yoffdiag(ij), Ydiag(jj)];
        constraintset = constraintset+...
            ([eye(2), Xij; Xij', [Ydiag(ii), Yoffdiag(ij); Yoffdiag(ij), Ydiag(jj)]]>=0);
        constraintset = constraintset+...
            ([1, dij(ij); dij(ij), Dij(ij)]>=0);
        veci = [ii,jj, n-m+ij,n-m+Edges+ij,n-m+2*Edges+ij,n-m+3*Edges+vex];
        
        %Using Schur complement (Sedumi is not as smart as you think!)
        vectS = [Ydiag(ii)', Ydiag(jj)', Yoffdiag(ij)', Dij(ij)', dij(ij)', vec(Xij)']' - zij(veci, ij);
        
        resij = lambdaij(veci, ij)+rho*(Yijold(veci) - zij(veci, ij)); 
        constraintset = constraintset + ([eye(length(vectS)), vectS; vectS', t(ij)]>=0) + (t(ij)>=0);
        costfunction = costfunction+resij'*[Ydiag(ii)', Ydiag(jj)', Yoffdiag(ij)', Dij(ij)', dij(ij)', vec(Xij)']'+rho/2*t(ij);
        lambdaij(veci, ij) = resij;
    end
    
    if EdgesA
        ee=1:EdgesA;
        ii = i3(ee); kk = n-m+i4(ee); ik=ee;
        
        costfunction = costfunction+sum(Eik)-2*eik'*diag(dnoise(ii,kk));
        M = Ydiag(ii)+diag(PP(1:2,kk)'*PP(1:2,kk)) - diag(2*X(1:2,ii)'*PP(1:2,kk)) - Eik(ik);
        constraintset = constraintset+(M(:)==0);
        
        for ee=1:EdgesA
            ik=ee;
            constraintset = constraintset+...
                ([1, eik(ik); eik(ik), Eik(ik)]>=0);
        end
    end
    Yi = [Ydiag', Yoffdiag', Dij', dij', vec(X)'];
    pi = [Eik', eik'];
    
%%
%Solving using sedumi..

%Set accuracy level, for safety anything below 1e-5, do as default
if accuracy>1e-5
        options = sdpsettings('verbose',0, 'solver','sedumi','sedumi.eps',accuracy);
    
else
        options = sdpsettings('verbose',0);
end
        
    solvesdp(constraintset, costfunction, options);
 
%%   
    %Output the solution:
    Yi = double(Yi); pi = double(pi);
        
end