function [PP,dd]=generateD_revised_eps(m,N,R,max_deg,c,epsilon)

load('test10-500.mat','PP')
n=size(PP,1);
I=eye(n*(N+m));

PP=[PP(:,491:491+m-1) PP(:,1:N)]; % set PP with the required number of sensors for generateD to select the transformed sensor
PP=PP+epsilon*reshape(I(:,c+n*m),n,N+m); % transform the required coordinate of a non-anchor

[PP,dd]=generateD(PP,m,N+m,R,0,max_deg);
dd=dd'+dd;