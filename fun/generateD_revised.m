function [PP,dd]=generateD_revised(m,N,R,max_deg)
%% INFORMATION

% DESCRIPTION: truncates the Benchamrk 500 network and generates sesnor 
%              adjacency matrix of the network with noised distances. Based
%              on the code generateD from https://web.stanford.edu/~yyye/Col.html

% INPUTS:
% * m - number of anchor sensors
% * N - total number of non-anchor sensors
% * R - radius
% * max_deg - maximum degree of neighbors for non-anchors (set it high enough so it has no effect, since maximum degree is dictated by the radius)

% OUTPUTS:
% * PP - real positions where last m columns are anchors
% * dd - sesnor adjacency matrix of the network with distances (noised or not)

load('test10-500.mat','PP') % loads the Benchmark 500 network

PP=[PP(:,491:491+m-1) PP(:,1:490)]; % moves m from the last columns of PP to the be the first m columns (a limitation of generateD) 

%% generateD function

% Input:
% * PP: 2xn matrix representing all point on 2D
% * N+m : the number of total points
% * m: the number of anchor points; the first m columns of PP
% * R: the radio range
% * sigma: noise factor (in this implementation should be set to 0, since
%          noises are genrated in create_realization)
% * max_deg: same as above

% Output
% * PP:2xn matrix representing all point on 2D, the last m columns are anchors
% * dd: the uppertriangle distance matrix from PP(:i) to PP(:,j), i<j.

[PP,dd]=generateD(PP,m,N+m,R,0,max_deg);

%%
dd=dd'+dd; % transofrms dd to a symmetric matrix (rather then upper triangular)