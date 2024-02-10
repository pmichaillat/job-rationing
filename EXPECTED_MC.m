%%=============================================
%% Construct sequence of vectors of conditional expectations of a variable following a Markov chain 
%% with state space A (row) and transition matrix PI
%% Column t is expectation at t-1, for t<=num+1
%% Row s is conditional expectation in state s
%=============================================

function [EM]=EXPECTED_MC(PI,A,num);

EM=[A'];%current state
for i= 1:num
  EM=[EM,PI*EM(:,end)];
end