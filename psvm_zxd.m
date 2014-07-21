function [w,gamma,nu]=psvm_zxd(A,d)
% version 1.1
% last revision: 01/24/03
%===============================================================================
% Usage:    [w,gamma,trainCorr, testCorr,cpu_time,nu]=psvm(A,d,k,nu,output,bal)
%
% A and d are both required, everything else has a default
% An example: [w gamma train test time nu] = psvm(A,d,10);
%
% Input:
% A is a matrix containing m data in n dimensions each.
% d is a m dimensional vector of 1's or -1's containing
% the corresponding labels for each example in A.
% k is k-fold for correctness purpose
% nu - the weighting factor.
%      -1 - easy estimation
%      0  - hard estimation
%      any other value - used as nu by the algorithm
%      default - 0
% output - indicates whether you want output
%
% If the input parameter bal is 1
% the algorithm weighs the classes depending on the
% number of points in each class and balance them.
% It is useful when  the number of point in each class
% is very unbalanced.
%
% Output:
% w,gamma are the values defining the separating
% Hyperplane w'x-gamma=0 such that:
%
% w'x-gamma>0 => x belongs to A+
% w'x-gamma<0  => x belongs to A-
% w'x-gamma=0 => x can belongs to both classes
% nu - the estimated or specified value of nu
%
% For details refer to the paper:
% "Proximal Support Vector Machine Classifiers"
% available at: www.cs.wisc.edu/~gfung
% For questions or suggestions, please email:
% Glenn Fung, gfung@cs.wisc.edu
% Sept 2001.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(A);
r=randperm(size(d,1));d=d(r,:);A=A(r,:);    % random permutation
%move one point in A a little if perfectly balanced
ma=A(d==1,:); mb=A(d==-1,:);
[s1 ~]=size(ma);
c1=sum(ma)/s1;
[s1 ~]=size(mb);
c2=sum(mb)/s1;
if (c1==c2)
    A(3,:)=A(3,:)+0.01*norm(A(3,:)-c1,inf)*ones(1,n);
end
% default values for input parameters
bal=0;
[H,v]=HV(A,d,bal);  % calculate H and v
nu = EstNuLong(H,d,m);
[w, gamma] = core(H,v,nu,n);
%%%%%%%%%%%%%%%% core function to calcuate w and gamma %%%%%%%%
function [w, gamma]=core(H,v,nu,n)
v=(speye(n+1)/nu+H'*H)\v;
w=v(1:n);gamma=v(n+1);
return
%%%%%%%%%%%%% EstNuLong %%%%%%%%%%%%
% use to estimate nu
function lamda=EstNuLong(H,d,m)
if m<201
    H2=H;d2=d;
else
    r=rand(m,1);
    [~,s2]=sort(r);
    H2=H(s2(1:200),:);
    d2=d(s2(1:200));
end
lamda=1;
[vu,u]=eig(H2*H2');u=diag(u);p=length(u);
yt=d2'*vu;
lamdaO=lamda+1;
cnt=0;
while (abs(lamdaO-lamda)>10e-4)&&(cnt<100)
    cnt=cnt+1;
    nu1=0;pr=0;ee=0;waw=0;
    lamdaO=lamda;
    for i=1:p
        nu1= nu1 + lamda/(u(i)+lamda);
        pr= pr + u(i)/(u(i)+lamda)^2;
        ee= ee + u(i)*yt(i)^2/(u(i)+lamda)^3;
        waw= waw + lamda^2*yt(i)^2/(u(i)+lamda)^2;
    end
    lamda=nu1*ee/(pr*waw);
end
return
%%% function to calculate H and v %%%%%%%%%%%%%
function [H,v]=HV(A,d,bal)
[m,~]=size(A);e=ones(m,1);
if (bal==0)
    H=[A -e];
    v=(d'*H)';
else
    H=[A -e];
    mm=e;
    m1=find(d==-1);
    mm(m1)=(1/length(m1));
    m2=find(d==1);
    mm(m2)=(1/length(m2));
    mm=sqrt(mm);
    N=spdiags(mm,0,m,m);
    H=N*H;
    %keyboard
    v=(d'*N*H)';
end