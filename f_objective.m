function [f,R]=f_objective(x,snr)
%  x:   power assignment
%  f:   function value

N=length(x)+1;
snr=10.^(snr./10);

% multi-user interference
sum_alpha=zeros(1,N);
for i=1:N-1
    sum_alpha(i+1)=sum_alpha(i)+x(i);
end

% data rate of each user
R=zeros(1,N);
for i=1:N-1
    R(i)=log2(1+x(i)*snr(i)/(sum_alpha(i)*snr(i)+1));
end
R(N)=log2(1+(1-sum_alpha(N))*snr(N)/(sum_alpha(N)*snr(N)+1));

% objective_function
f=sum(log(R));
% f=var(R,1);
end