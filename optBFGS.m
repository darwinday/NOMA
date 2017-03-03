function [xk,fk,R,k,m] = optBFGS(myFun,x0,snr,maxIter,m)
% 
% Function optBFGS performs multivariate local optimization using the BFGS method.
% 
% Input
%   myFun:      the optimized function handle
%   x0:         vector of initial start
%   maxIter:    max number of iteration  
% Output
%   x1:         optimized variable
%   f1:         optimized function value
%   k:          iteration number
%
% Example
%   [optx,optf,k] = optBFGS(@myfun,x0,1000,)
% 

gradToler = 1e-3;  % tolerance for the norm of the slope
xToler = 1e-10;     % tolerance for the variables' refinement
k = 1;
n = length(x0);
Sk = zeros(n,1);
Yk = zeros(n,1);
Hk = eye(n);
xk = x0;

alpha_old = 0;

%     fprintf('%5s %15s %15s %15s\n', 'iter','step','fvalue','norm(g)');
%     fprintf('%5s %15s %15s %15s %15s %15s %15s\n', 'iter','step','a1','a2','a3','fvalue','norm(g)');
    [fk,gk]=feval(myFun,xk,snr);
    while true
      if k > maxIter
%           fprintf('Exceed the maxIter!\n');
          m = m+1;
        break;
      end
      
      % search direction       
      pk = - Hk*gk;
      % line search
      [alpha,fk1,gk1] = strongwolfe(myFun,pk,xk,snr,fk,gk,alpha_old,k);
      alpha_old = alpha;
      
      xk1 = xk+alpha*pk;
      sk = xk1-xk;
      yk = gk1-gk;
      Hk1 = getHk_bfgs(sk,yk,Hk);

      xk = xk1;
      gk = gk1;
      fk = fk1;
      Hk = Hk1;
      % save calculation
%       fprintf('%5d %15.4f %15.6f %15.6e\n',k,alpha,-fk,norm(gk));
%       fprintf(fid,'%.6f,',-fk);
%       fprintf('%5d %15.4f %15.6f %15.6f %15.6f %15.6f %15.4e\n',k,alpha,xk(1),xk(2),1-sum(xk),-fk,norm(gk));
      
      gnorm = norm(gk1);
      xnorm = norm(sk);
      if gnorm < gradToler
%           fprintf('Function Derivative below gradToler!\n');
          m = m;
          break;
      end
      if xnorm < xToler
%           fprintf('Function Value below xToler!\n');
          m = m;
        break;
      end
      k = k + 1;
    end
    [~,R]=f_objective(xk,snr);
    
end