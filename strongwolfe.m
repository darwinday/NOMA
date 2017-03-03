function [alphas,fs,gs] = strongwolfe(myFun,d,x0,snr,fx0,gx0,alpha_old,k)
% 
% Function strongwolfe performs Line search satisfying strong Wolfe conditions
% 
% Input
%   myFun:      the optimized function handle
%   d:          the direction we want to search
%   x0:         vector of initial start
%   fx0:        the function value at x0
%   gx0:        the gradient value at x0
% 
% Output
%   alphas:     step size
%   xs          x0+alphas*d
%   fs:         the function value at x0+alphas*d
%   gs:         the gradient value at x0+alphas*d

maxIter = 5;

c1 = 1e-4;
c2 = 0.9;
gx0 = gx0'*d;
fxp = fx0;
gxp = gx0;

i=1;

temp = x0./d;
temp1 = max(temp(temp<0));
temp2 = (1-sum(x0))/sum(d);
if temp1 <= 0
    alphamax = min(abs(temp1),abs(temp2));
else
    alphamax = abs(temp2);
end

alphamin = 0;

if (abs(gx0) > 1/5)|(k == 1)
    alphax = alphamax*0.8;     
%     alphax = sum(x0)/abs(sqrt(abs(gx0)))/4;   
elseif (abs(gx0) > 1/10)
    alphax = min(alpha_old*0.5,alphamax*0.8);
else
    alphax = min(1,alphamax*0.8);
end
    


% alphap is alpha_{i-1}
% alphax is alpha_i
% alphas is what we want.
    while true
        xx = x0+alphax*d;
        [fxx,gxx] = feval(myFun,xx,snr);
        fs = fxx;
        gs = gxx;
        gxx = gxx'*d;
        if (fxx > fx0 + c1*alphax*gx0) | ((i>1)&(fxx>=fxp))
            [alphas,fs,gs] = Zoom(myFun,x0,snr,d,0,alphax);
            return;
        end
        if abs(gxx) <= -c2*gx0,
            alphas = alphax;
            return;
        end
        if gxx >= 0,
            [alphas,fs,gs] = Zoom(myFun,x0,snr,d,alphax,0);
            return;
        end
        fxp = fxx;
        gxp = gxx;
        
        if i > maxIter
            alphas = alphax;
%             fprintf('Exceed the maxIter in Strongwolfe!');
            return
        end  
        % r = rand(1);%randomly choose alphax from interval (alphap,alpham)
%         alphax = min(alphax*2,0.8*alphamax);
        if alphax*1.2 >= alphamax
            alphax = alphax+0.8*(alphamax-alphax);
        else
            alphax = alphax*1.2;
        end
        i = i+1;
        end
end