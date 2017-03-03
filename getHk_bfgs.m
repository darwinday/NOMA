function Hk = getHk_bfgs(S,Y,h)
% 
% This function returns the approximate inverse Hessian
% 
% Input
%   S:      Memory matrix (n by k), s{i} = x{i+1}-x{i}
%   Y:      Memory matrix (n by k), df{i} = df{i+1}-df{i}
%   h:      value of initial Hessian diagonal elements
% 
% Output
%   Hk:     the the approximate inverse Hessian
% 

    n = length(S);
    ro = 1/sum(Y.*S);
    
    temp = S*Y';
    Hk = (eye(n)-ro*temp) * h * (eye(n)-ro*temp') + ro*S*S';  

end % end of getHk_bfgs

