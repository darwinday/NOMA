function [alphas,fs,gs] = Zoom(myFun,x0,snr,d,alphal,alphah)
c1=1e-4;
c2=0.9;
i=0;
maxIter=5;

    while true
       % interpolate
       xl = x0 + alphal*d;
       xh = x0 + alphah*d;
       [fx0,gx0] = feval(myFun,xl,snr);
       [fx1,gx1] = feval(myFun,xh,snr);
       gx0 = gx0'*d;
       gx1 = gx1'*d;
       
%        alphax = -(gx0*(alphah^2))/(2*(fx1-fx0-alphah*gx0));
       alphax = (gx0*alphah-gx1*alphal)/(gx0-gx1);       
       if (alphax > max(alphal,alphah))|(alphax < min(alphal,alphah))
           alphax = 0.5*(alphal+alphah);    
       end
%        alphax = 0.5*(alphah+alphal);
       alphas = alphax;
       xx = x0 + alphax*d;
       [fxx,gxx] = feval(myFun,xx,snr);
       fs = fxx;
       gs = gxx;
       gxx = gxx'*d;

       if ((fxx > fx0+c1*alphax*gx0) | (fxx>=fx0)),
          alphah = alphax;
       else
          if abs(gxx) <= -c2*gx0,
            alphas = alphax;
            return;
          end
          if gxx*(alphah-alphal) >= 0
            alphah = alphal;
          end
          alphal = alphax;
       end
         i = i+1;
       if i > maxIter
          alphas = alphax;
%           fprintf('Exceed the maxIter in Zoom1!');
          return
       end
    end
end