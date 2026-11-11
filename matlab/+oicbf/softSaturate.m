function [y,dy] = softSaturate(v,L,U,symmetric)
%SOFTSATURATE Algebraically identical to source softplus, evaluated stably.
% The 1e-6 slope regularization is retained for Python parity.
k=30;
if symmetric
    offset=0; gamma=1/(tanh(.5*k*U)+1e-6);
else
    offset=(softplus(k*L)-softplus(-k*U))/k;
    s0=sigmoid(k*(-offset-L))-sigmoid(k*(-offset-U));
    gamma=1/(s0+1e-6);
end
a=k*(gamma*v-offset-L); b=k*(gamma*v-offset-U);
% Avoid subtracting two large softplus values on the positive tail.
if b>0
    y=U+(log1p(exp(-a))-log1p(exp(-b)))/k;
else
    y=L+(softplus(a)-softplus(b))/k;
end
dy=gamma*(sigmoid(a)-sigmoid(b));
end

function y=softplus(x)
y=max(x,0)+log1p(exp(-abs(x)));
end
function y=sigmoid(x)
if x>=0, y=1/(1+exp(-x)); else, e=exp(x); y=e/(1+e); end
end
