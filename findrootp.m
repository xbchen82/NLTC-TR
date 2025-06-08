function x = findrootp(y, lambda, p)
[n1,n2] = size(y);
y1 = abs(y);
x = zeros(n1,n2);
if p == 1
    idx = find(y1>lambda);
    tmp = sign(y).*(y1-lambda);
    x(idx) = tmp(idx);
elseif p == 0
    idx1 = find(y1>sqrt(2*lambda));
    x(idx1)=y(idx1);
else    
    gst = (2*lambda*(1-p))^(1/(2-p))+lambda*p*(2*lambda*(1-p))^((p-1)/(2-p));
    idx2 = find(y1>gst);
    x(idx2) = y1(idx2);
    for k = 1:10
        tmp1 = x(idx2);
        x(idx2) = y1(idx2)-lambda*p*tmp1.^(p-1);
    end
    x = sign(y).*x;
end
end


% function sigma = findrootp(alpha, lambda, p)
% [n1,n2] = size(alpha);
% a = (lambda*p*(1-p)/2)^(1/(2-p))+eps;
% b = 2*a-2*alpha+lambda*p*a^(p-1);
% idx = find(b<0);
% idx1 = find(b>=0);
% sigma = zeros(n1,n2);
% % sigma1 = zeros(n1,n2);
% % ob1 = zeros(n1,n2);
% sigma(idx) = 2*a;
% for i = 1:10
%     tmp = exp((p-1)*log(sigma));
%     f = 2*(sigma-alpha)+lambda*p*tmp;
%     g = 2+lambda*p*(p-1)*(tmp./sigma);
%     sigma = sigma-f./g;
% end
% sigma1 = sigma;
% ob1 = sigma.^2-2*sigma.*alpha+lambda*(abs(sigma)).^p;%exp((p)*log((abs(sigma))));          %(abs(sigma)).^p;
% sigma1(idx1) = 1;
% ob1(idx1) = inf;
% 
% a = -a;
% b = 2*a-2*alpha+lambda*p*abs(a)^(p-1);
% idx = find(b>0);
% idx1 = find(b<=0);
% % sigma2 = zeros(n1,n2);
% % ob2 = zeros(n1,n2);
% sigma(idx) = 2*a;
% for i = 1:10
%     tmp = exp((p-1)*log((abs(sigma))));
%     f = 2*(sigma-alpha)-lambda*p*tmp;%exp((p-1)*log((abs(sigma))));          %(abs(sigma)).^(p-1);
%     g = 2+lambda*p*(p-1)*(tmp./(abs(sigma)));%exp((p-2)*log((abs(sigma))));          %(abs(sigma)).^(p-2);
%     sigma = sigma-f./g;
% end
% sigma2 = sigma;
% ob2 = sigma.^2-2*sigma.*alpha+lambda*(abs(sigma)).^p;%exp((p)*log((abs(sigma))));          %(abs(sigma)).^p;
% sigma2(idx1) = 1;
% ob2(idx1) = inf;
% 
% A = cat(3,zeros(size(alpha)),sigma1,sigma2);
% B = cat(3,zeros(size(alpha)),ob1,ob2);
% [~,idx] = min(B,[],3);
% T_A = (reshape(A,[n1*n2,3]))';
% t_idx = (reshape(idx,[n1*n2,1]))';
% num = (0:n1*n2-1)*3;
% tem5 = t_idx+num;
% tem6 = T_A(tem5);
% sigma = reshape(tem6,[n1,n2]);
% end