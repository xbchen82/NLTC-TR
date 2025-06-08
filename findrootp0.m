function x = findrootp0(a, r, p)
[n1,n2] = size(a);
x = zeros(n1,n2);
if p == 1
    idx = find(a>r);
    x(idx) = a(idx)-r;
elseif p == 0
    idx1 = find(a>sqrt(2*r));
    x(idx1) = a(idx1);
else
    v = (r*p*(1-p))^(1/(2-p))+eps;
    v1 = v+r*p*v^(p-1);
    ob0 = 0.5*a.^2;
    idx2 = find(a<=v1);
    x = a;
    for k = 1:10
        tmp = x.^(p-1);
        f = (x-a) + r*p*tmp;
        g = 1-r*p*(1-p)*(tmp./x);
        x = x-f./g;
    end
    ob1 = 0.5*(x-a).^2 + r*x.^p;
    A = cat(2,zeros(size(a)),x);
    B = cat(2,ob0,ob1);
    [~,idx3] = min(B,[],2);
    A = A';
    idx3 = idx3';
    num = 0:n1-1;
    idx4 = idx3+num*2;
    x = A(idx4);
    x(idx2) = 0;
end
end