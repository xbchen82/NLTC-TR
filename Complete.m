function [iter,X,errList] = Complete(Tsr,w,lambda,beta,p,rho,epsilon,maxIter)
X = Tsr;
Omega = ~isnan(Tsr);
X(logical(1-Omega)) = mean(Tsr(Omega));
errList = zeros(maxIter,1);
dim = size(Tsr);
n1 = dim(1);
n2 = dim(2);
n3 = dim(3);
Tdim = ndims(Tsr);
M = cell(Tdim, 1);
Y1 = M;
for i = 1:Tdim
    M{i} = X;
    Y1{i} = zeros(dim);
end
Y4 = zeros(dim);
H = -eye(n2);
for i = 2:n2-1
    a = [i-1 i+1];
    H(i,a) = [1/2 1/2];
end
H(1,2) = 1;
H(n2,n2-1) = 1;
% H = dctmtx(n2)*H;
Msum = zeros(dim);
Y1sum = zeros(dim);
for iter = 1: maxIter
    rho =rho * 1.05;
    % update Mn
    ss = 0;
    Msum = 0*Msum;
    Y1sum = 0*Y1sum;
    for i = 1:Tdim
        mat1 = Unfold(X+Y1{i}/rho,dim,i);
        mat2 = Unfold(M{i},dim,i);
        tmp1 = dct2(mat1);
        tmp2 = dct2(mat2);
        M1_d = Fold(tmp1,dim,i);
        M2_d = Fold(tmp2,dim,i);
        Temp_X2 = Unfold(X,dim,i); %展开成矩阵，i维度展开
        temp = dct2(Temp_X2);
        X2 = Fold(temp,dim,i);%折回张量
        for u = 1:dim(i)
            if i == 1
                [U, M_S, V] = svd(squeeze(M1_d(u, :, :)), 'econ');
                d_m = svd(squeeze(M2_d(u, :, :)));
                [~, X2_D, ~] = svd(squeeze(X2(u, :, :)), 'econ');
            elseif i == 2
                [U, M_S, V] = svd(squeeze(M1_d(:, u, :)), 'econ');
                d_m = svd(squeeze(M2_d(:, u, :)));
                [~, X2_D, ~] = svd(squeeze(X2(:, u, :)), 'econ');
            else
                [U, M_S, V] = svd(M1_d(:, :, u), 'econ');
                d_m = svd(M2_d(:, :, u));
                [~, X2_D, ~] = svd(squeeze(X2(:, :, u)), 'econ');
            end
            s_m = diag(M_S);
            
            b = diag(X2_D);%取出X2_D从大到小的列向量
            ss = ss+sum(b);
            
            
            f = lambda*p*d_m.^(p-1)./(1+d_m.^p);
            f(d_m<10^-9) = Inf;
            temp = s_m-w(i)/rho*f;
            temp(temp<0) = 0;
            M_shrink = U * diag(temp) * V';
            if i == 1
                M_tmp(u, :, :) = M_shrink;
            elseif i == 2
                M_tmp(:, u, :) = M_shrink;
            else
                M_tmp(:, :, u) = M_shrink;
            end
        end
        mat = Unfold(M_tmp,dim,i);
        tmp = idct2(mat);
        M{i} = Fold(tmp,dim,i);
        Y1sum = Y1sum + Y1{i};
        Msum = Msum + M{i};
    end
    % update D
    tmp_1 = beta*H'*H+rho*eye(size(H'*H));
    tmp_2 = Unfold(rho*X+Y4,dim,2);
    D_mat = tmp_1\tmp_2;
    D = Fold(D_mat,dim,2);
    % update X
    oldX = X(logical(1-Omega));
    X = (Msum-Y1sum/rho+D-Y4/rho)/4;
    X(Omega) = Tsr(Omega);
    newX = X(logical(1-Omega));        
    % update Y
    for i = 1:Tdim
        Y1{i} = Y1{i} + rho*(X-M{i});
    end
    Y4 = Y4 + rho*(X-D);
    errList(iter) = norm(newX-oldX) / norm(oldX);
    if errList(iter) < epsilon
        break;
    end
end
    