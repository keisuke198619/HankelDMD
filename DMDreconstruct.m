function [time_dmd,Xdm,err,cumerror,VAF,idx] = DMDreconstruct(X,Phi,lambda,t,idx)
% coded by Keisuke Fujii

x1 = X(:, 1);
b = pinv(Phi)*x1;  
r = length(lambda) ;
dt = t(2) - t(1);
omega = log(lambda)/dt;
T = size(X,2) ;
time_dmd = zeros(T-1,r);
for iter = 1:T-1
    for p = 1:r
        time_dmd(iter,p) = b(p)*(exp(omega(p)*t(iter)));% time dynamics
        Xdm(:,iter,p) = real(Phi(:,p)*(b(p).*exp(omega(p)*t(iter))));% reconstruct
    end
end
    
X_dmd = sum(Xdm,3) ;
err = mean(mean(sqrt(abs(X(:,1:end-1)-X_dmd).^2),2)) ; % absolute error

if nargout > 3
    for p = 1:r
        VAF(p,1) = 1-sum(sum(abs(X(:,1:end-1)-Xdm(:,:,p)).^2))/sum(X(:).^2) ;
    end
    if isempty(idx)
        [val,idx] = sort(VAF,'descend') ;
    end
    time_dmd = time_dmd(:,idx);
    VAF = VAF(idx);
    Xdm = Xdm(:,:,idx) ;
    
    cum_Xdm = cumsum(Xdm,3) ;
    for p = 1:r
        cumerror(p,1) = mean(mean(sqrt(abs(X(:,1:end-1)-cum_Xdm(:,:,p)).^2),2)) ; % absolute error
    end
end

