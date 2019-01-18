function [lambda,Phi,Xaug,Norms,Phif,PhiP,PhiR] = Hankel_DMD( Data,n,m,param,type)
%HANKEL-Dynamic Mode Decompsition as presented in
% "Ergodic theory, dynamic mode decomposition and computation of Koopman
% spectral properties" by H. Arbabi and I. Mezic, arXiv:1611.06664

% modified by Keisuke Fujii

% inputs : 
% Data - the data matrix: each row is a time series data on an observable
% m,n: number of columns and rows, respectively, in the data Hankel
% matrices: size(Data,2)>m+n and preferably m>>n 

% or 

% ( Data,n,m,tol )- with Tol (optional) being the threshold for filtering thru SVD - the
% default value is 1e-10




% outputs: 
% 1 - HModes- Dynamic modes which approximate the Koopman eigenfunctions

% 2- HEvalues - Koopman or Dynamic Eigenvalues: the exponents and frequencies
% that make up the time-varaiation of data in time

% 3- Norms - The L2-norm of Koopman eigenfunction contribution to the
% observable

index1 = 1:n;
index2 = n:n+m-1;

X = []; Y=[];

for ir = 1:size(Data,1)
   
    % Hankel blocks ()
    c = Data(ir,index1).'; r = Data(ir,index2);
    H = hankel(c,r).';
    c = Data(ir,index1+1).'; r = Data(ir,index2+1);
    UH= hankel(c,r).';
    
    % scaling of Hankel blocks
%     if ir>1
%         alpha = norm(H(:,1))/norm(X(:,1));
%         H = alpha*H;
%         UH= alpha* UH;
%     end
    
    % the data matrices fed to exact DMD
    if strcmp(type,'row')
        X=[X,H]; Y=[Y,UH];
    elseif strcmp(type,'column')
        X=[X;H]; Y=[Y;UH];
    end
    
end
% X = X' ; Y = Y' ;
Xaug = X ;


%% DMD
[U, S, V] = svd(X, 'econ');
eval   = diag(S); 
[~,IX] = sort(eval,'descend');
if isempty(param.r_SVD)
    IX = IX(abs(eval)>param.eps_SVD);
else
    IXr = IX ; 
    IX = IX(1:param.r_SVD);
    IXr(IX) = [];
end

U_r = U(:, IX); % truncate to rank-r
S_r = S(IX, IX);
V_r = V(:, IX);

Atilde = U_r' * Y * V_r *diag(1./diag(S_r)); % low-rank dynamics
[W_r, D, Z] = eig(Atilde);
Phi = Y * V_r *diag(1./diag(S_r)) * W_r * diag(1./diag(D)); % Exact DMD modes
PhiP = U_r*W_r ; % Projected DMD modes
lambda = diag(D); % discrete-time eigenvalues

if nargout>3 % norm
    if 1 % scaling Phi
        Ahat = (S_r^(-1/2)) * Atilde * (S_r^(1/2));
        [What, D, Z] = eig(Ahat);
        W_r2 = S_r^(1/2) * What;
        Phif = Y*V_r/S_r*W_r2;
    end
    Cnorm = sqrt(sum(abs(Phif).^2,1))';  %% Euclidean norm of Koopman Modes
    [Norms, ind]=sort(Cnorm,'descend');  %% sorting based on the norm
    % ind     = ind(abs(Norms)>eps_norm); % absolute tolerance
    if isempty(param.r_SVD)
        ind     = ind(abs(Norms)/abs(Norms(1))>param.eps_norm); % relative tolerance
    end
    if isempty(ind)
        error('empty')
    end
    lambda= lambda(ind);           %% Koopman Eigenvalues
    Phi = Phi(:,ind);                %% Koopman Modes
    Phif = Phif(:,ind); 
    Norms = Norms(1:length(ind));
    PhiP = PhiP(:,ind); 

    if strcmp(type,'row')
        PhiR = pinv(Phi)*X ;
    else ; PhiR = [];
    end
end

end

