function [lambda,Phi,Norms,Phif] = DMD(X1,X2,eps_SVD,eps_norm) 
% Code is provided by: 
% Kutz et al. Dynamic Mode Decomposition: Data-Driven Modeling of Complex Systems, 2017

% and modified by Keisuke Fujii

% function [Phi,omega,lambda,b,Xdmd] = DMD(X1,X2,r,dt)
% Computes the Dynamic Mode Decomposition of X1, X2
%
% INPUTS: 
% X1 = X, data matrix
% X2 = X', shifted data matrix
% Columns of X1 and X2 are state snapshots 
% r = target rank of SVD
% dt = time step advancing X1 to X2 (X to X')
%
% OUTPUTS:
% Phi, the DMD modes
% omega, the continuous-time DMD eigenvalues
% lambda, the discrete-time DMD eigenvalues
% b, a vector of magnitudes of modes Phi
% Xdmd, the data matrix reconstrcted by Phi, omega, b

%% DMD
[U, S, V] = svd(X1, 'econ');
eval   = diag(S); %  
[~,IX] = sort(eval,'descend');
IX     = IX(abs(eval)>eps_SVD); %  

U_r = U(:, IX); % truncate to rank-r
S_r = S(IX, IX);
V_r = V(:, IX);
Atilde = U_r' * X2 * V_r *diag(1./diag(S_r)); % low-rank dynamics
[W_r, D, Z] = eig(Atilde);
Phi = X2 * V_r *diag(1./diag(S_r)) * W_r *diag(1./diag(D)); % DMD modes

lambda = diag(D); % discrete-time eigenvalues
% omega = log(lambda)/dt/2/pi; % continuous-time eigenvalues

if nargout>2 % norm
    if 1 % scaling Phi
        Ahat = (S_r^(-1/2)) * Atilde * (S_r^(1/2));
        [What, D] = eig(Ahat);
        W_r = S_r^(1/2) * What;
        Phif = X2*V_r/S_r*W_r;
    end
    Cnorm = (diag(Phif'*Phif));
%     Cnorm = sqrt(sum(abs(Phi).^2,1));  %% Euclidean norm of Koopman Modes
    [Norms, ind]=sort(Cnorm,'descend');      %% sorting based on the norm
    % ind     = ind(abs(Norms)>eps_norm); % absolute tolerance
    ind     = ind(abs(Norms)/abs(Norms(1))>eps_norm); % relative tolerance
    lambda= lambda(ind);           %% Koopman Eigenvalues
    Phi = Phi(:,ind);                %% Koopman Modes
    Phif = Phif(:,ind);  
    Norms = Norms(1:length(ind));
end

