function [err_DMD,optM,optR,optM_err] = HDMD_val(X_val,T_val,Fs,param,rep)
% coded by Keisuke Fujii

D = size(X_val{1},2) ;
M = param.M ;
[M indM] = sort(M) ;%
[R indR] = sort(M) ;%
I = size(X_val,1) ;
Nty = 2 ;
Nn = 1 ;
err_DMD = NaN(length(M),length(R),length(Nn),I,Nty) ;

% Hankel DMD
Type = {'column','row'} ;
D2 = [D 1] ;
for i = 1:I
    N = T_val{i}(1:Nn) ;
    for ty = 1:Nty%
        for n = 1:Nn
            n_hankel = N(n) ;
            for m = 1:length(M)
                m_hankel = M(m) ;
                if (~(ty ==1 && m_hankel>param.limColumn))&&(~(ty ==2 && m_hankel>n_hankel*(rep-2)))
                t{n,i} = 0:1/Fs:(n_hankel+1)/Fs ;
                for r = 1:length(R)
                    if R(r)<= min(m_hankel*D2(ty),n_hankel)
                        param.r_SVD = R(r) ;
                        [lamn{m,n}{i}{r,ty} ,xi{m,n}{i}{r,ty},Xaug{m,n}{i}{r,ty}] ...
                            = Hankel_DMD(X_val{i}',n_hankel,m_hankel,param,Type{ty}) ;
                    end
                end
                end
            end
        end
    end
end

for i = 1:I
    N = T_val{i}(1:Nn) ;
    for ty = 1:Nty%
        for n = 1:length(N)
            n_hankel = N(n) ;
            for m = 1:length(M)
                m_hankel = M(m) ;
                if (~(ty ==1 && m_hankel>param.limColumn))&&(~(ty ==2 && m_hankel>n_hankel*(rep-2)))
                for r = 1:length(R)
                    if R(r)<= min(m_hankel*D2(ty),n_hankel)
                        r_SVD = R(r) ;
                        if ty == 1 % column type HDMD
                            [~,~,err_DMD(m,r,n,i,ty)] = ...
                                DMDreconstruct(Xaug{m,n}{i}{r,ty},xi{m,n}{i}{r,ty},lamn{m,n}{i}{r,ty},t{n,i}) ; % DMDreconstruct
                        elseif ty == 2  % row type HDMD
                            for d = 1:D
                                int = (d-1)*n_hankel+1:d*n_hankel ;
                                Xaug2{m,n}{i,d}{r,ty} = Xaug{m,n}{i}{r,ty}(:,int) ;
                                [~,~,errd(m,r,n,i,ty,d)] = ...
                                    DMDreconstruct(Xaug2{m,n}{i,d}{r,ty},xi{m,n}{i}{r,ty},lamn{m,n}{i}{r,ty},t{n,i}) ;
                            end
                            err_DMD(m,r,n,i,ty) = mean(errd(m,r,n,i,ty,:)) ;
                        end
                        
                    end
                end
                end
            end
        end
    end
end

% Cross validation

ind0 = 1:I;
for i = 1:I
    ind = ind0(ind0~=i) ;
    for n = 1:length(N)
        for ty = 1:Nty
            % optM
            tmp = squeeze(nanmin(nanmin(err_DMD(:,:,n,ind,ty),[],2),[],4)) ;
            [err_CV_m(n,i,ty),ind_m(n,i,ty)] = nanmin(tmp,[],1) ;
        end
    end
end
for n = 1:length(N)
    for ty = 1:Nty
        M2 = unique(ind_m(n,:,ty)) ;
        for m = 1:length(M2)
            [err_Vm{n,ty}(m),indR2{n,ty}(m)] = nanmin(nanmean(err_DMD(M2(m),:,n,:,ty),4),[],2) ;
        end
        [optM_err(n,ty),tmpind] = min(err_Vm{n,ty});
        optM(n,ty) = M(M2(tmpind)) ;
        optR(n,ty) = R(indR2{n,ty}(tmpind)) ;
    end
end

