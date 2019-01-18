function [Phin,lamn,err,n_hankel,PhiP,recDMD,VAF,cumerror,timeDyn,Normn] = ...
    HDMD_batch(X_ori,T_te,methodstr,Fs,datType,Fig,optM,optR,rep)
N = length(X_ori) ; 
Nty = 2 ;
eps_SVD = 0 ; % for exact DMD
eps_norm = 0 ; % for exact DMD and companion-matrix DMD
eps_DMDerror = 1e-4 ; % for all DMD

nn = 1 ; nnn = 1 ;

for n = 1:N
    if 1
        n_hankel = T_te{n}(1) ;
        data = X_ori{n}(1:n_hankel,:)'; % transpose or not
        [D,T] = size(data) ;
        D3 = D ;
        X = data(:,1:end-1) ;
        Y = data(:,2:end) ;
            
        t = 1/Fs:1/Fs:(T+1)/Fs ;
        [lamn{1}{n,1},Phin{1}{n,1},Normn{n,1}] = DMD(X,Y,eps_SVD,eps_norm) ;
        [lamn{2}{n,1},Phin{2}{n,1},Normn{n,2}] = CompanionMatrix_DMD(data,eps_norm) ;
        nc = 2;% number of conventional DMDs
        D2 = [D 1] ;
        % Hankel DMD
        Type = {'column','row'} ;
        for ty = 1:Nty%
            if isempty(optM)
                m_hankel = n_hankel;%optM(ty) ;%
                m_hankel_1 = n_hankel;
            else ; m_hankel = optM(ty) ;%
                m_hankel_1 = optM(1);
            end
            if isempty(optR)
                param.r_SVD = 50;
            else; param.r_SVD = min([m_hankel*D2(ty),n_hankel,optR(ty)]) ;
            end
            [lamn{nc+ty}{n,1} ,xi{nc+ty}{n,1},Xaug{ty}{n,1},Normn{n,nc+ty},Phin{nc+ty}{n,1},PhiP{nc+ty}{n,1},PhiR{nc+ty}{n,1}] ...
                = Hankel_DMD(X_ori{n}',n_hankel,m_hankel,param,Type{ty}) ;
        end
        
        if nargout > 2
            for ty = 1:Nty+nc
                clear idx
                if ty <= nc
                    [timeDyn{n,ty},recDMD{n,ty},err(n,ty),cumerror{n,ty},VAF{n,ty},idx] = ...
                    DMDreconstruct(data,Phin{ty}{n,1},lamn{ty}{n,1},t,[]) ; 
                elseif ty == nc + 1 % column type HDMD
                    [timeDyn{n,ty},recDMD{n,ty},err(n,ty),cumerror{n,ty},VAF{n,ty},idx] = ...
                    DMDreconstruct(Xaug{ty-nc}{n},xi{ty}{n,1},lamn{ty}{n,1},t,[]) ;
                elseif ty == nc + 2  % row type HDMD
                    for d = 1:D
                        if d == 1 ; idx = []; end
                        int = (d-1)*n_hankel+1:d*n_hankel ;
                        Xaug2{d} = Xaug{ty-nc}{n}(:,int) ;
                        [timeDynd{n,ty,d},recDMDd{n}(:,:,:,d),errd(n,ty,d),cumerrord{n}(:,d),VAFd{n}(:,d),idx] = ...
                        DMDreconstruct(Xaug2{d},xi{ty}{n,1},lamn{ty}{n,1},t,idx) ;
                    end
                    err(n,ty) = mean(errd(n,ty,:)) ;
                end
                
                dd = 1 ; % row-type HDMD
                if ty <= nc + 2
                    if 1 
                        lamn{ty}{n,1} = lamn{ty}{n,1}(idx) ;
                        Normn{n,ty} = Normn{n,ty}(idx);              
                        if ty <= nc + 1 % column type 
                            if ty == nc + 1
                                PhiP{ty}{n,1} = PhiP{ty}{n,1}(:,idx) ;
                                Phin{ty}{n,1} = xi{ty}{n,1}(:,idx) ; % Phin=xi
                            else Phin{ty}{n,1} = Phin{ty}{n,1}(:,idx) ;    
                            end
                        elseif ty == nc + 2 % row type
                            timeDyn{n,ty} = timeDynd{n,ty,dd};%(:,tmp);
                            VAF{n,ty} = mean(VAFd{n},2);
                            recDMD{n,ty} = mean(recDMDd{n},4) ;
                            cumerror{n,ty} = mean(cumerrord{n},2) ;
                            PhiP{ty}{n,1} = Phin{ty}{n,1}(:,idx) ;
                            Phin{ty}{n,1} = PhiR{ty}{n,1}(idx,:)' ; % Phin=PhiR                    
                        end
                    end
                end
            end
        end
        disp(['Hankel dmd n=',num2str(n)]);
        
        % basic frequency (prior knowledge)
        if strcmp(datType,'pendulum') % double pendulum
            omTrue = sqrt([2-sqrt(2),2+sqrt(2)])*sqrt(9.81)/2/pi ;
        elseif strcmp(datType,'walk') % walk
            omTrue(1) = 1/(size(X_ori{n},1)/rep/Fs);
            omTrue(2:5) = omTrue(1).*[2:5] ;
        end
        
        % visualize time series
        if Fig
            for ty = 1:Nty+nc
                Mode{ty}{n,1} = [] ;
                if ty <=  nc
                    for r = 1:size(Phin{ty}{n,1},2)
                        Mode{ty}{n,1} = [Mode{ty}{n,1},abs(Phin{ty}{n,1}(:,r))',zeros(1,2)] ; % abs/real
                    end
                    if strcmp(datType,'walk')
                        R = min(size(Phin{ty}{n,1},2),6) ;
                    else ; % R = 4 ;
                        if ty == 1 ; R = 2 ; 
                        else ; R = 4 ;
                        end
                    end
                    D = size(Phin{ty}{n,1},1);
                else
                    D = D3 ;
                    if strcmp(datType,'pendulum') % double pendulum
                        R = 4 ;
                    elseif strcmp(datType,'walk') % walk
                        R = 6 ;
                    end
                    for r = 1:R
                        tmpP = [];
                        for d = 1:D
                            if ty == 3
                                tmpP(1,d) = abs(Phin{ty}{n,1}((d-1)*m_hankel_1+1,r)) ; % w/o delay Phin=xi
                            elseif ty == 4
                                tmpP(1,d) = abs(Phin{ty}{n,1}((d-1)*n_hankel+1,r)) ; % w/o delay Phin=PhiR
                            end
                        end
                        Mode{ty}{n,1} = [Mode{ty}{n,1},tmpP,zeros(1,2)] ; % abs/real
                    end
                end
                SP = 5 ;
                
                TimeDyn = real(timeDyn{n,ty}) ;
                figure(100+n) % compare with SVD
                nnn = ty+1 ;% compare with SVD
                subplot(5,SP,(nnn-1)*SP+1) % only DMD
                bar(Mode{ty}{n,1}(1:R*(D+2))) ;
                set(gca,'xtick',3:(D+2):(R-1)*(D+2)+3,'xticklabel',1:R);
                box off
                ylabel(methodstr{ty})  % compare with SVD
                if ty == 3 ; xlabel('mode'); end
                subplot(5,SP,(nnn-1)*SP+2)
                plot(TimeDyn(:,1),'b-'); hold on ;
                box off
                if ty == 3 ; xlabel('time (frame)'); end
                subplot(5,SP,(nnn-1)*SP+3)
                if R >= 3 ; plot(TimeDyn(:,3),'b-'); hold on ; end% r
                box off
                if ty == 3 ; xlabel('time (frame)'); end
                subplot(5,SP,(nnn-1)*SP+4)
                if R >= 5 ; plot(TimeDyn(:,5),'b-'); hold on ; end
                box off
                if ty == 3 ; xlabel('time (frame)'); end
            end
        end
        
        % Visualize eigenvalues and DMD modes
        if Fig
            methodstr2 = {'Exact DMD','Comp. DMD','column-type HDMD','row-type HDMD'};
            for ty = 1:Nty+nc
                colormap jet
                figure(n*1000)
                subplot(2,6,(ty-1)*3+1)
                % eigenvalue
                plot(real(lamn{ty}{n,1}),imag(lamn{ty}{n,1}),'bo') ;
                rectangle('Position', [-1 -1 2 2], 'Curvature', 1, ...
                    'EdgeColor', 'k', 'LineStyle', '--');
                axis([-1.2 1.2 -1.2 1.2]);
                axis square;
                xlabel('Re(\lambda)') ;
                ylabel('Im(\lambda)') ;
                title(methodstr2{ty}) 
                
                subplot(2,6,(ty-1)*3+2:(ty-1)*3+3)
                omega{ty,n} = log(lamn{ty}{n,1})*Fs/2/pi ;
                f = abs(imag(omega{ty,n}));
                maxP = max(Normn{n,ty}) ;
                P = (Normn{n,ty}./maxP);
                stem(f, P, 'b'); hold on
                
                % DFT
                if ty <= nc
                    XX = X_ori{n}' ;
                else XX = Xaug{ty-nc}{n};
                end
                timesteps = size(XX, 2);
                nelectrodes = size(XX, 1);
                NFFT = 2^nextpow2(timesteps);
                ff = Fs/2*linspace(0, 1, NFFT/2+1);
                clear fftp
                for c = 1:nelectrodes
                    fftp(c,:) = fft(XX(c,:), NFFT);
                end
                tmp = 2*abs(mean(fftp(c,1:NFFT/2+1), 1)) ;
                plot(ff, tmp./max(tmp), ...
                    'Color', 0.6*[1 1 1]);
                xlabel('Frequency(Hz)') ;
                ylabel('norm(normalized)') ;
                ylim([0 1]);
                if strcmp(datType,'pendulum') % double pendulum
                    line([omTrue(1) omTrue(1)], [0 1], 'Color', 'k', 'LineStyle', '--');
                    line([omTrue(2) omTrue(2)], [0 1], 'Color', 'k', 'LineStyle', '--');
                    xlim([-0.1 2]);
                elseif strcmp(datType,'walk') % walk
                    for o = 1:length(omTrue)
                    line([omTrue(o) omTrue(o)], [0 1], 'Color', 'c', 'LineStyle', '--');
                    end
                    xlim([-1 6]);
                end
                
            end
        end
        
    end
end


