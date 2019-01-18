% run_expr_doublePendulum.m
% perform row-type and column-type Hankel DMDs 
% Keisuke Fujii

clear; close all
dbstop if error
dataDir = './' ; 
ResDir = './' ;
Fs = 40 ; % sampling frequency

% parameter --------
rep = 6 ;  % length of cycles to analyze (max)
nval = 10 ; % number of sequences for validation
ntest = 1  ; % number of sequences for test, example:1 batch:10
Fig = 1 ; % example:1 batch:0
SS = 1 ; % number of participants: 1
C = 1 ; % number of conditions: 1

% load data
load([dataDir,'doublePendulum'])  

% preprocessing
for n = 1:50
    Ang_Seg{1,1}{n} = Ang(:,(n-1)*Fs*4+1:n*Fs*4)' ; % 160 time stamps
end

% segmentate data
for s = 1:SS
    for c = 1:C
        N = length(Ang_Seg{s,c}) ;
        nnn = 1 ;
        
        clear Ang Angrep
        % preprocessing
        for n = 1:N
            Ang{n} = Ang_Seg{s,c}{n} ;
            label_str{s,c}{n,1} = 'double pendulum';
        end
        nn = 1 ; nnv = 1 ;
        for n = 1:N
            if n <= N-rep
                if rep == 6 ; Angrep{n,1} = [Ang{n};Ang{n+1}(2:end,:);Ang{n+2}(2:end,:);Ang{n+3}(2:end,:);Ang{n+4}(2:end,:);Ang{n+5}(2:end,:)] ;
                    Trep{n,1} = cumsum([size(Ang{n},1),size(Ang{n+1},1)-1,size(Ang{n+2},1)-1,size(Ang{n+3},1)-1,size(Ang{n+4},1)-1,size(Ang{n+5},1)-1]) ;
                end
                if n <= nval
                    Xval{s,c}{nnv,1} = Angrep{n,1} ;
                    tmp = Xval{s,c}{nnv} ;  
                    Xval{s,c}{nnv,1} = tmp - repmat(mean(tmp,1),size(tmp,1),1) ;
                    Tval{s,c}{nnv,1} = Trep{n,1} ;
                    nnv = nnv + 1 ;
                else
                    if n < nval+rep+ntest && n >= nval+rep
                        X{s,c}{nn,1} = Angrep{n} ;
                        tmp = X{s,c}{nn} ; 
                        Ang_mean{nn} = tmp - repmat(mean(tmp,1),size(tmp,1),1) ;
                        X{s,c}{nn,1} = Ang_mean{nn} ;
                        Tte{s,c}{nn,1} = Trep{nn,1} ;
                        nn = nn + 1 ;
                    end
                end
            end
        end
    end
end



% DMD
datType = 'pendulum' ;
param.M = [1:2:5 8:4:20 25:5:40 50:10:80 100] ; % parameter m
param.limColumn = 20 ; % limit of column parameter
% HDMD for validation and convergence
if 0 % we did not perform validation procedure for pendulum data
    for s = 1:SS  
        for c = 1:C
            [errVal{s,c},optM{s,c},optR{s,c},CVerr{s,c}] = HDMD_val(Xval{s,c},Tval{s,c},Fs,param,rep) ;
            disp(['Validation of Hankel DMD s = ',num2str(s) ' c = ', num2str(c), ' finished' ]);
        end
    end
    save([dataDir,'Pendulum_val'],'errVal','optM','optR','CVerr','param')
else load([dataDir,'Pendulum_val'])
end

if 1% HDMD for all trials
    for s = 1:SS%%%
        for c = 1:C %
            optM{s,c} = [20 100] ;
            optR{s,c} = [20 20] ;
            methodstr = {'Exact DMD','Comp. DMD','Column HDMD','Row HDMD'};
            [Phin{s,c},Lamn{s,c},err{s,c},n_hankel(s,c),PhiP{s,c},recDMDs{c,1},...
                VAF{s,c},cumerror{s,c},timeDyns{c,1},Normn{s,c}] ...
                = HDMD_batch(X{s,c},Tte{s,c},methodstr,Fs,datType,Fig,optM{s,c},optR{s,c},rep) ;
            disp(['Test of Hankel DMD s = ',num2str(s) ' c = ', num2str(c), ' finished' ]);
            
        end
        recDMD = recDMDs ;
        save([ResDir,'Pendulum_Result_recDMD_',num2str(s)],'recDMD','-v7.3');
        timeDyn = timeDyns ;
        save([ResDir,'Pendulum_Result_timeDyn_',num2str(s)],'timeDyn','-v7.3');
    end
    save([dataDir,'Pendulum_result'],'Phin','PhiP','Lamn','err','n_hankel','VAF','cumerror','Normn')
    
else load([dataDir,'Pendulum_result'])
end



