function initExpES(IndivSize,net,P,T)
%% Initialization
    for Mu=10:10:10
        for Rho=10:10:10
            for Lambda=10:10:50
                for Sigma=1:1:1
                    for MaxGens=100:100:100
                        for TRec=1:1:2
                            for TMut=1:1:5
                                tic;    
                                WCandidate=Sigma*randn(IndivSize,Mu);
                                z=Sigma*ones(IndivSize,Mu);
                                F=zeros(1,Mu);
                                G=2;   %number of constant gens
                                Gs=0;   %number of succesful gens
                                Ps=0;   %probability of sucess Gs/g
                                g=1;    %current generation    
                                psi=zeros(1,MaxGens);   %progress rate 
                                npsi=zeros(1,MaxGens);  %normalized progress rate 
                                psigma=zeros(1,MaxGens);%progress sigmma rate 
                                npsigma=zeros(1,MaxGens);%normalized progress sigmma rate 
                                bF=zeros(1,MaxGens);    %best fitness of generation 
                                mF=zeros(1,MaxGens);    %mean fitness of generation
                                tg=zeros(1,MaxGens);    %execution time for generation
                                if TRec==1
                                    TypeRec='inter';    %recombination type: 'inter'=intermediate-rho recomb.
                                else
                                    TypeRec='domin';    %                    'domin'=dominant-rho recomb.
                                end
                                switch TMut
                                    case 1
                                        TypeMut='isotr_1/5';%mutation type: 'isotr_1/5'=isotropic. 1/5th rule
                                    case 2
                                        TypeMut='isotr_mul';%               'isotr_mul'=iso. multiplicative
                                    case 3
                                        TypeMut='isotr_tpr';%               'isotr_tpr'=iso. two point rule
                                    case 4
                                        TypeMut='nonis_mul';%               'nonis_mul'=non-iso multiplicative
                                    case 5
                                        TypeMut='nonis_sco';%               'nonis_sco'=selective correlation
                                end
                                a=(rand*0.15)+0.85;    %param 0.85<a<1 for isotropic mutation sigma adaptation
                                tao=1/sqrt(IndivSize); %learning rate for isotropic multiplicative sigma adapt 1/sqrt(N) or 1/sqrt(2*N)
                                c=1;    %constant for non-isotropic multiplicative self adaptation
                                tao_ni=c/sqrt(2*IndivSize);         %learning parameters for non-isotropic
                                tao0_ni=c/sqrt(2*sqrt(IndivSize));  %multiplicative self adaptation 1/sqrt(2*N)
                                                                    %and 1/sqrt(2*sqrt(N)
                                CorrSize=100;       %size of the selective correlation    
                                sc=round(1+(rand(1,CorrSize)*(IndivSize-1)));   %initial selective components
                                %sC=corr(z(sc,:)')'; 
                                sC=eye(CorrSize);   %initial selective correlation matrix

                                %test fitness
                                for l=1:Mu
                                    %update input weights
                                    Layer1Dim=net.layers{1}.dimensions;
                                    IWs=zeros(size(net.IW{1}));
                                    for input=1:net.inputs{1}.size        
                                        IWs(:,input)=WCandidate((Layer1Dim*(input-1))+1:(Layer1Dim*input),l);
                                    end
                                    net.IW{1}=IWs;
                                    lastIdx=Layer1Dim*input;
                                    %update internal weights
                                    for layer=1:net.numlayers-1
                                        [r c]=size(net.LW{layer+1,layer});
                                        net.LW{layer+1,layer}=reshape(WCandidate(lastIdx+1:(r*c)+lastIdx,l),r,c);
                                        lastIdx=(r*c)+lastIdx;        
                                    end
                                    %update biases
                                    for layer=1:net.numlayers
                                        net.b{layer}=WCandidate(lastIdx+1:net.layers{layer}.dimensions+lastIdx,l);
                                        lastIdx=net.layers{layer}.dimensions+lastIdx;
                                    end
                                    Y=sim(net,P); e=T-Y; F(l)=sum(abs(e));
                                end    
                                bF(1)=min(F); mF(1)=mean(F); tg(1)=toc;
                                [OrdFit IdxFit]=sort(F);
                                sParams=struct('Mu',Mu,'Rho',Rho,'Lambda',Lambda,'Sigma',Sigma,'z',z(:,IdxFit),'G',G,...
                                    'Gs',Gs,'Ps',Ps,'g',g,'TypeRec',TypeRec,'TypeMut',TypeMut,'a',a,...
                                    'tao',tao,'tao_ni',tao_ni,'tao0_ni',tao0_ni,'CorrSize',CorrSize,...
                                    'sC',sC,'sc',sc,'MaxGens',MaxGens,'psi',psi,'npsi',npsi,...
                                    'psigma',psigma,'npsigma',npsigma,'bF',bF,'mF',mF,'tg',tg);
                                BMu={WCandidate(:,IdxFit), sParams, F(:,IdxFit)};



                            %% Main part  
                                while stop_criteriaES(BMu)
                                    tic;
                                    for l=1:Lambda           
                                       C{l}=reproductionES(BMu,Rho);
                                       s{l}=s_recombinationES(C{l},Rho);
                                       Hats{l}=s_mutationES(s{l});
                                       y{l}=y_recombinationES(C{l},Rho);
                                       Haty{l}=y_mmutationES(y{l},Hats{l});
                                       HatF{l}=FitnessES(Haty{l},net,P,T);
                                       Hats{l}=Hats{l}{1};  %save variances
                                    end
                                    HatBLambda={Haty,Hats,HatF};
                                    BMu=selectionES(HatBLambda,BMu);        
                                    BMu{2}.g=BMu{2}.g+1;        
                                    BMu{2}.bF(BMu{2}.g)=BMu{3}(1);
                                    BMu{2}.mF(BMu{2}.g)=mean(BMu{3});
                                    BMu{2}.tg(BMu{2}.g)=toc;        
                                end
                                fname=['M-' num2str(Mu) '-R-' num2str(Rho) '-L-' num2str(Lambda) '-S-' num2str(Sigma) '-MG-' num2str(MaxGens) '-TR-' num2str(TRec) '-TM-' num2str(TMut) '-.png'];
                                h1=figure(1);
                                subplot(1,2,1), plot(BMu{2}.mF(1:BMu{2}.g),'b.-'), hold on;
                                plot(BMu{2}.bF(1:BMu{2}.g),'r'), hold off;
                                h = legend('Mean Fitness','Best Fitness',2,'Location','NorthEast'); set(h,'Interpreter','none')
                                subplot(1,2,2), plot(BMu{2}.npsi(1:BMu{2}.g),BMu{2}.npsigma(1:BMu{2}.g),'.');        
                                title(fname); saveas(h1,fname); pause(0.1);
                                nnet=SetFinalWeightsES(BMu{1}(:,1),net);
                                fname=['Mu-' num2str(Mu) '-Rho-' num2str(Rho) '-Lambda-' num2str(Lambda) '-Sigma-' num2str(Sigma) '-MaxGens-' num2str(MaxGens) '-TRec-' num2str(TRec) '-TMut-' num2str(TMut) '-.mat'];
                                save(fname, 'nnet', 'BMu'); close all;
                                clear nnet BMu C s Hats y Haty HatF Hats
                            end
                        end
                    end
                end
            end
        end
    end                 
end