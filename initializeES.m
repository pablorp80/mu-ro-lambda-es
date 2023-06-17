function BMu=initializeES(IndivSize, net, P, T)
    tic;
    Mu=10;
    Rho=10;
    Lambda=30;
    Sigma=1;
    WCandidate=Sigma*randn(IndivSize,Mu);
    z=Sigma*ones(IndivSize,Mu);
    F=zeros(1,Mu);
    G=2;   %number of constant gens
    Gs=0;   %number of succesful gens
    Ps=0;   %probability of sucess Gs/g
    g=1;    %current generation
    MaxGens=100;    %total maximum of generations to stop
    psi=zeros(1,MaxGens);   %progress rate 
    npsi=zeros(1,MaxGens);  %normalized progress rate 
    psigma=zeros(1,MaxGens);%progress sigmma rate 
    npsigma=zeros(1,MaxGens);%normalized progress sigmma rate 
    bF=zeros(1,MaxGens);    %best fitness of generation 
    mF=zeros(1,MaxGens);    %mean fitness of generation
    tg=zeros(1,MaxGens);    %execution time for generation
    TypeRec='domin';    %recombination type: 'inter'=intermediate-rho recomb.
                        %                    'domin'=dominant-rho recomb.
    TypeMut='isotr_1/5';%mutation type: 'isotr_1/5'=isotropic. 1/5th rule
                        %               'isotr_mul'=iso. multiplicative
                        %               'isotr_tpr'=iso. two point rule
                        %               'nonis_mul'=non-iso multiplicative
                        %               'nonis_sco'=selective correlation
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
    s=struct('Mu',Mu,'Rho',Rho,'Lambda',Lambda,'Sigma',Sigma,'z',z(:,IdxFit),'G',G,...
        'Gs',Gs,'Ps',Ps,'g',g,'TypeRec',TypeRec,'TypeMut',TypeMut,'a',a,...
        'tao',tao,'tao_ni',tao_ni,'tao0_ni',tao0_ni,'CorrSize',CorrSize,...
        'sC',sC,'sc',sc,'MaxGens',MaxGens,'psi',psi,'npsi',npsi,...
        'psigma',psigma,'npsigma',npsigma,'bF',bF,'mF',mF,'tg',tg);
    BMu={WCandidate(:,IdxFit), s, F(:,IdxFit)};
end