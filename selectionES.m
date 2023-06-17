function NewBMu=selectionES(HatBLambda,BMu)
    %plus selection (Mu+Lambda)
    TotSols=zeros(size(BMu{1},1),(BMu{2}.Mu+BMu{2}.Lambda));
    TotSigm=zeros(size(BMu{2}.z,1),(BMu{2}.Mu+BMu{2}.Lambda));
    TotFitn=zeros(size(BMu{3},1),(BMu{2}.Mu+BMu{2}.Lambda));
    TotSols(:,1:BMu{2}.Mu)=BMu{1};
    TotSigm(:,1:BMu{2}.Mu)=BMu{2}.z;
    TotFitn(:,1:BMu{2}.Mu)=BMu{3};
    TotSols(:,(BMu{2}.Mu+1):(BMu{2}.Mu+BMu{2}.Lambda))=cell2mat(HatBLambda{1});
    TotSigm(:,(BMu{2}.Mu+1):(BMu{2}.Mu+BMu{2}.Lambda))=cell2mat(HatBLambda{2});
    TotFitn(:,(BMu{2}.Mu+1):(BMu{2}.Mu+BMu{2}.Lambda))=cell2mat(HatBLambda{3});
    [OrdFit IdxFit]=sort(TotFitn);
    NewBMu=BMu;
    NewBMu{1}=TotSols(:,IdxFit(1:BMu{2}.Mu));
    NewBMu{2}.z=TotSigm(:,IdxFit(1:BMu{2}.Mu));
    NewBMu{3}=TotFitn(:,IdxFit(1:BMu{2}.Mu));
    
    %parameter update
    if min(NewBMu{3})<min(BMu{3})
        NewBMu{2}.Gs=NewBMu{2}.Gs+1;
    end
    NewBMu{2}.Ps=NewBMu{2}.Gs/NewBMu{2}.g;
    NewBMu{2}.sc=round(1+(rand(1,NewBMu{2}.CorrSize)*(size(BMu{1},1)-1)));
    NewBMu{2}.sC=corr(NewBMu{2}.z(NewBMu{2}.sc,:)')';
    R=mean(sum((BMu{1}-repmat(NewBMu{1}(:,1),1,size(BMu{1},2))).^2).^0.5);
    delta=R-mean(sum((NewBMu{1}-repmat(NewBMu{1}(:,1),1,size(NewBMu{1},2))).^2).^0.5);
    NewBMu{2}.psi(NewBMu{2}.g)=delta*(size(NewBMu{1},1)/R);
    NewBMu{2}.npsi(NewBMu{2}.g)=NewBMu{2}.psi(NewBMu{2}.g)*(size(BMu{2}.z,1)/R);
    NewBMu{2}.psigma(NewBMu{2}.g)=mean(mean(NewBMu{2}.z));
    NewBMu{2}.npsigma(NewBMu{2}.g)=NewBMu{2}.psigma(NewBMu{2}.g)*(size(BMu{2}.z,1)/R);    
end
