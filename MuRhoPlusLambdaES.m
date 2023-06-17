function nnet=MuRhoPlusLambdaES(IndivSize,net,P,T)
    BMu=initializeES(IndivSize,net,P,T);
    Lambda=BMu{2}.Lambda; Rho=BMu{2}.Rho;
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
        figure(1)
        plot(BMu{2}.bF(1:BMu{2}.g),'r'), hold on;
        plot(BMu{2}.mF(1:BMu{2}.g),'b'), hold off;        
    end
    nnet=SetFinalWeightsES(BMu{1}(:,1),net);
end