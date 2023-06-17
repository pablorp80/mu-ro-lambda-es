function yl=y_recombinationES(Cl,Rho)
    yl=zeros(size(Cl{1},1),1);
    switch Cl{2}.TypeRec
        case 'inter'
            yl=mean(Cl{1},2);
        case 'domin'
            [r c]=size(Cl{1});
            ParIdx=round(1+(rand(r,1)*(Rho-1)));
            for l=1:r
                yl(l)=Cl{1}(l,ParIdx(l));
            end
    end
end