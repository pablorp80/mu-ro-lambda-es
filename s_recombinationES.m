function sl=s_recombinationES(Cl,Rho)
    sl=zeros(size(Cl{2}.z,1),1);
    switch Cl{2}.TypeRec
        case 'inter'
            sl=mean(Cl{2}.z,2);
        case 'domin'
            [r c]=size(Cl{2}.z);
            ParIdx=round(1+(rand(r,1)*(Rho-1)));
            for l=1:r
                sl(l)=Cl{2}.z(l,ParIdx(l));
            end
    end
    sl={sl,Cl{2}};
end