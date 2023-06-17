function BitCont=stop_criteriaES(BMu)
    BitCont=1;
    if any(BMu{3}==0)
        BitCont=0;  %nice! solution found
    elseif BMu{2}.g>=BMu{2}.MaxGens
        BitCont=0;
    end
end