function Cl=reproductionES(BMu,Rho)
    sel=round(1+rand(1,Rho)*(BMu{2}.Mu-1));    
    WCandidate=BMu{1}(:,sel);
    F=BMu{3}(sel);
    s=BMu{2};
    s.z=s.z(:,sel);
    Cl={WCandidate, s, F};
end