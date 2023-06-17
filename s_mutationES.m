function Hatsl=s_mutationES(sl)
    nsl=zeros(size(sl{1},1),1);
    z=zeros(size(sl{1},1),1);
    switch sl{2}.TypeMut
        case 'isotr_1/5'
            %the one fifth rule
            if sl{2}.g>=sl{2}.G
                if sl{2}.Ps>(1/5)
                    s=sl{1}/sl{2}.a;                
                elseif sl{2}.Ps<(1/5)
                    s=sl{2}.a.*sl{1};
                else
                    s=sl{1};
                end
            else
                s=sl{1};
            end
            nsl=s;
            z=randn(size(nsl)).*nsl;
        case 'isotr_mul'
            %multiplicative (log-normal operator)
            s=sl{1};
            nsl=s.*exp(sl{2}.tao.*randn(size(s)));
            z=randn(size(nsl)).*nsl;
        case 'isotr_tpr'
            %two point rule (flip a coin)
            if rand<=0.5
                nsl=sl{1}.*(1+sl{2}.tao);
            else
                nsl=sl{1}./(1+sl{2}.tao);
            end
            z=randn(size(nsl)).*nsl;
        case 'nonis_mul'
            %non isotropic multiplicative (extended log-normal operator)
            s=sl{1};
            s=exp(sl{2}.tao0_ni*randn).*(sl{1}.*exp(sl{2}.tao_ni.*randn(size(s))));
            nsl=s;
            z=randn(size(nsl)).*nsl;
        case 'nonis_sco'
            %z =randn(2,15000)'*(1e150.*R);
            s=sl{1};
            s=exp(sl{2}.tao0_ni*randn).*(sl{1}.*exp(sl{2}.tao_ni.*randn(size(s))));            
            nsl=s;            
            z=randn(size(nsl)).*nsl;
            z(sl{2}.sc)=sl{2}.sC*z(sl{2}.sc);
    end    
    Hatsl={nsl,z};
end