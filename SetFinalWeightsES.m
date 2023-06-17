function nnet=SetFinalWeightsES(Hatyl,net)
    %update input weights
    Layer1Dim=net.layers{1}.dimensions;
    IWs=zeros(size(net.IW{1}));
    for input=1:net.inputs{1}.size        
        IWs(:,input)=Hatyl((Layer1Dim*(input-1))+1:(Layer1Dim*input));
    end
    net.IW{1}=IWs;
    lastIdx=Layer1Dim*input;
    %update internal weights
    for layer=1:net.numlayers-1
        [r c]=size(net.LW{layer+1,layer});
        net.LW{layer+1,layer}=reshape(Hatyl(lastIdx+1:(r*c)+lastIdx),r,c);
        lastIdx=(r*c)+lastIdx;        
    end
    %update biases
    for layer=1:net.numlayers
        net.b{layer}=Hatyl(lastIdx+1:net.layers{layer}.dimensions+lastIdx);
        lastIdx=net.layers{layer}.dimensions+lastIdx;
    end
    nnet=net;

end