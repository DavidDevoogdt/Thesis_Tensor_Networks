function err = error_eigenvalue(mpo_1, mpo_2,type_02,M,d,opts)
    %type= "array"
    %type= "MPO"
    
    %opts= "normed"
    
    
    p = inputParser;
    addParameter(p,'ref',0)
    addParameter(p,'cyclic',0);
    parse(p,opts)
    
    a = mpo_1.contract_mpo(M-1,1,p.Results.cyclic);
    prefact_1 = mpo_1.nf;
    

    switch type_02
        case "array"
            b = mpo_2;
            prefact_2=1;
        case "MPO"
            b = mpo_2.contract_mpo(M,0);
            prefact_2 = mpo_2.nf;
        otherwise
            error("invalid type_01 for error_eigenvalue")
    end
    

    
    if p.Results.ref == 2
        b=b*prefact_2^M;
        a =a*prefact_1^M;
        
    else 
        
       error("todo not implemented") ;
    end
        
    %try
    
    [~,S1,~] = svd(a-b);
    
    sum1 = 0;
    
    for i=1:size(S1,1)
        sum1= sum1+S1(i,i);   
    end
    
    
    [~,S2,~] = svd(b);
    
    sum2 = 0;
    
    for i=1:size(S2,1)
        sum2= sum2+S2(i,i);   
    end
    
    
    if p.Results.ref ~= 0
        err = sum1/sum2;
        return
    else
        
       error("not impl"); 
    end


end