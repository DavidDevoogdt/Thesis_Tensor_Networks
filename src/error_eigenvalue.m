function err = error_eigenvalue(mpo_1, mpo_2,type_02,M,d,opts)
    %type= "array"
    %type= "MPO"
    
    %opts= "normed"
    
    a = mpo_1.contract_mpo(M-1,1,0);
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
    

    p = inputParser;
    addParameter(p,'ref',0)
    parse(p,opts)
    
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

   
%     catch
%         warning('error = inf');
%         err = Inf;
%     end
    


%    tr_b = trace(b);
%    tr_a = trace(a);
% 
%    mag_2=tr_b*prefact_2;
%    mag_1=tr_a*prefact_1;
% 
%    
%    if mag_2 > mag_1 || p.Results.ref == 2
%        prefact = prefact_2*tr_b;
%        b = b/tr_b; %b normed to 1
%        a = a*(  prefact_1/ prefact);
%    else
%        prefact = prefact_2*tr_b;
%        a = a/tr_a; %b normed to 1
%        b = b*(  prefact_2/ prefact);
%    end
%                        
%     matrix= (a-b);
%     
%     try
%     	if p.Results.ref ~= 0
%             err = eigs(matrix,1);
%             return;
%         end
% 
%         err = prefact*eigs(matrix,1);
%     catch
%         warning('error = inf');
%         err = Inf;
%     end
    
    
   
   
end



function matrix = generate_cycle(mpo_list,d)
    M = size(mpo_list,2);
   
    leg_list = cell(1,M);
    for i = 1:M
       leg_list{i} = [i,-(i+1), -(i+M+1)  ,i+1  ] ;
    end

    leg_list{1}(1) = -1;
    leg_list{M}(4) = -(2*M+2);
    
    temp = ncon( mpo_list,leg_list   );
    
    matrix = reshape( temp, [ d^M,d^M ]);
    
end