function err = error_eigenvalue(mpo_1,type_01, mpo_2,type_02,M,d,opts)
    %type= "array"
    %type= "list" 
    %type= "O"
    
    %opts= "normed"
    
    switch type_01
        case "array"
            a = mpo_1;
            prefact_1=1;
        case "list"
            a= generate_cycle(mpo_1,d);
            prefact_1=1;
        case "O"
            N1 = mpo_1{1};
            O1 = mpo_1{2};
            
            list = generateList(O1,M,-1,0);
            a = generate_cycle(list,d);
            prefact_1=N1^M;
            
        otherwise
            error("invalid type_01 for error_eigenvalue")
    end
    
      
    switch type_02
        case "array"
            b = mpo_2;
            prefact_2=1;
        case "list"
            b= generate_cycle(mpo_1,d);
            prefact_2=1;
        case "O"
            N2 = mpo_2{1};
            O2 = mpo_2{2};
            list = generateList(O2,M,-1,0);
            b = generate_cycle(list,d);
            prefact_2=N2^M;
            
        otherwise
            error("invalid type_01 for error_eigenvalue")
    end
    

    p = inputParser;
    addParameter(p,'ref',0)
    parse(p,opts)
    

    
    
    
    
    
   tr_b = trace(b);
   tr_a = trace(a);

   mag_2=tr_b*prefact_2;
   mag_1=tr_a*prefact_1;

   
   if mag_2 > mag_1 || p.Results.ref == 2
       prefact = prefact_2*tr_b;
       b = b/tr_b; %b normed to 1
       a = a*(  prefact_1/ prefact);
   else
       prefact = prefact_2*tr_b;
       a = a/tr_a; %b normed to 1
       b = b*(  prefact_2/ prefact);
   end
                       
    matrix= (a-b);
    
    try
    	if p.Results.ref ~= 0
            err = eigs(matrix,1);
            return;
        end

        err = prefact*eigs(matrix,1);
    catch
        warning('error = inf');
        err = Inf;
    end
    
    
   
   
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