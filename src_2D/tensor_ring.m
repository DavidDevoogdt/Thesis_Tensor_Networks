function [elem_list,err] = tensor_ring(target_MPS, dim,opts )
    
    %dims shouold be all the same

    physical_dims = size(target_MPS);
    physical_dim = physical_dims(1);
    
    num_elem = size(physical_dims,2);
    elem_list = cell(1,num_elem);
    con_list = cell(1,num_elem);
    
    sizes = [dim,physical_dim,dim];
    
    for i = 1:num_elem
        
        elem_list{i} = rand(sizes);
        
        con_list{i} = [ i,-i, i+1 ];
        
    end
    
    con_list{num_elem}(3)=1;
    
    function y= mod_index(n,max)
        y= mod(n,max);
        if y == 0
            y=max;
        end
    end

    
    function C= rotate_rhs(B,num)
        B = reshape(B,physical_dims);
        
        perm = 1:num_elem;
        perm(num)=[];
        perm = [perm,num];
        
        B = permute(B, perm );
        C=reshape(B,[],physical_dim);
    end

    function [A,x0]= get_A(num)
        final_order = -1:-1: -num_elem;
        temp_list = elem_list;
        x0 = temp_list{num};
        temp_list(num) = [];
        
        prev = mod_index(num-1,num_elem);
        next = mod_index(num+1,num_elem);
        
        con_list_cpy = con_list;
        con_list_cpy{prev}(3) = -(num_elem+1);
        con_list_cpy{next}(1) = -(num_elem+2);
        con_list_cpy( num) = [];
        
        
        final_order(num)= [];
        final_order = [final_order,-(num_elem+1),-(num_elem+2)];
        
        seq = 1:num_elem;
        a = num;
        b = next;
        l = max(a,b);
        k = min(a,b);
        
        seq(l)=[];
        seq(k)=[];
        
        A = ncon( temp_list,con_list_cpy,seq,final_order  );
        
        A = reshape(A, [],dim^2);
        x0 = reshape( permute(x0,[1,3,2]),[dim^2,physical_dim]);
    end
    

    function y = diff
       y = ncon(elem_list,con_list ) - target_MPS;
    end

    err = Inf;
    for i=1:opts.maxit
 
        curr_elem = mod_index(i,num_elem);

        if opts.print_level ==1
           if curr_elem==1
                fprintf("round %.4d err %.4e \n", (i-1)/num_elem,err);
           end
        end
        
        B = rotate_rhs(target_MPS,curr_elem);
        [A,~] = get_A(curr_elem); 
        
        x=lsqminnorm(A,B);
        elem_list{curr_elem} = permute(reshape(x, [dim,dim,physical_dim]),[1,3,2]);
        
        err = sum(reshape( A*x-B,[],1).^2).^(0.5);
        
        if err<opts.tol
            break
            
        end
        
       
    end

end


