function [elem_list,err] = tensor_ring(elem_list,map,target_MPS, opts )
    
  if ~isfield(opts,'optim')
      opts.optim = 1:map.N;
  end
     
  num_optim = size( opts.optim,2 );

    physical_dims = size(target_MPS);
    physical_dim = physical_dims(1);

    
    function y= mod_index(n,max)
        y= mod(n,max);
        if y == 0
            y=max;
        end
    end

    
    function C= rotate_rhs(B,num)
        B = reshape(B,physical_dims);
        
        perm = 1:map.N;
        perm(num)=[];
        perm = [perm,num];
        
        B = permute(B, perm );
        C=reshape(B,[],physical_dim);
    end

    function [A,x0_shape]= get_A(num)
        
        temp_list = elem_list;
        
        x0_shape = size(temp_list{num},[1,2,3,4,5,6]); %properly padded with ones
        
        temp_list(num) = []; %remove x0
        
        l=map.h_bond_l_lookup{num};
        r=map.h_bond_r_lookup{num};
        u=map.v_bond_u_lookup{num};
        d=map.v_bond_d_lookup{num};
        
        
        con_list_cpy = map.leg_list;
        
        ii = 0;
       
        
        
        if ~isempty(l)
            ii = ii+1; 
            pair = map.h_bonds{l};
            other = pair(1);
            con_list_cpy{other}(2+3) = -(map.external_legs+ii);
            
        end
        
        if ~isempty(u)
            ii = ii+1; 
            pair = map.v_bonds{u};
            other = pair(1);
            con_list_cpy{other}(2+4) = -(map.external_legs+ii);
        end
        
        if ~isempty(r)
            ii = ii+1; 
            pair = map.h_bonds{r};
            other = pair(2);
            con_list_cpy{other}(2+1) = -(map.external_legs+ii);
        end
        
        if ~isempty(d)
            ii = ii+1; 
            pair = map.v_bonds{d};
            other = pair(2);
            con_list_cpy{other}(2+2) = -(map.external_legs+ii);
        end
        
        num_legs_cpy = con_list_cpy{num};
        
        con_list_cpy( num) = [];
        

        final_order = -1:-1: -(map.external_legs+ii);
        seq = 1:map.internal_legs; %contraction sequence
        
        
        final_external_missing =sort( - num_legs_cpy(num_legs_cpy<0),'descend');
        final_internal_missing = sort( num_legs_cpy(num_legs_cpy>0),'descend');
        
        for l = 1:size(final_external_missing,2)
            final_order(  final_external_missing(l)  )= [];
        end
        
        for l = 1:size(final_internal_missing,2)
            seq(  final_internal_missing(l)  )= [];
        end
      
        A = ncon( temp_list,con_list_cpy,seq,final_order  );
      
        %put in site ordering
        perm_vector = [site_ordering_permute(map.N-1); ((2*map.N-1):size(size(A),2)).']; 

        
        A = reshape( permute(A, perm_vector)  , [],prod(x0_shape(3:end)));
    end
    
    err = Inf;
    for i=1:opts.maxit
 
        curr_ind = mod_index(i,num_optim);
        curr_elem = opts.optim(curr_ind);

        if opts.print_level ==1
           if curr_ind==1
                fprintf("round %.4d err %.4e \n", (i-1)/num_optim ,err);
           end
        end
        
        B = rotate_rhs(target_MPS,curr_elem);
        [A,x0_shape] = get_A(curr_elem); 
        
        x=lsqminnorm(A,B);
        elem_list{curr_elem} = permute(reshape(x, [x0_shape(3),x0_shape(4),x0_shape(5),x0_shape(6),x0_shape(1),x0_shape(2)]),[5,6,1,2,3,4]); %put phys dims back in front
        
        err = sum(reshape( A*x-B,[],1).^2).^(0.5);
        
        if err<opts.tol
            break
            
        end
        
       
    end

end


function p = site_ordering_permute(n)
% changes from |left i1 i2 ... j1 j2.. right> to |left i1 j1 i2 j2 ...
% right>
p = zeros(2*n,1);
%p(1)=1;
%p(2*n+2)=2*n+2;
for i = 1:n
    p(2*i-1)=i;
end
for i = 1:n
    p(2*i) = n+i ;
end
end

