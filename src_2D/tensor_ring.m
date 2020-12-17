function [elem_list,err] = tensor_ring(elem_list,map,target_MPS, opts )
    
  if ~isfield(opts,'optim')
      opts.optim = 1:map.N;
  end
  
  %mapping if same element appears multiple times
  if ~isfield(opts,'get_elem_num')
       opts.get_elem_num = 1:map.N;
  end
  
  %'fsolve' or 'matrix_inv'
  if ~isfield(opts,'solve_type')
       error("provide solve type")
  end
  

 %select elements, including identical ones and subst num_elem for x;
    function elems = get_elem_list(num,x)
        elems = cell(1,map.N);
        if nargin<1
        
            for ii=1:map.N
                elems{ii} = elem_list{ opts.get_elem_num(ii) } ;
            end    
        else
            for ii=1:map.N
                elem_num = opts.get_elem_num(ii);
               
                if elem_num == opts.get_elem_num(num)
                    elems{ii} = x;
                else
                    elems{ii} = elem_list{ elem_num } ;
                end

            end
            
        end
    end
  
  optim_lookup = zeros(map.N,1);
  optim_lookup( opts.optim ) = 1;
  
  
  num_optim = size( opts.optim,2 );
  current_norms = zeros( map.N,1 );

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

    function [A,x0_shape]= get_A(num,x)

        if nargin< 2  
            temp_list = get_elem_list();
        else
            temp_list = get_elem_list(num,x);
        end
        
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
    
    function [F,G] = get_residual(num,duplicates,target,x)

        temp_list_1 = get_elem_list(num,x);

        A1  = ncon( temp_list_1,map.leg_list);
        
        F = reshape(  permute(A1,site_ordering_permute(map.N)), size(target)) -target;

        if nargout > 1 %calculate gradient
            
            size_x_red = size(x,[1,2,3,4,5,6]);
            d2 = size_x_red(1)^2;
            size_x_red(1) = d2;
            size_x_red(2)=[];

            grad_total_size = [size(F),size_x_red];
            Grad_total = zeros( grad_total_size );
            
            
            for ii = 1:size(duplicates,1)
                num = duplicates(ii);
                
                index_before = num-1;
                index_after = map.N-num;
                
                Grad_total = reshape(Grad_total, [d2^index_before, d2,d2^index_after, size_x_red ]    );
                
                [Ai,~]= get_A(num,x);
                
                Ai_res = reshape(Ai, [d2^index_before,1,d2^index_after,1, size_x_red(2:5)  ] );
                
                for iii =1:size_x_red(1)
                    %delta ij
                    Grad_total(:,iii,:,iii,:,:,:,:) = Grad_total(:,iii,:,iii,:,:,:,:)+ Ai_res(:,1,:,1,:,:,:,:);
                    
                end
                
                
            end
            
            
            G = reshape(Grad_total, numel(F), numel(x )  );
            
        end
        
    end

    function err = get_current_error()
      
        temp_list_1 = get_elem_list();
        A1  = ncon( temp_list_1,map.leg_list);
        a= reshape(  permute(A1,site_ordering_permute(map.N)), size(target_MPS)) -target_MPS; 
        
        
        err = sum(reshape( a,[],1).^2).^(0.5);
    end

    function  optimize_norms()
        norms_red = current_norms (current_norms>0 );
        
        [largest, argmax] = max(norms_red);
        [smallest, argmin] = min(norms_red);
        
        if largest/smallest >2 
            index_arr = [opts.optim(argmax), opts.optim(argmin)];
            value_arr = [largest,smallest];
            
            for iindex = 1:2
                l_index = index_arr(iindex);
                
                l=map.h_bond_l_lookup{l_index};
                r=map.h_bond_r_lookup{l_index};
                u=map.v_bond_u_lookup{l_index};
                d=map.v_bond_d_lookup{l_index};

                neighbours = [l_index];
                othernorms = [value_arr(iindex)];

                 if ~isempty(l)
                    pair = map.h_bonds{l};
                    other = pair(1);
                    if optim_lookup(other) == 1
                         neighbours=[ neighbours,other];
                         othernorms = [othernorms, current_norms(other)  ];
                    end
                end

                if ~isempty(u)
                    pair = map.v_bonds{u};
                    other = pair(1);
                    if optim_lookup(other) == 1
                         neighbours=[ neighbours,other];
                         othernorms = [othernorms, current_norms(other)  ];
                    end
                end

                if ~isempty(r)
                    pair = map.h_bonds{r};
                    other = pair(2);
                    if optim_lookup(other) == 1
                         neighbours=[ neighbours,other];
                         othernorms = [othernorms, current_norms(other)  ];
                    end
                end

                if ~isempty(d)
                    pair = map.v_bonds{d};
                    other = pair(2);
                    if optim_lookup(other) == 1
                         neighbours=[ neighbours,other];
                         othernorms = [othernorms, current_norms(other)  ];
                    end
                end

                ns = size(othernorms,2);
                avg_norm = prod(othernorms)^(1/ns);

                for m=1:ns
                   num =  neighbours(m);
                   n =   current_norms(num);
                   fact = avg_norm/n;
                   current_norms(num) = avg_norm;
                   elem_list{opts.get_elem_num(num)} = elem_list{opts.get_elem_num(num)}*fact;
                end
            end
            
            optimize_norms(); %call recursive
        end
    end

    err = Inf;
    for i=1:opts.maxit
 
        curr_ind = mod_index(i,num_optim);
        curr_elem = opts.optim(curr_ind);
        
        elem_index = opts.get_elem_num(curr_elem);
       
        if  elem_index ==  curr_elem

            switch opts.solve_type{curr_elem}
                case 'fsolve'

                    %for fsolve only unique subprobs should be solved


                        options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',true,'Algorithm','levenberg-marquardt' );

                        duplicates =  find( opts.get_elem_num  == elem_index  );

                        x = fsolve( @(x) get_residual(curr_elem, duplicates,target_MPS,x),  elem_list{ elem_index }, options );
                        normx = norm( reshape(x,[],1));

                        F = get_residual(curr_elem,duplicates,target_MPS,x);

                        %err = sum(reshape( F,[],1).^2).^(0.5);

                        current_norms(opts.get_elem_num(curr_elem)) = normx;
                        elem_list{opts.get_elem_num(curr_elem) } = x;



                case 'matrix_inv'

                    B = rotate_rhs(target_MPS,curr_elem);
                    [A,x0_shape] = get_A(curr_elem); 

                    x=lsqminnorm(A,B);

                    %F= A*x-B;

                    normx = norm(x);

                    x = permute(reshape(x, [x0_shape(3),x0_shape(4),x0_shape(5),x0_shape(6),x0_shape(1),x0_shape(2)]),[5,6,1,2,3,4]); %put phys dims back in front

                    %err = sum(reshape( F,[],1).^2).^(0.5);

                    current_norms(opts.get_elem_num(curr_elem)) = normx;
                    elem_list{opts.get_elem_num(curr_elem) } = x;


                otherwise 
                    error("unknown solve type")
            end
         end
 
        err = get_current_error();
        
        if curr_ind==num_optim
            if opts.print_level ==1
                
                fprintf("round %.4d err %.4e \n", (i)/num_optim ,err);
            end
           optimize_norms();
        end
        
        
        if err <opts.tol && i>num_optim
            break
            
        end
    
    end
    optimize_norms();

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

