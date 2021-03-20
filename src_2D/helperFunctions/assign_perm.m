function obj = assign_perm(obj, pat,rot_mask)


    if nargin<3
       rot_mask = [0,0,0,0]; 
    end
    
    orig_perms = 1:4;

    %orig = [a, a, 0, 0];
    all_perms = perms( orig_perms(~rot_mask)  );
    perm_pat = pat(all_perms);
    [~,ia,~] = unique(perm_pat,'rows');
    
    perm_vect = repmat(orig_perms, numel(ia),1) ;
    perm_vect(:,~rot_mask) = all_perms(ia,:);

    perm_val = repmat(pat, numel(ia),1) ;
    perm_val(:,~rot_mask) = perm_pat(ia,:);
    
    if sum(rot_mask)~=0 %also include  rotations
        pvect = perm_vect;
        pval = perm_val;
        
        for k =1:3
            perm_val = cat(1, perm_val, circshift(pval, k,2));
            perm_vect = cat(1, perm_vect, circshift(pvect, k,2));
        end
    end
    
    
    
    tensor0 = obj.PEPO_cell{pat(1) + 1, pat(2) + 1, pat(3) + 1, pat(4) + 1};
   
    tensor = sym_com(tensor0, pat,rot_mask);
   
    
    for i=1:size(perm_vect,1)
       vect = [1,2, perm_vect(i,:)+2];
       rpat = perm_val(i,:);
        
        
        eq_tensor = permute( tensor,   vect );

        %if isempty(obj.PEPO_cell{rpat(1) + 1, rpat(2) + 1, rpat(3) + 1, rpat(4) + 1})
            obj.PEPO_cell{rpat(1) + 1, rpat(2) + 1, rpat(3) + 1, rpat(4) + 1} = eq_tensor;
        %end

    
    end
end

function tens = sym_com(tens, pat,rot_mask)
    unique_nums = unique(pat(~rot_mask));
    orig_perm = 1:4;
    
    orig_tens = tens;
    
    for i=1:numel(unique_nums)
       n =  unique_nums(i);
       
     
       mask = find(pat == n);

       tens0 = tens;
       tens = zeros( size(tens) );


       [all_perms] =  perms( mask );
       np = size(all_perms,1);

       for ii=1:np
           this_perm = orig_perm;
           this_perm(mask) = all_perms(ii,:);

           this_perm = [1,2,this_perm+2];

           tens = tens + permute( tens0, this_perm)./np;
       end
       
       
       
    end
    
    %sym_err =svds( reshape( orig_tens-tens, [],1)  )
end

