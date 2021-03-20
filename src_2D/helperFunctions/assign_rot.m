function obj = assign_rot(obj, pat)
    %orig = [a, a, 0, 0];
    all_perms = perms([1,2,3,4]);
    perm_pat = pat(all_perms);
    [~,ia,~] = unique(perm_pat,'rows');
    perm_vect = all_perms(ia,:);
    perm_val = perm_pat(ia,:);
    
    
    tensor = obj.PEPO_cell{pat(1) + 1, pat(2) + 1, pat(3) + 1, pat(4) + 1};
    
    for i=1:size(perm_vect,1)
       vect = [1,2, perm_vect(i,:)+2];
       rpat = perm_val(i,:);
        
        
    eq_tensor = permute( tensor,   vect );

        if isempty(obj.PEPO_cell{rpat(1) + 1, rpat(2) + 1, rpat(3) + 1, rpat(4) + 1})
            obj.PEPO_cell{rpat(1) + 1, rpat(2) + 1, rpat(3) + 1, rpat(4) + 1} = eq_tensor;
        end

    
    end
end