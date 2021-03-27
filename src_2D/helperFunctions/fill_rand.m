function obj = fill_rand(obj, patterns, start_size,force)
    if nargin < 3
        start_size = 1e-1 / exp(obj.nf);
    end
    
    if nargin<4
       force = false; 
    end

    H = obj.virtual_level_sizes_horiz;
    V = obj.virtual_level_sizes_vert;
    d = obj.dim;

    %guess initial value
    for i = 1:size(patterns, 2)

        pat = patterns{i} + 1;

            if isempty(obj.PEPO_cell{pat(1), pat(2), pat(3), pat(4)}) || force == true
                pattern_s = [d, d, H(pat(1)), V(pat(2)), H(pat(3)), V(pat(4))];

                if obj.complex==true
                    m = rand(pattern_s)+1i*rand(pattern_s);
                else
                    m= rand(pattern_s);
                end

                obj.PEPO_cell{pat(1), pat(2), pat(3), pat(4)} = m * start_size;
            end
        
    end
end
