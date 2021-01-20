function bool = same_pattern(leg1, leg2)
    bool = 0;

    if size(leg1) == size(leg2)

        if leg1 == leg2
            bool = 1;
        end
    end
end

