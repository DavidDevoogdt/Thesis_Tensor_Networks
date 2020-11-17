function err = error_eigenvalue(mpo_1, mpo_2, type_02, M, d, opts)
%type= "array"
%type= "MPO"

%opts= "normed"


pa = inputParser;
addParameter(pa, 'ref', 0)
addParameter(pa, 'cyclic', 0);
parse(pa, opts)

a = mpo_1.contract_mpo(M-1, 1, pa.Results.cyclic);
prefact_1 = mpo_1.nf;


switch type_02
    case "array"
        b = mpo_2;
        prefact_2 = 1;
    case "MPO"
        b = mpo_2.contract_mpo(M-1, 1, pa.Results.cyclic);
        prefact_2 = mpo_2.nf;
    otherwise
        error("invalid type_01 for error_eigenvalue")
end


if pa.Results.ref == 2
    b = b * prefact_2^M;
    a = a * prefact_1^M;

else

    error("todo not implemented");
end

p = 2;

[~, S1, ~] = svds(a-b, 30);


sum1 = (sum(diag(S1).^p))^(1 / p);


[~, S2, ~] = svds(b, 30);


sum2 = (sum(diag(S2).^p))^(1 / p);

if pa.Results.ref ~= 0
    err = sum1 / sum2;
    return
else

    error("not impl");
end


end