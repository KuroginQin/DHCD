function [base_init, coef_init] = NMF_rand_init(mat, K)
%Randam initialization of NMF
%mat: matrix to be factorized
%K: dimension of the hidden space
%base_init: initialized base matirx
%coef_init: initialized coefficient matrix

    %====================
    %Get the size of matrix to be factorized
    [N, M] = size(mat);
    %==========
    %Initialize the base matrix
    base_init = rand(N, K);
    norms = sqrt(sum(base_init.^2, 1));
    norms = max(norms, 1e-10);
    base_init = base_init./repmat(norms, N, 1);
    %==========
    %Initialize the coefficient matrix
    coef_init = rand(M, K);
    coef_init = coef_init/sum(sum(coef_init));
end

