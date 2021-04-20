function [Mapi,Mapm,Mapd] = compute_normalized_stats_Multiple_inputs(A,Bm,Var_A,Var_Bm,d,Var_d)

M = size(A,1);
Mapi = zeros(M);
Jm = size(Bm,3);
Mapm = zeros(M,M,Jm);
Mapd = zeros(M,1);
for m = 1:M
    for n = 1:M
        if m ~= n
            Mapi(m,n) = A(m,n)/sqrt(Var_A(m,n));
        end
        Mapm(m,n) = Bm(m,n)/sqrt(Var_Bm(m,n));
    end
end
for j = 1:Jm
    for m = 1:M
        for n = 1:M
            Mapm(m,n,j) = Bm(m,n,j)/sqrt(Var_Bm(m,n,j));
        end
    end
end

for m = 1:M
    Mapd(m) = d(m)/sqrt(Var_d(m));
end
    