function zm = build_pathestmatrix(z)

T = length(z);
L = max(z);

zm = zeros(T,L)

for i=1:T
    zm(i,z(i)) = 1;
end
