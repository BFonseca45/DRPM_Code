function [Rd, Id_hist] = NR(Vk, Vm, Is, fi, Vt, Ns, Rs, Iph)

NRmax = 1000;
tol = 0.001;
for k = 1:NRmax
        
        Geq = (Is/(fi*Vt))*exp(Vk/(Ns*fi*Vt));
        Rd = 1/Geq;
        Id = Is*(exp((Vk/Ns)/(fi*Vt)) - 1);   % ATENTAR, NÃO É EQUAÇÃO (2) DO ARTIGO
%        Id_hist = Id - Vk/Rd; % Equação corrigida
        Id_hist = Id - Vk/(Ns*Rd); % Equação corrigida. Ver equção (37) do artigo
       
%        Vk_new = ( Vm/(Ns*Rs) + Iph - Id_hist) / ( 1/Rd + 1/(Ns*Rs));
        Vk_new = ( Vm/(Ns*Rs) + Iph - Id_hist) / ( 1/(Ns*Rd) + 1/(Ns*Rs)); % Equação corrigida. Ver equção (37) do artigo
        Geq = (Is/(fi*Vt))*exp(Vk_new/(Ns*fi*Vt));
        Rd = 1/Geq;
    
        f = Id + (Vk_new - Vk)/(Ns*Rs) - Iph;
        
        if (abs(f) < tol)
            break
        end
        
        Vk = Vk_new;
end