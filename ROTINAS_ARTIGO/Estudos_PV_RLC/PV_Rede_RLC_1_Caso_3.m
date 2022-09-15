% set(0,'DefaultFigureWindowStyle','docked')
clc
close all

%% Vetores de Entrada
%G = [1000 800 1100];
%Tamb = [25 25 25];

G = 1000;
Tamb = 25;

%% Constantes:
Kb = 1.38064852e-23; 
q = 1.60217662e-19; 

%% Dados do Datasheet
Voc = 37.7;
Isc = 9.23; 
Vmp = 30.6;
Imp = 8.66;
Np = 1;
Ns = 60;
Pmax = 265; 
alpha = 0.053;
beta = -0.31;
NOCT = 45;
Isc_ref = Isc/Np;
Voc_ref = Voc/Ns;
Vmp_ref = Vmp/Ns;
Imp_ref = Isc/Np;

IB = Isc;
VB = Voc;
PB = Vmp*Imp;


%% Definição do numero maximo de iterações
Xmax = length(G);

%**************** DADOS DE SIMULAÇÃO *****************
dt = 0.1E-6;
tsim = 6E-3; % 3.5tempo de simulaçao em segundos
t_evento = tsim;

tmax = tsim; 
nmax = fix(tsim/dt); % Numéro de pontos avaliados na análise transitória
nmax = nmax*Xmax;

t = dt:dt:(tmax*Xmax); % tempo para verificar comportamento da resposta

% ****ELEMENTOS R,L e C DO CIRCUITO E CONDIÇOES INICIAS DOS ELEMENTOS****
Rlp = 0.25*(Vmp/Imp); Rcp = 1e-3;
L = 100e-6; C = 1e-6;
Rc = dt/(2*C); Rl = (2*L)/dt;
Rch_close = 1E-3;
Rch_open = 100E6;

% Alocação de variáveis
Il_hist = zeros(1,nmax);
Ic_hist = zeros(1,nmax);
V = zeros(5,nmax);
ikm_l = zeros(1,nmax);
ikm_c = zeros(1,nmax);
Im = zeros(1, nmax); 

Rd_hist_plt = zeros(1,nmax);
Id_hist_plt = zeros(1,nmax);

n = 2;

%% Inicio da rotina 
for u = 1:Xmax

%% Correção da Temperatura e Irradiancia
%% Cálculo dos parâmetros
fi = 1;   

%% Definição da irradiância referencia
Gref = 1000;

%% Calculo da temperatura da célula
T = 273 + Tamb(u);
T_ref = 273 + 25;
Tc = T+((NOCT-20)/0.8)*G(u)/Gref;
Tc_ref = T_ref+((NOCT-20)/0.8)*G(u)/Gref;

%% Aplicação das correções em Voc
varphi = 1+((beta)*(Tc-Tc_ref)/Ns);
phi = log(G(u)./Gref)/Ns;
Voc = Voc_ref*varphi+phi;

%% Aplicação das correções em Isc
psi = (1+alpha.*(Tc-Tc_ref)/Ns)*(G(u)./Gref);
Isc_new = Isc_ref*psi;
Iph = Isc_new;

%% Calculo da tensao terminca
Vt = (Kb*T)/q;

%% Calculo de Io e sua correcao
d = 7.02e-4;
Eg = 1.16-d.*(Tc.^2./(Tc-1108));
b = Eg*q./(fi*Vt);
tau = (Tc/Tc_ref).^(3/fi) * exp(b.*((Tc./Tc_ref) - 1));
Is_1 = Iph/(exp((Voc)/(fi*Vt)) - 1);
Is = Is_1*tau;
 
%% Calculo da resistencia RS
Xv = (Is/(fi*Vt))*exp(Voc_ref/(fi*Vt));
dVdI_Voc = - (0.43/Ns);
Rs = -dVdI_Voc-1/Xv;

while (n*dt <= tsim) 
    
    % Etapa Newton-Raphson (NR)
    Vk = V(1,n-1); % Tensão no diodo
    Vm = V(2,n-1); % Tensão terminal do módulo
    [Rd, Id_hist] = NR(Vk, Vm, Is, fi, Vt, Ns, Rs, Iph);
    
    % armazenamento para plotagem 
    Rd_hist_plt(n-1) = Rd;
    Id_hist_plt(n-1) = Id_hist;
    
    % Solução da rede discreta equivalente: RLC + Gerador Fotovoltaico 
    
    % Fontes de correntes históricas
    Ic_hist(n-1) = -V(3,n-1)/Rc - ikm_c(n-1);
    Il_hist(n-1) = V(5,n-1)/Rl  + ikm_l(n-1);
    
    % Matriz de condutância nodal 
    
%     0 - 1 ms: Chave aberta
%     1 - 2 ms: Chave fechada com Rcarga = 1,0 Ohms 
%     1,5 - 2,5ms: Redução de carga Rcarga = 10 Ohms
%     2,5 - 3,5ms: Redução de carga Rcarga = 0,1 Ohms

    if t(n-1) <= 1E-3
        Rch = Rch_open;
    elseif t(n-1)> 1E-3 && (t(n-1) <= 2E-3)
        Rch = Rch_close;
    elseif t(n-1)> 2E-3 && (t(n-1) <= 3E-3)
        Rlp = 0.5*(Vmp/Imp);
    elseif t(n-1)> 3E-3 && (t(n-1) <= 4E-3)
        Rlp = 0.75*(Vmp/Imp);
    elseif t(n-1)> 4E-3 && (t(n-1) <= 5E-3)
        Rlp = Vmp/Imp;     % Resistência que define o MPP na STC
    else
        Rlp = 1.25*Vmp/Imp;     % Resistência que define o MPP na STC
    end
    
    GG = [1/(Ns*Rd)+1/(Ns*Rs)  -1/(Ns*Rs) 0 0 0
        -1/(Ns*Rs) 1/(Ns*Rs)+1/(Rcp)+1/(Rch) -1/Rcp -1/Rch 0
        0 -1/Rcp 1/Rcp+1/Rc 0 0
        0 -1/Rch 0 1/Rlp+1/Rch -1/Rlp
        0 0 0 -1/Rlp 1/Rlp+1/Rl];  
    
    % Vetor fonte de corrente 
    I = [Iph-Id_hist; 0; -Ic_hist(n-1); 0; -Il_hist(n-1)];
  
    % Cálculo das tensões nodais
    V(:,n) = GG\I;
    
    % Atualização das correntes nos ramos
    ikm_l(n) = V(5,n)/Rl + Il_hist(n-1);
    ikm_c(n) = V(3,n)/Rc + Ic_hist(n-1);
    Im(n) = ikm_l(n) + ikm_c(n);
          
    n = n+1; 
end

tsim = tsim + t_evento;

end
t = t*1E3; %Tempo em ms


figure(1)
plot(t, V(2,:)/VB, '--k', t,ikm_l/IB,'-r', t, (V(2,:).*ikm_l)/PB,'g','LineWidth',1)
xlabel('Tempo (ms)')
ylabel('PU')
legend('V_M', 'I_L', 'P_{M}' )
grid on
box on
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7.5 4.5]); %Set the paper to have width 5 and height 5.
% print(gcf, 'C:\Users\bfons\Desktop\Figura1.pdf', '-dpdf','-r600'); close %Save figure

figure(2)
semilogy(t, Ns*Rd_hist_plt,'b', t,V(2,:)./ikm_l,'r')
ylabel('(Ohms)')
xlabel('Tempo (ms)')
legend('R_d','R_{pv}')
grid on
box on
ax = gca; 
ax.FontSize = 14; 
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7.5 4.5]); %Set the paper to have width 5 and height 5.
print(gcf, 'C:\Users\bfons\OneDrive\Desktop\Figura2.pdf', '-dpdf','-r600'); close %Save figure


% figure(1)
% %subplot(2,2,1)
% plot(t,V(2,:),'r',t,V(4,:),'b')
% xlabel('Tempo (ms)','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% ylabel('Tensao (V)','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% legend('V_2','V_4')
% %legend('V_2','V_4','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% %box off
% grid
% fig=1;
% saveas(fig,'V','png')
% 
% figure(2)
% %subplot(2,2,2)
% plot(t, Im, 'k', t,ikm_l,'r',t,ikm_c,'b')
% xlabel('Tempo (ms)','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% ylabel('Correntes na rede (A)','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% legend('I_M','I_L','I_C')
% %legend('I_M','I_L','I_C','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% %box off
% grid
% fig=2;
% saveas(fig,'I','png')
% 
% % figure(3)
% % %subplot(2,2,3)
% % plot(V(2,100:end), Im(100:end))
% % axis([0 41 0 (10)])
% % xlabel('Tensão Terminal (V)','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% % ylabel('Correntes de saída (A)','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% % set(gcf,'color','white')
% % %box off
% % hold off
% % grid
% % fig=3;
% % saveas(fig,'V_I','png')
% 
% figure(4)
% hold on
% semilogy(t, Ns*Rd_hist_plt)
% semilogy( t,V(2,:)./ikm_l)
% ylabel('(Ohms)')
% xlabel('Tempo (ms)')
% legend('N_s \times R_d','R_{pv}')
% hold off
% box off
% 
% figure(5)
% plot(t,Id_hist_plt,'r')
% ylabel('Id_{hist} (A)')
% xlabel('Tempo (ms)')
% box off
% 
% figure(6)
% plot( t,V(2,:).*ikm_l)
% ylabel('P_{pv} (W)')
% xlabel('Tempo (ms)')
% box off
% 
% % figure(7)
% % %subplot(2,2,1)
% % plot(t,V(2,:),'r',t,ikm_l,'--b')
% % xlabel('Tempo (ms)','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% % ylabel('Tensao (V)','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% % legend('V_M','I_L')
% % %legend('V_2','V_4','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% % %box off
% % grid
% % fig=1;
% % saveas(fig,'V','png')
% 
% % figure(8)
% % plot( V(2,:),ikm_l)
% % ylabel('I_{L} (A)')
% % xlabel('V_{M}(V)')
% % box off
% 
% 
% figure(8)
% plot(t, V(2,:)/VB, '--k', t,ikm_l/IB,'-r', t, (V(2,:).*ikm_l)/PB,'g','LineWidth',1)
% xlabel('Tempo (ms)')
% ylabel('PU')
% legend('V_M', 'I_L', 'P_{M}' )
% grid on
% box on
% 
% % figure(8)
% % hold on, grid on, box on
% % [hAx, hLine1, hLine2] = plotyy(t, V(2,:), t, ikm_l, 'plot','plot');
% % 
% % xlabel ('Tempo (ms)','FontSize',12)
% % % xlim(hAx(1),[6 22])
% % % xlim(hAx(2),[6 22])
% % 
% % ylabel(hAx(1),'V_M (V)','FontSize',12)    % left y-axis 
% % ylabel(hAx(2),'I_L (A)','FontSize',12)    % right y-axis
% % hLine1.LineStyle = '--';
% % hLine2.LineStyle = '-';
% % grid on
% % box on


