set(0,'DefaultFigureWindowStyle','docked')
clc
close all
clear

%% Vetores de Entrada
Vin = 0: 0.1: 37.7; 

G = 1000*ones(length(Vin));
Tamb = 25*ones(length(Vin));
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
RB = Vmp/Imp;
%% Definição do numero maximo de iterações
Xmax = length(G);

%**************** DADOS DE SIMULAÇÃO *****************
dt = 50E-6;
tsim = 1.0E-3; % tempo de simulaçao em segundos
t_evento = tsim;

tmax = tsim; 
nmax = fix(tsim/dt); % Numéro de pontos avaliados na análise transitória
nmax = nmax*Xmax;

t = dt:dt:(tmax*Xmax); % tempo para verificar comportamento da resposta

% ****ELEMENTOS R,L e C DO CIRCUITO E CONDIÇOES INICIAS DOS ELEMENTOS****
Rlp = 1; Rcp = 1e-3;
L = 100e-6; C = 1e-6;
Rc = dt/(2*C); Rl = (2*L)/dt;
Rch_close = 1E-3;
Rch_open = 1E6;

% Alocação de variáveis
V = zeros(1,nmax);
Im = zeros(1, nmax); 
Rd_hist_plt = zeros(1,nmax);
Id_hist_plt = zeros(1,nmax);
Vm_plt = zeros(1,nmax);
I_Rd_plt = zeros(1,nmax);

n = 2;

%% Inicio da rotina 
for u = 1:Xmax

fi = 1;  % fator de idealidade 

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
    Vk = V(n-1); % Tensão no diodo
    Vm = Vin(u); % Tensão terminal do módulo
    [Rd, Id_hist] = NR(Vk, Vm, Is, fi, Vt, Ns, Rs, Iph);
    
    % armazenamento para plotagem 
    Rd_hist_plt(n) = Rd;
    Id_hist_plt(n) = Id_hist;
    Vm_plt(n) = Vm;
    I_Rd_plt(n) = Vk/(Ns*Rd);
    
    GG = 1/(Ns*Rs) + 1/(Ns*Rd);
    
    % Vetor fonte de corrente 
    I = Iph - Id_hist + Vm/(Ns*Rs);
  
    % Cálculo das tensão nodal no diodo 
    V(n) = GG\I;
    
    % Corrente de saída do módulo fotovoltaico 
     Im(n) = (V(n) - Vm)/(Ns*Rs);
          
    n = n+1; 
end

tsim = tsim + t_evento;

end

load ('Irrad.mat')
t = t*1E3; %Tempo em ms


figure(1)
hold on
plot(t,V(:),'r')
plot(t,Vm_plt,'b')
xlabel('Tempo (ms)')
ylabel('V_d (V)')
legend('V_d','V_M')
hold off
grid on
box on

figure(2)
%subplot(2,2,2)
plot(t, Im, 'k')
xlabel('Tempo (ms)')
ylabel('I_M(A)')
%legend('I_M','I_L','I_C')
%legend('I_M','I_L','I_C','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
%box off
grid on
box on

% figure(3)
% %subplot(2,2,3)
% plot(V(2,100:end), Im(100:end))
% axis([0 41 0 (10)])
% xlabel('Tensão Terminal (V)','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% ylabel('Correntes de saída (A)','Interpreter','LaTex','FontSize',13,'FontName','Times New Roman')
% set(gcf,'color','white')
% %box off
% hold off
% grid
% fig=3;
% saveas(fig,'V_I','png')

% Rbase = max(Ns*Rd_hist_plt);
% Pbase = max(Im.*Vm_plt);

figure(4)
%hold on
%semilogy(t, Rd_hist_plt)
semilogy(Vm_plt/VB, Ns*Rd_hist_plt/RB, 'k', ...
         Vm_plt/VB, (Im.*Vm_plt)/PB, '--b', ...
         Vm_plt/VB, (Vm_plt./Im)/RB, ':r','LineWidth',2)
ylabel('Valor em PU')
legend({'R_d','P_{pv}','R_{pv}'},'Location','southwest','Orientation','vertical')
%legend('R_d','P_{pv}','R_{pv}')
%xlabel('Tempo (ms)')
xlabel('V_M (PU)')
xlim([0 1.1])
grid on
box on
ax = gca; 
ax.FontSize = 14; 
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7.5 4.5]); %Set the paper to have width 5 and height 5.
print(gcf, 'C:\Users\bfons\OneDrive\Desktop\Figura_XX.pdf', '-dpdf','-r600'); close %Save figure

figure(5)
% hold on
% plot(t,Id_hist_plt,'r')
% plot(t,I_Rd_plt,'--b')
% plot(t, I_Rd_plt+Id_hist_plt, '.-g' )
% ylabel('Id_{hist} (A)')
% xlabel('Tempo (ms)')
% legend('Id_{hist}','I_{Ns \times Rd}','I_d')
% grid on
% box on

hold on
plot(Vm_plt,Id_hist_plt,'r')
plot(Vm_plt,I_Rd_plt,'--b')
plot(Vm_plt, I_Rd_plt+Id_hist_plt, '.-g' )
ylabel('Corrente(A)')
xlabel('V_M (V)')
legend('Id_{hist}','I_{ Rd}','I_d','Location','northwest')
grid on
box on
ax = gca; 
ax.FontSize = 14; 
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7.5 4.5]); %Set the paper to have width 5 and height 5.
print(gcf, 'C:\Users\bfons\OneDrive\Desktop\Figura_YY.pdf', '-dpdf','-r600'); close %Save figure


figure(6)
hold on
plot(Vm_plt(20:end), Im(20:end),'r')
plot(Xrad1000,Yrad1000,'o','Linewidth',1,'Color',[0,0,0])
% plot(Xrad800,Yrad800,'o','Linewidth',1,'Color',[0,0,0])
% plot(Xrad600,Yrad600,'o','Linewidth',1,'Color',[0,0,0])
% plot(Xrad400,Yrad400,'o','Linewidth',1,'Color',[0,0,0])
legend('Modelo','Datasheet')
ylabel('I_M (A)')
xlabel('V_M (V)')
grid on
box on
% set(gcf, 'PaperPosition', [0 0 7.5 4.5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7.5 4.5]); %Set the paper to have width 5 and height 5.
% print(gcf, 'C:\Users\bfons\Desktop\Curva_completa.pdf', '-dpdf','-r600'); close %Save figure


figure(7)
hold on
plot(Vm_plt, Vm_plt.*Im,'r')
plot(Xrad1000,Xrad1000.*Yrad1000,'o','Linewidth',1,'Color',[0,0,0])
% plot(Xrad800,Xrad800.*Yrad800,'o','Linewidth',1,'Color',[0,0,0])
% plot(Xrad600,Xrad600.*Yrad600,'o','Linewidth',1,'Color',[0,0,0])
% plot(Xrad400,Xrad400.*Yrad400,'o','Linewidth',1,'Color',[0,0,0])
legend('Modelo','Datasheet','Location','northwest')
ylabel('P_M (W)')
xlabel('V_M (V)')
grid on
box on
% set(gcf, 'PaperPosition', [0 0 7.5 4.5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7.5 4.5]); %Set the paper to have width 5 and height 5.
% print(gcf, 'C:\Users\bfons\Desktop\Curva_Pot.pdf', '-dpdf','-r600'); close %Save figure