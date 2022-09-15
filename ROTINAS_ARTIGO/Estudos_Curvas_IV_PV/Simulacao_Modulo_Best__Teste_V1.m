set(0,'DefaultFigureWindowStyle','docked') 
clc, clear, close all

%% Vetores de Entrada
G = [1000 800 600 400];
Tamb = [25 25 25 25];

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
%Voc_ref = Voc/Ns;
Voc_ref = Voc;
%Vmp_ref = Vmp/Ns;
Vmp_ref = Vmp;
Imp_ref = Isc/Np;

%% Definição do numero maximo de iterações
Xmax = length(G);

%% Inicio da rotina 
for u = 1:Xmax

%% Fator de idealidade 
n = 0.9; 
%n = (Vmp_ref-Voc_ref)/(Vt*log(1-(Imp_ref/Isc_ref)));

%% Definição da irradiância referencia
Gref = 1000;

%% Calculo da temperatura da célula
T = 273 + Tamb(u);
T_ref = 273 + 25;
Tc = T+((NOCT-20)/0.8)*G(u)/Gref;
Tc_ref = T_ref + ((NOCT-20)/0.8)*(G(u)/Gref);

%% Aplicação das correções em Voc
% varphi = 1+((beta)*(Tc-Tc_ref)/Ns);
varphi = 1+((beta)*(Tc-Tc_ref));
% phi = log(G(u)./Gref)/Ns;
phi = log(G(u)./Gref);

%Voc = Voc_ref*varphi+phi;
Voc = (Voc_ref*varphi+phi);


%% Aplicação das correções em Isc
%psi = (1+alpha.*(Tc-Tc_ref)/Ns)*(G(u)./Gref);
psi = (1+alpha.*(Tc-Tc_ref))*(G(u)./Gref);
Isc_new = Isc_ref*psi;
Iph = Isc_new;

%% Calculo da tensao terminca
Vt = (Kb*T)/q;

%% Calculo de Io e sua correcao
d = 7.02e-4;
Eg = 1.16-d.*(Tc.^2./(Tc-1108));
%b = Eg*q./(n*Vt);    
b = Eg./(n*Vt);                                                             
tau = (Tc/Tc_ref).^(3/n) * exp(b.*((Tc./Tc_ref) - 1));
% Is_1 = Iph/(exp((Voc)/(n*Vt)) - 1);
Is_1 = Iph/(exp((Voc)/(Ns*n*Vt)) - 1);
Is = Is_1*tau;
 
%% Calculo da resistencia RS
%Xv = (Is/(n*Vt))*exp(Voc_ref/(n*Vt));
Xv = (Is/(n*Vt))*exp(Voc_ref/(Ns*n*Vt));

dVdI_Voc = - (0.43/Ns);
Rs = -dVdI_Voc-1/Xv;

% Variaveis auxiliares
Vsaida = []; Isaida = [];
iteration_number = 0;

for k = 0:20000
    
%    V = Voc*Ns*k/20000; % Resistencia simulando a carga
    V = Voc*k/20000; % Resistencia simulando a carga
    
    %Rv = 10e3; % Descomente para testar o circuito para uma carga especÃ­fica
    vd = 30; id = 0.001; % Chutes iniciais
    
    while(true)    
        % Determina os parametros Geq e Ieq:
        Geq = (Is/(n*Vt))*exp((vd/Ns)/(Vt*n));
        Ieq = id - Geq*(vd/Ns);
    
        % Encontra os novos valores de vd e id
        auxvd = (Iph - Ieq + V/(Ns*Rs))/(1/(Ns*Rs) + Geq/Ns); % Obtido diretamente atravÃ©s da Lei de Kirchhoff
        auxid = Is*(exp((auxvd/Ns)/(n*Vt)) - 1); % Obtido diretamente atravÃ©s da equaÃ§Ã£o de Shockley
    
        % Testa se o erro eh pequeno o suficiente
        if (abs(auxvd - vd) < 0.000001)
        
            % Atualiza os valores de id e vd
            vd = auxvd; id = auxid;
        
            % Sai do loop caso a convergencia seja alcancada
             break;
        end
    
        % Atualiza os valores de id e vd
        vd = auxvd; id = auxid;        
    
        % Contabiliza o numero de iteracoes
        iteration_number = iteration_number + 1;
    end    
    I2 = Iph - id;
    %V = (Rv/(Rv + Rs))*vd;
    Vsaida = [Vsaida; V]; 
    Isaida = [Isaida; I2];
end

V_Saida(:,u) = Vsaida;
I_Saida(:,u) = Isaida;


end

%% Carregamento dos Dados Experimentais

load ('Irrad.mat')

%% GRAFICOS 

% figure(1)
% plot(V_Saida(:,1),I_Saida(:,1),'Linewidth',1.5)
% set(gcf,'color','white')
% hold on
% plot(V_Saida(:,2),I_Saida(:,2),'Linewidth',1.5)
% plot(V_Saida(:,3),I_Saida(:,3),'Linewidth',1.5)
% plot(V_Saida(:,4),I_Saida(:,4),'Linewidth',1.5)
% axis([0 41 0 11])
% hold off
% grid
% xlabel('Tensão de Saída (V)','Interpreter','LaTex','FontSize',16)
% ylabel('Corrente Saída (A)','Interpreter','LaTex','FontSize',16)
% %title('Curva I x V do módulo fotovoltaico');
% %legend('Experimental','Simulação','Location','northwest','Orientation','vertical')
% text(2,9.7,'1000 $W/m^2$','Interpreter','LaTex','FontSize',12)
% text(2,7.84,'800 $W/m^2$','Interpreter','LaTex','FontSize',12)
% text(2,6,'600 $W/m^2$','Interpreter','LaTex','FontSize',12)
% text(2,4.2,'400 $W/m^2$','Interpreter','LaTex','FontSize',12)
% grid on
% set(gcf, 'PaperPosition', [0 0 7.5 4.5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7.5 4.5]); %Set the paper to have width 5 and height 5.
% fig=1;
% saveas(fig,'Curva','png')


figure(1)
plot(V_Saida(:,1),I_Saida(:,1),'Linewidth',1.5,'Color',[0,0,0])
axis([0 40 0 11])
hold on
plot(Xrad1000,Yrad1000,'o','Linewidth',1,'Color',[0,0,0])
set(gcf,'color','white')
plot(V_Saida(:,2),I_Saida(:,2),'Linewidth',1.5,'Color',[0,0,0])
plot(V_Saida(:,3),I_Saida(:,3),'Linewidth',1.5,'Color',[0,0,0])
plot(V_Saida(:,4),I_Saida(:,4),'Linewidth',1.5,'Color',[0,0,0])

plot(Xrad400,Yrad400,'o','Linewidth',1,'Color',[0,0,0])
plot(Xrad600,Yrad600,'o','Linewidth',1,'Color',[0,0,0])
plot(Xrad800,Yrad800,'o','Linewidth',1,'Color',[0,0,0])
hold off
grid
text(2,9.7,'1000 $W/m^2$','Interpreter','LaTex','FontSize',12)
text(2,7.84,'800 $W/m^2$','Interpreter','LaTex','FontSize',12)
text(2,6,'600 $W/m^2$','Interpreter','LaTex','FontSize',12)
text(2,4.2,'400 $W/m^2$','Interpreter','LaTex','FontSize',12)
xlabel('Tensão de Saída (V)','Interpreter','LaTex','FontSize',12)
ylabel('Corrente Saída (A)','Interpreter','LaTex','FontSize',12)
legend({'Modelo','Datasheet'},'Interpreter','LaTex','FontSize',14,'FontName','Times New Roman','Location','southwest','Orientation','vertical')
fig=1;
saveas(fig,'Curva_completa','png')

figure(2)
plot(V_Saida(:,1),V_Saida(:,1).*I_Saida(:,1),'Linewidth',1.5,'Color',[0,0,0])
axis([0 40 0 300])
hold on
plot(Xrad1000,Xrad1000.*Yrad1000,'o','Linewidth',1,'Color',[0,0,0])
set(gcf,'color','white')
plot(V_Saida(:,2),V_Saida(:,2).*I_Saida(:,2),'Linewidth',1.5,'Color',[0,0,0])
plot(V_Saida(:,3),V_Saida(:,3).*I_Saida(:,3),'Linewidth',1.5,'Color',[0,0,0])
plot(V_Saida(:,4),V_Saida(:,4).*I_Saida(:,4),'Linewidth',1.5,'Color',[0,0,0])

plot(Xrad400,Xrad400.*Yrad400,'o','Linewidth',1,'Color',[0,0,0])
plot(Xrad600,Xrad600.*Yrad600,'o','Linewidth',1,'Color',[0,0,0])
plot(Xrad800,Xrad800.*Yrad800,'o','Linewidth',1,'Color',[0,0,0])
hold off
grid
text(28,285,'1000 $W/m^2$','Interpreter','LaTex','FontSize',12)
text(28,230,'800 $W/m^2$','Interpreter','LaTex','FontSize',12)
text(28,180,'600 $W/m^2$','Interpreter','LaTex','FontSize',12)
text(28,125,'400 $W/m^2$','Interpreter','LaTex','FontSize',12)
xlabel('Tensão de Saída (V)','Interpreter','LaTex','FontSize',12)
ylabel('Potência de Saída (A)','Interpreter','LaTex','FontSize',12)
legend({'Modelo','Datasheet'},'Interpreter','LaTex','FontSize',14,'FontName','Times New Roman','Location','northwest','Orientation','vertical')
fig=2;
saveas(fig,'Curva_pot','png')