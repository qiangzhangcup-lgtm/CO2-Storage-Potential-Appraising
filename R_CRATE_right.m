clc
clear all

t_w=xlsread('CO2PTA.xlsx','sheet2','A2:A8920'); 
qg_ECL=xlsread('CO2PTA.xlsx','sheet2','E2:E8920');  %注入
qw_ECL=xlsread('CO2PTA.xlsx','sheet2','E2:E8920'); 
Wp=xlsread('CO2PTA.xlsx','sheet2','F2:F8920'); % bbl
pwin_ECL_w=xlsread('CO2PTA.xlsx','sheet2','D2:D8920'); % ECL pave
pbar_ecl_w=xlsread('CO2PTA.xlsx','sheet2','B2:B8920'); % ECL pave 
Sw_ECL_w=xlsread('CO2PTA.xlsx','sheet2','C2:C8920'); 
pbar_MBE(i)=xlsread('CO2PTA.xlsx','sheet2','G2:G8920'); 
Sw_MBE(i)=xlsread('CO2PTA.xlsx','sheet2','H2:H8920'); 
Sw_ECL_w=Sw_ECL_w./100; 
% 1 means gas, 2 means water
t=t_w;
aa=0;
    Control=2    %2为定压，1为定产  
    Pi=pbar_ecl_w(1); % (in put) psi initial condition 
    Cw=1E-7;%3.17E-5; 
    Cp=4.7E-6;%4.7E-6;%1e-8
    Cgi=(Pi)^-1;   % for7 g10as model
    k=100;  
    h=50; %1ft=0.3048m
    r=5000;
    rw_set=5;
    fai=0.1; % fracture porosity
    Bwi=1.02; %^1.03;
    ga=8E-6;
    uw=1.25;  
    ug(1)=(0.0254.*log(Pi)-0.1266).*1.001;%0.0357.*lexp(0.0002.*Pi);
   
    Bgi=48.45338.*5.615./Pi;
Q1=7000; %max(Q,Q_g)
tspt_MBEa=zeros(Q1,1);
tspt_MBE=zeros(Q1,1);        
  %%     物质平衡法求平均裂缝压力
    Swi= Sw_ECL_w(1);
    Sw_MBE(1)=Swi;  
    Sw_MBEa(1)=Swi;
    Sg_MBE(1)=1-Swi;
    dSw(1)=0;
    pbar_MBE(1)=Pi;
    pbar_MBEa(1)=pbar_MBE(1);
    Cgefmi=(1-Swi).*(Cgi+Cp)-dSw(1);
    Vfi=3.14*r^2*h*fai; % ft3  
    Krw_endpoint=1;
    Krg_endpoint=1;
    Swa=0.2;
    Sgc=0;
    krwn=3;%相渗指数
    krgn=1.5;%相渗指数
    Krg(1)=Krg_endpoint.*((1-Swi-Sgc)/(1-Swa-Sgc)).^krgn; 
    Krw(1)=Krw_endpoint.*((Swi-Swa)/(1-Swa-Sgc)).^krwn;       
   %%%%%%%%%%%%%%%%%%%平均压力
              
   %% 水相 特征曲线 
      for i=2:Q1 
       Swi=Sw_ECL_w(1);   
       ug(i)=(0.0254.*log(pbar_MBE(i))-0.1266).*1.001;%0.0357.*exp(0.0002.*pbar_MBE(i));
       Cg(i)=1./pbar_MBE(i);
       Krg(i)=Krg_endpoint.*((1-Sw_MBE(i)-Sgc)/(1-Swa-Sgc)).^krgn;%0.16^2.*(1-Sw_MBE(i));%3.67.*Sw_MBE(i).^2-5.3.*Sw_MBE(i)+1.87; 
       Krw(i)=Krw_endpoint.*((Sw_MBE(i)-Swa)/(1-Swa-Sgc)).^krwn; 
       Bg(i)=48.45338.*5.615./pbar_MBE(i);     
       M(i)=Krw(i)./uw./(Krg(i)./ug(i));  
     end
   %% 定产拟压力 
   for i=2:Q1    
       N=1:1000; 
   if qg_ECL(100)==qg_ECL(102);%定产 
       
%          co_P_Sw=polyfit(pbar_MBE_BL(3000:Q1),Sw_BL_Ava(3000:Q1),1);
                                     % for gas model
        pdiscrete=Pi+(pwin_ECL_w(i)-Pi)/(length(N)).*N; % discreticize pressure from pwf to pi
        Sw_discrete=-0.79.*log(pdiscrete)+6.625;% co_P_Sw(1)*pdiscrete+co_P_Sw(2);
%         Sw_discrete=co_P_Sw(1).*pdiscrete+co_P_Sw(2);
        Sw_discrete(Sw_discrete > 0.8) = 0.798;
        %         Sw_discrete=(-0.000006.*pdiscrete+75.432)/100;%-3E-06.*pdiscrete+0.7042;%=-0.00065.*pdiscrete+1.642;%-2E-5.*pdiscrete+0.4234;%=-2E-13.*pdiscretexp3+3E-10.*pdiscretexp2+3E-07.*pdiscrete+0.3992;%-6E-15.*pdiscretexp4 + 6E-11.*pdiscretexp3 - 2E-07.*pdiscretexp2 + 0.0003.*pdiscrete + 0.2273;%-2E-05.*pdiscrete+0.4231(平均压力);
        Bg_d=48.45338.*5.615./pdiscrete;
        ug_d=(0.0254.*log(pdiscrete)-0.1266).*1.001;%0.0357.*exp(0.0002.*pdiscrete);
        Krg_d= Krg_endpoint.*((1-Sw_discrete-Sgc)/(1-Swa-Sgc)).^krgn;%0.16^2.*(1-Sw_MBE(i));%3.67.*Sw_MBE(i).^2-5.3.*Sw_MBE(i)+1.87; 
        innerppi=Krg_d./ug_d./Bg_d.*exp(-(ga).*(Pi-pdiscrete));
        ppoi(i)=ug(1).*Bgi.*trapz(pdiscrete,innerppi);      
        
%         co_P_Sw_2=polyfit(pbar_MBE_BL(1000:Q1),Sw_MBE(1000:Q1),1);
                                     % for gas model
        pdiscrete_FMB=pwin_ECL_w(i)+(pwin_ECL_w(i)-pbar_MBE(i))/(length(N)).*N;
%        Sw_discrete_2=co_P_Sw_2(1)*pdiscrete_FMB+co_P_Sw_2(2);
%        Sw_discrete_2=0.66;%
%         Sw_discrete_2=co_P_Sw_2(1).*pdiscrete_FMB.^3+co_P_Sw_2(2).*pdiscrete_FMB.^2+co_P_Sw_2(3).*pdiscrete_FMB+co_P_Sw_2(4);
         Sw_discrete_2=-0.79.*log(pdiscrete_FMB)+6.525;;
        Bg_FMB=48.45338.*5.615./pdiscrete_FMB;
        ug_FMB=(0.0254.*log(pdiscrete_FMB)-0.1266).*1.001;%0.0357.*exp(0.0002.*pdiscrete_FMB);
        Krg_FMB=Krg_endpoint.*((1-Sw_discrete_2-Sgc)/(1-Swa-Sgc)).^krgn;%0.16^2.*(1-Sw_MBE(i));%3.67.*Sw_MBE(i).^2-5.3.*Sw_MBE(i)+1.87; 
        innerppi_FMB=Krg_FMB./ug_FMB./Bg_FMB.*exp(-(ga)*(pbar_MBE(1)-pdiscrete_FMB));
        ppoi_MBE(i)=ug(1).*Bgi.*trapz(pdiscrete_FMB,innerppi_FMB);
   else
        ug2=(0.0254.*log(pbar_MBE)-0.1266).*1.001;%0.0357.*exp(0.0002.*pbar_MBE');
        Bg2=48.45338.*5.615./pbar_MBE;
        Krg2= Krg_endpoint.*((1-Sw_MBE-Sgc)./(1-Swa-Sgc)).^krgn;
        Krw2=Krw_endpoint.*((Sw_MBE-Swa)./(1-Swa-Sgc)).^krwn;
        Cg2=1./pbar_MBE;
        Ce2=(Cg2.*(1-Sw_MBE)+Cw.*Sw_MBE);
        innerppi_g=Krg2./ug2./Bg2.*exp(-(ga).*(Pi-pbar_MBE)); % Bg=17.553343*p^(-1)
        ppoi=trapz(pbar_MBE,innerppi_g).*ug(1)*Bgi+zeros(Q1,1); %和张老师模型保持一致

        pdiscrete_FMB=pbar_MBE(i)+(pwin_ECL_w(i)-pbar_MBE(i))/(length(N)).*N;
        %         Sw_discrete_2=(-6E-04.*pdiscrete_FMB +75.857)/100; 
        co_P_Sw=polyfit(pbar_MBE(1:Q1),Sw_MBE(1:Q1),1);
        intercept_P_Sw= co_P_Sw(2);  
        Sw_discrete_2=co_P_Sw(1)*pdiscrete_FMB+intercept_P_Sw;
        Bg_FMB=48.45338.*5.615./pdiscrete_FMB;
        ug_FMB=(0.0254.*log(pdiscrete_FMB)-0.1266).*1.001;%0.0357.*exp(0.0002.*pdiscrete_FMB);
        Krg_FMB= Krg_endpoint.*((1-Sw_discrete_2-Sgc)/(1-Swa-Sgc)).^krgn;
        innerppi_FMB=Krg_FMB./ug_FMB./Bg_FMB.*exp(-(ga)*(pbar_MBE(1)-pdiscrete_FMB));
        ppoi_MBE(i)=ug(1).*Bgi.*trapz(pdiscrete_FMB,innerppi_FMB);
   end 
        RNP(i)=ppoi(i)./qg_ECL(i);
        PNP(i)=qg_ECL(i)./ppoi(i);
        g_FMB(i)=pi.*r^2.*fai.*h.*(ppoi(i)-ppoi_MBE(i))./ppoi(i);
   end
    for i=2:Q1-1 
        dSw_MBE1(i)=-(Sw_MBE(i)-Sw_MBE(i-1))./7./(pwin_ECL_w(i)-pwin_ECL_w(i-1));
        dSw_MBE1(1)=-(Sw_MBE(2)-Sw_MBE(1))./7./(pwin_ECL_w(2)-pwin_ECL_w(1));
        Ce(i)=(1-Sw_MBE(i)).*(Cp+Cg(i))+dSw_MBE1(i);
        Cei=(1-Swi).*(Cp+Cgi)+dSw_MBE1(1);
        
        innerint_MBE(i)=Krg(i)./ug(i).*exp(-(ga-Cp).*(Pi-pwin_ECL_w(i)))./Ce(i);      
        tp_MBE(i,1)=ug(1).*Cei.*trapz(t(1:i),innerint_MBE(1:i));% 
    if qg_ECL(100)==qg_ECL(102);
        tspt_MBE(i,1)=tp_MBE(i,1);
    else
        for j1=2:i
        tspt_MBE(i,1)=tspt_MBE(i,1)+(qg_ECL(j1)-qg_ECL(j1-1))*(tp_MBE(i)-tp_MBE(j1-1))./(qg_ECL(i));
        end
    end
    end
for i=2:Q1-1  
    DRNPw(i,1)=(RNP(i+1)- RNP(i-1))./(log(tspt_MBE(i+1))-log(tspt_MBE(i-1)));         
end 
   
%        tp_MBE=tspt_MBE;

      hold on  
      figure (3)    % disgnostic plot to see half slope
      loglog(tspt_MBE(2:Q1-1),DRNPw(2:Q1-1),'.b')
      xlabel('等效拟时间');
      ylabel('DRNPw');
      legend('Tr_shutin vs  DRNP')%,'gas tp ANA MBE   vs  RNP gas');
      title('焖井阶段水诊断曲线');        
%% 拟合直线 做判断曲线
     RL=2000; % for variable q 
     Q=6990;
     Flow=2
     tspt_MBE=tspt_MBE';%tspt_MBE;
    %%%for IALF
   if Flow==0
        trp=tspt_MBE.^0.5;  
        co_IALF=polyfit(trp(RL:Q),RNP(RL:Q),1);        % for gas model
        slope_IALF=co_IALF(1);                                  % for gas model
        intercept_IALF=co_IALF(2);                                 % for gas model
    for i=2:Q1
        regress_IALF(i)=slope_IALF.*trp(i)+ intercept_IALF;    % for gas model  拟合判断直线？？？
    end
        k_IALF=(14.18*Bgi*ug(1).^0.5/wf/h/fai^0.5/Cei^0.5/slope_IALF)^2
            hold on
              figure (6) % specialty plot
              plot(trp(2:Q1),RNP(2:Q1),'.r',trp(2:Q1),regress_IALF(2:Q1),'-b')
              xlabel('拟时间^0.5');
              ylabel('RNPg');
              legend('tp ANA MBE  vs  RNP','straight line');
              title('线性流气相特征曲线');
   else 
        co_BDF=polyfit(tspt_MBE(RL:Q),RNP(RL:Q),1);        % for gas model
        co_FMB=polyfit(g_FMB(RL:Q),PNP(RL:Q),1);
        slope_BDF=co_BDF(1);                                  % for gas model
        intercept_BDF=co_BDF(2);                              % for gas model
        slope_FMB=co_FMB(1);                                  % for gas model
        intercept_FMB=co_FMB(2);                              % for gas model 
        for i=2:Q1
        regress_BDF(i)=slope_BDF.* tspt_MBE(i)+intercept_BDF;    % for gas model  拟合判断直线？？？
        regress_FMB(i)=slope_FMB.* g_FMB(i)+intercept_FMB;  
        end    
        rm_BDF=sqrt(Bgi/fai/h/Cei/3.14/slope_BDF) % for gas model
        km_BDF=157.62*Bgi*ug(1)/h/(2*3.14)*(log(rm_BDF/rw_set)-3/4)/intercept_BDF
        
       rm_FMB=sqrt(5.615*Bgi.*(-intercept_FMB./slope_FMB)/fai/h./3.14)
        km_FMB=157.62*intercept_FMB*Bgi*ug(1)*(log(rm_FMB/rw_set)-3/4)/h/(2*pi)
        
        V_CO2=pi.*rm_FMB^2.*fai.*h.*0.6.*(ppoi(3000)-ppoi_MBE(3000))/5.615/1000;
        
        hold on
        figure (7) % specialty plot
        plot((tspt_MBE(2:Q1)),RNP(2:Q1),'.r',(tspt_MBE(2:Q1)),regress_BDF(2:Q1),'-b')
        xlabel('拟时间');
        %%
        % *BOLD TEXT*
        ylabel('RNPg');
        legend('tp ANA MBE  vs  RNP','straight line');
        title('边界流气相特征曲线');
        figure (8) % specialty plot
        plot((g_FMB(2:Q1)),PNP(2:Q1),'.r',(g_FMB(2:Q1)),regress_FMB(2:Q1),'-b')
        xlabel('物质平衡时间');
        ylabel('PNPg');
        legend('tp ANA MBE  vs  RNP','straight line');
        title('边界流气相特征曲线');
        
   end 
              





