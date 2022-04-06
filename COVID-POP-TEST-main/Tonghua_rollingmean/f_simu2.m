
clear 
load Data_TH.mat %data 
load chain2.mat
load results2.mat
xdata=data.xdata;


CIFcn = @(x,p)std(x,'omitnan')/sqrt(sum(~isnan(x))) * tinv(abs([0,1]-(1-p/100)/2),sum(~isnan(x))-1) + mean(x,'omitnan'); 

theta = mean(chain2,1);
theta

Nchain2 = length(chain2);

City = 'Tonghua';

K=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,theta(2)];
Nsimu = 100;
isample = randsample(Nchain2,Nsimu);

for i = 1:length(K)
    res_Theta = zeros(size(data.ydata,1),Nsimu);
    for j = 1:Nsimu
        theta = chain2(isample(j),:);
        theta(2) = K(i);
        res_Theta(:,j)=fsimu_deterministic(size(data.ydata,1),theta,xdata);
    end
    Mean = mean(res_Theta,2);
    UpCI = zeros(size(data.ydata,1),1);
    DownCI = zeros(size(data.ydata,1),1);
    for ii = 1:size(data.ydata,1)
        tmp = CIFcn(res_Theta(ii,:),95);
        UpCI(ii,1) = tmp(2);
        DownCI(ii,1) = tmp(1);
    end
    T = table(DownCI,Mean,UpCI, 'VariableNames', {'DownCI', 'Mean','UpCI'} );
    writetable(T, strcat('./simulation_result/',City,'_ContactTracingEffect_currentPopulationTesting.xlsx'),'Sheet',strcat('k=',sprintf('%.6f',K(i))),'WriteVariableNames',true)
end

%--No population level testing----

tau(1:size(data.ydata,1)) = 0.01;
for i = 1:length(K)
    res_Theta = zeros(size(data.ydata,1),Nsimu);
    for j = 1:Nsimu
        theta = chain2(isample(j),:);
        theta(2) = K(i);
        res_Theta(:,j)=fsimu_deterministic(size(data.ydata,1),theta,xdata,tau);
    end
    Mean = mean(res_Theta,2);
    UpCI = zeros(size(data.ydata,1),1);
    DownCI = zeros(size(data.ydata,1),1);
    for ii = 1:size(data.ydata,1)
        tmp = CIFcn(res_Theta(ii,:),95);
        UpCI(ii,1) = tmp(2);
        DownCI(ii,1) = tmp(1);
    end
    T = table(DownCI,Mean,UpCI, 'VariableNames', {'DownCI', 'Mean','UpCI'} );
    writetable(T, strcat('./simulation_result/',City,'_ContactTracingEffect_No_PopulationTesting.xlsx'),'Sheet',strcat('k=',sprintf('%.6f',K(i))),'WriteVariableNames',true)
end


%-------------------------------------------------------------------------------------

%-----Using the current settings from Tonghua, How many round of population testing

theta = mean(chain2,1);
 
ContactTracing = [0.2:0.02:0.5,theta(2),0]
Break = [2,1,0]
for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(2) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_deterministic(N,theta,xdata,tau);
        end
        writematrix(res_test,strcat('./simulation_result/',City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.4f',ContactTracing(j)),".csv"))
    end
end


theta = mean(chain2,1);
 
ContactTracing = [0.66:0.02:0.9];
Break = [2];
for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(2) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_deterministic(N,theta,xdata,tau);
        end
        writematrix(res_test,strcat('./simulation_result/',City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.4f',ContactTracing(j)),".csv"))
    end
end

theta = mean(chain2,1);
 
ContactTracing = [0.52:0.02:0.64];
Break = [2];
for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(2) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_deterministic(N,theta,xdata,tau);
        end
        writematrix(res_test,strcat('./simulation_result/',City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.4f',ContactTracing(j)),".csv"))
    end
end



theta = mean(chain2,1);
 
ContactTracing = [0.1:0.02:0.18];
Break = [2];
for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(2) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_deterministic(N,theta,xdata,tau);
        end
        writematrix(res_test,strcat('./simulation_result/',City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.4f',ContactTracing(j)),".csv"))
    end
end


%-------------------------------------------------------------------------------------------------------------



%----relative transimmision effect----------------------------------------

tmp = mean(chain2,1);
Ratio = [1:0.1:3.4];
Baseline = tmp(1);

ContactTracing = [0.35,0.4,0.45];
Break = [2];
theta = mean(chain2,1);


for r = 1:length(Ratio)
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(2) = ContactTracing(j);
            R0 = 3.2*Ratio(r);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_R_deterministic(N,theta,xdata,tau,R0);
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_Ratio",sprintf('%.4f',Ratio(r)),".csv"))
        end
    end
end

%------------------------------------------------------------------------------



%-----asympotomaic proportion effect--------------------------
theta = mean(chain2,1);
 
ContactTracing = [0.35,0.4,0.45];
Break = [2];
Prob = [0:0.02:0.9];
for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(2) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        for m =1:length(Prob)
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_P_deterministic(N,theta,xdata,tau,Prob(m));
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.4f',ContactTracing(j)),"_AymProb",sprintf('%.4f',Prob(m)),".csv"))
        end
        
    end
end


%----------------------------------------------------------------------------------------





%----P.1------------------------

theta = mean(chain2,1);
 
ContactTracing = [0.2:0.02:0.5,theta(2),0];
Break = [2];

theta = mean(chain2,1);
tmp = mean(chain2,1);
R0=1.7*3.2;


for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(2) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_R_deterministic(N,theta,xdata,tau,R0);
        end
        writematrix(res_test,strcat('./simulation_result/',City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.4f',ContactTracing(j)),"_P.1.csv"))
    end
end

%--------------------------------------



%----B.1.1.7------------------------

theta = mean(chain2,1);
 
ContactTracing = [0.2:0.02:0.5,theta(2),0];
Break = [2];

theta = mean(chain2,1);
tmp = mean(chain2,1);
R0=1.5*3.2;


for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(2) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_R_deterministic(N,theta,xdata,tau,R0);
        end
        writematrix(res_test,strcat('./simulation_result/',City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.4f',ContactTracing(j)),"_B.1.1.7.csv"))
    end
end

%--------------------------------------

%----Delta------------------------

theta = mean(chain2,1);
 
ContactTracing = [0.2:0.02:0.5,theta(2),0];
Break = [2];

theta = mean(chain2,1);
R0=2*3.2;


for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(2) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_R_deterministic(N,theta,xdata,tau,R0);
        end
        writematrix(res_test,strcat('./simulation_result/',City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.4f',ContactTracing(j)),"_Delta.csv"))
    end
end

%--------------------------------------

%----Omicron------------------------

theta = mean(chain2,1);
 
ContactTracing = [0.2:0.02:0.5,theta(2),0];
Break = [2];

theta = mean(chain2,1);
tmp = mean(chain2,1);
R0=10;


for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(2) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_R_deterministic(N,theta,xdata,tau,R0);
        end
        writematrix(res_test,strcat('./simulation_result/',City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.4f',ContactTracing(j)),"_Omicron.csv"))
    end
end

%--------------------------------------




%---middle size outbreak---------------

theta=mean(chain2,1);

theta

theta(2) = 0;
N=356;
tau = zeros(N,1);
MiddleSize=5000;
%%%% HS(i),HLa(i),HIa(i),HLs(i),HLp(i),HI(i),HR(i),i
Middle_value = fsimu_MiddleSize_deterministic_v2(N,theta,xdata,tau,MiddleSize);


Inital_Middle.HS = Middle_value.HS(size(Middle_value.HS,2));
Inital_Middle.HE = Middle_value.HE(size(Middle_value.HS,2));
Inital_Middle.HA = Middle_value.HA(size(Middle_value.HS,2));
Inital_Middle.HP = Middle_value.HP(size(Middle_value.HS,2));
Inital_Middle.HI = Middle_value.HI(size(Middle_value.HS,2));
Inital_Middle.HR = Middle_value.HR(size(Middle_value.HS,2));
Inital_Middle.HQ = Middle_value.HQ(size(Middle_value.HS,2));
Inital_Middle.HT = Middle_value.HT(size(Middle_value.HS,2));
Inital_Middle.HC = Middle_value.HC(size(Middle_value.HS,2));


theta = mean(chain2,1);
 
ContactTracing = [0.2:0.02:0.5,theta(2),0]
Break = [2]
for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(2) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_runFromMiddle_deterministic_v2(N,theta,xdata,tau,Inital_Middle);
        end
        writematrix(res_test,strcat('./simulation_result/',City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_Tracing",sprintf('%.4f',ContactTracing(j)),"_Middle.csv"))
    end
end


%----vaccine-----



%--------------Vaccine--------------------



%----P.1-----------------------------


theta = mean(chain2,1);
R0=1.7*3.2;

ContactTracing = [0.15,0.2,0.25,0.30,0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.77;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(2) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_v2(N,theta,xdata,tau,Vac_eff,R0);
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_Case_Test_break",int2str(Break(i)),"_Tracing",sprintf('%.2f',ContactTracing(j)),"_Vac",sprintf('%.2f',Vac),"_P.1.csv"))
        end
    end
end


theta=mean(chain2,1);
R0=1.7*3.2;
theta(1) = theta(1)*2;

p=0.284+0.25;

ContactTracing = [0.15,0.2,0.25,0.30,0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.77;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(2) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_asym_v2(N,theta,xdata,tau,Vac_eff,R0,p);
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_Case_Test_break",int2str(Break(i)),"_Tracing",sprintf('%.2f',ContactTracing(j)),"_Vac",sprintf('%.2f',Vac),"_Asym",sprintf('%.2f',p),"_P.1.csv"))
        end
    end
end




%--------------------------------------------------------------------------------

%----------------primary infection------------


theta=mean(chain2,1);
R0=3.2;
ContactTracing = [0.15,0.2,0.25,0.30,0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.88;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(2) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_v2(N,theta,xdata,tau,Vac_eff,R0);
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_Case_Test_break",int2str(Break(i)),"_Tracing",sprintf('%.2f',ContactTracing(j)),"_Vac",sprintf('%.2f',Vac),"_Primary.csv"))
        end
    end
end





theta=mean(chain2,1);
p=0.284+0.25;
R0=3.2;
ContactTracing = [0.15,0.2,0.25,0.30,0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.88;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(2) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_asym_v2(N,theta,xdata,tau,Vac_eff,R0,p);
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_Case_Test_break",int2str(Break(i)),"_Tracing",sprintf('%.2f',ContactTracing(j)),"_Vac",sprintf('%.2f',Vac),"_Asym",sprintf('%.2f',p),"_Primary.csv"))
        end
    end
end




%--------------------------------------------------------------------

%----reinfection with B.1.1.7-----------------------------


theta=mean(chain2,1);
R0=1.5*3.2;
ContactTracing = [0.15,0.2,0.25,0.30,0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.86;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(2) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_v2(N,theta,xdata,tau,Vac_eff,R0);
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_Case_Test_break",int2str(Break(i)),"_Tracing",sprintf('%.2f',ContactTracing(j)),"_Vac",sprintf('%.2f',Vac),"_B.1.1.7.csv"))
        end
    end
end



theta=mean(chain2,1);
R0=1.5*3.2;

p=0.284+0.25;

ContactTracing = [0.15,0.2,0.25,0.30,0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.86;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(2) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_asym_v2(N,theta,xdata,tau,Vac_eff,R0,p);
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_Case_Test_break",int2str(Break(i)),"_Tracing",sprintf('%.2f',ContactTracing(j)),"_Vac",sprintf('%.2f',Vac),"_Asym",sprintf('%.2f',p),"_B.1.1.7.csv"))
        end
    end
end

%--------------------------------------------------------------------------------
%----reinfection with Delta-----------------------------


theta=mean(chain2,1);
R0=2*3.2;

ContactTracing = [0.15,0.2,0.25,0.30,0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.65;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(2) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_v2(N,theta,xdata,tau,Vac_eff,R0);
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_Case_Test_break",int2str(Break(i)),"_Tracing",sprintf('%.2f',ContactTracing(j)),"_Vac",sprintf('%.2f',Vac),"_Delta.csv"))
        end
    end
end



theta=mean(chain2,1);
R0=2*3.2;

p=0.284+0.25;
ContactTracing = [0.15,0.2,0.25,0.30,0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.65;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(2) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_asym_v2(N,theta,xdata,tau,Vac_eff,R0,p);
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_Case_Test_break",int2str(Break(i)),"_Tracing",sprintf('%.2f',ContactTracing(j)),"_Vac",sprintf('%.2f',Vac),"_Asym",sprintf('%.2f',p),"_Delta.csv"))
        end
    end
end
%--------------------------------------------------------------------------------
%----reinfection with Omicron-----------------------------


theta=mean(chain2,1);
R0=10;

ContactTracing = [0.15,0.2,0.25,0.30,0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.44;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(2) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_v2(N,theta,xdata,tau,Vac_eff,R0);
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_Case_Test_break",int2str(Break(i)),"_Tracing",sprintf('%.2f',ContactTracing(j)),"_Vac",sprintf('%.2f',Vac),"_Omicron.csv"))
        end
    end
end



theta=mean(chain2,1);
R0=10;

p=0.284+0.25;
ContactTracing = [0.15,0.2,0.25,0.30,0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.44;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(2) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_asym_v2(N,theta,xdata,tau,Vac_eff,R0,p);
            end
            writematrix(res_test,strcat('./simulation_result/',City,"_Case_Test_break",int2str(Break(i)),"_Tracing",sprintf('%.2f',ContactTracing(j)),"_Vac",sprintf('%.2f',Vac),"_Asym",sprintf('%.2f',p),"_Omicron.csv"))
        end
    end
end


