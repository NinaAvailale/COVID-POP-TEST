function pij = ContactTracing(A, P, I, tau, currentTime,rE,rA,rP,rI,pa,t0,tL,Sens1,Sens2)

if currentTime<=tL
	tL = currentTime;
	tL = tL-1;
end

%--the transition rate from compartment K-----
%--Note the tau is changing in reality. To make the calculation easy, the average tau is used here. 
tau = mean(tau);
rE_star = rE+tau*Sens1;
rA_star = rA+tau*Sens2;
rP_star = rP+tau*Sens2;
rI_star = rI+tau*Sens2;

%--fraction of individuals who leave J that reack K---
qEA = pa*rE/rE_star;
qEP = (1-pa)*rE/rE_star;
qPI = rP/rP_star;
qIR = rI/rI_star;
qAR = rA/rA_star;

%--transition probability function---------------
function a = P_AE(qEA,rE_star,rA_star,Time)
	a = qEA*rE_star*(exp(-rE_star*Time)-exp(-rA_star*Time))/(rA_star-rE_star);
end

function a = P_PE(qEP,rE_star,rP_star,Time)
	a = qEP*rE_star*(exp(-rE_star*Time)-exp(-rP_star*Time))/(rP_star-rE_star);
end

function a = P_IP(qPI,rI_star,rP_star,Time)
	a = qPI*rP_star*(exp(-rP_star*Time)-exp(-rI_star*Time))/(rI_star-rP_star);
end

function a = P_IE(qPI,qEP,rI_star,rE_star,rP_star,Time)
	a = qEP*qPI*rE_star*rP_star*((exp(-rE_star*Time)-exp(-rP_star*Time))/(rE_star-rP_star)-(exp(-rE_star*Time)-exp(-rI_star*Time))/(rE_star-rI_star))/(rP_star-rI_star);
end

function a = P_EE(rE_star,Time)
	a = exp(-rE_star*Time);
end

function a = P_AA(rA_star,Time)
	a = exp(-rA_star*Time);
end

function a = P_PP(rP_star,Time)
	a = exp(-rP_star*Time);
end

function a = P_II(rI_star,Time)
	a = exp(-rI_star*Time);
end

%----contact tracing removing probability-----

%---removing from compartment A-----
vAE = 0;
vAA = 0;
vAP = 0;
vAI = 0;
if tL>1 && A(currentTime)>0
	for i = 1:tL
    	vAE = vAE+P_EE(rE_star,i+t0)*P_AA(rA_star,i)*A(currentTime-i)/A(currentTime);
		vAA = vAA+P_AE(qEA,rE_star,rA_star,i+t0)*P_AA(rA_star,i)*A(currentTime-i)/A(currentTime);
		vAP = vAP+P_PE(qEP,rE_star,rP_star,i+t0)*P_AA(rA_star,i)*A(currentTime-i)/A(currentTime);
		vAI = vAI+P_IE(qPI,qEP,rI_star,rE_star,rP_star,i+t0)*P_AA(rA_star,i)*A(currentTime-i)/A(currentTime);
	end
else
	i = tL;
	vAE = P_EE(rE_star,i+t0)*P_AA(rA_star,i);
	vAA = P_AE(qEA,rE_star,rA_star,i+t0)*P_AA(rA_star,i);
	vAP = P_PE(qEP,rE_star,rP_star,i+t0)*P_AA(rA_star,i);
	vAI = P_IE(qPI,qEP,rI_star,rE_star,rP_star,i+t0)*P_AA(rA_star,i);

	vAE = 1;
	vAA = 1;
	vAP = 1;
	vAI = 1;
end


%---removing from compartment P-----
vPE = 0;
vPA = 0;
vPP = 0;
vPI = 0;
if tL>1 && P(currentTime)>0
	for i = 1:tL
    	vPE = vPE+P_EE(rE_star,i+t0)*P_PP(rP_star,i)*P(currentTime-i)/P(currentTime);
		vPA = vPA+P_AE(qEA,rE_star,rA_star,i+t0)*P_PP(rP_star,i)*P(currentTime-i)/P(currentTime);
		vPP = vPP+P_PE(qEP,rE_star,rP_star,i+t0)*P_PP(rP_star,i)*P(currentTime-i)/P(currentTime);
		vPI = vPI+P_IE(qPI,qEP,rI_star,rE_star,rP_star,i+t0)*P_PP(rP_star,i)*P(currentTime-i)/P(currentTime);
	end
else
	i = tL;
	vPE = P_EE(rE_star,i+t0)*P_PP(rP_star,i);
	vPA = P_AE(qEA,rE_star,rA_star,i+t0)*P_PP(rP_star,i);
	vPP = P_PE(qEP,rE_star,rP_star,i+t0)*P_PP(rP_star,i);
	vPI = P_IE(qPI,qEP,rI_star,rE_star,rP_star,i+t0)*P_PP(rP_star,i);

	vPE = 1;
	vPA = 1;
	vPP = 1;
	vPI = 1;

end
       
%---removing from compartment I-----

tmpIE1 = 0;
tmpIE2 = 0;
tmpIA1 = 0;
tmpIA2 = 0;
tmpIP1 = 0;
tmpIP2 = 0;
tmpII1 = 0;
tmpII2 = 0;

Indicator = P+I;
Zero = Indicator==0;
if tL>1 && I(currentTime)>0 && sum(Zero)==0
	for i = 1:tL
		tmpIE1 = tmpIE1+P_EE(rE_star,i+t0)*P_IP(qPI,rI_star,rP_star,i)*P(currentTime-i)^2/I(currentTime)/(P(currentTime-i)+I(currentTime-i)); %--K=P---
		tmpIE2 = tmpIE2+P_EE(rE_star,i+t0)*P_II(rI_star,i)*I(currentTime-i)^2/I(currentTime)/(P(currentTime-i)+I(currentTime-i)); %--K=I---
		tmpIA1 = tmpIA1+P_AE(qEA,rE_star,rA_star,i+t0)*P_IP(qPI,rI_star,rP_star,i)*P(currentTime-i)^2/I(currentTime)/(P(currentTime-i)+I(currentTime-i)); %--K=P---
		tmpIA2 = tmpIA2+P_AE(qEA,rE_star,rA_star,i+t0)*P_II(rI_star,i)*I(currentTime-i)^2/I(currentTime)/(P(currentTime-i)+I(currentTime-i)); %--K=I---
		tmpIP1 = tmpIP1+P_PE(qEP,rE_star,rP_star,i+t0)*P_IP(qPI,rI_star,rP_star,i)*P(currentTime-i)^2/I(currentTime)/(P(currentTime-i)+I(currentTime-i)); %--K=P---
		tmpIP2 = tmpIP2+P_PE(qEP,rE_star,rP_star,i+t0)*P_II(rI_star,i)*I(currentTime-i)^2/I(currentTime)/(P(currentTime-i)+I(currentTime-i)); %--K=I---
		tmpII1 = tmpII1+P_IE(qPI,qEP,rI_star,rE_star,rP_star,i+t0)*P_IP(qPI,rI_star,rP_star,i)*P(currentTime-i)^2/I(currentTime)/(P(currentTime-i)+I(currentTime-i)); %--K=P---
		tmpII2 = tmpII2+P_IE(qPI,qEP,rI_star,rE_star,rP_star,i+t0)*P_II(rI_star,i)*I(currentTime-i)^2/I(currentTime)/(P(currentTime-i)+I(currentTime-i)); %--K=I---
	end
	vIE = tmpIE1+tmpIE2;
	vIA = tmpIA1+tmpIA2;
	vIP = tmpIP1+tmpIP2;
	vII = tmpII1+tmpII2;
else
	i = tL;
	tmpIE1 = P_EE(rE_star,i+t0)*P_IP(qPI,rI_star,rP_star,i); %--K=P---
	tmpIE2 = P_EE(rE_star,i+t0)*P_II(rI_star,i); %--K=I---
	tmpIA1 = P_AE(qEA,rE_star,rA_star,i+t0)*P_IP(qPI,rI_star,rP_star,i); %--K=P---
	tmpIA2 = P_AE(qEA,rE_star,rA_star,i+t0)*P_II(rI_star,i); %--K=I---
	tmpIP1 = P_PE(qEP,rE_star,rP_star,i+t0)*P_IP(qPI,rI_star,rP_star,i); %--K=P---
	tmpIP2 = P_PE(qEP,rE_star,rP_star,i+t0)*P_II(rI_star,i); %--K=I---
	tmpII1 = P_IE(qPI,qEP,rI_star,rE_star,rP_star,i+t0)*P_IP(qPI,rI_star,rP_star,i); %--K=P---
	tmpII2 = P_IE(qPI,qEP,rI_star,rE_star,rP_star,i+t0)*P_II(rI_star,i); %--K=I---
	vIE = tmpIE1+tmpIE2;
	vIA = tmpIA1+tmpIA2;
	vIP = tmpIP1+tmpIP2;
	vII = tmpII1+tmpII2;

	vIE = 1;
	vIA = 1;
	vIP = 1;
	vII = 1;
end

%---rownames: A, P, I; colnames: E, A, P, I;------
tmp = vAE+vAA+vAP+vAI;
if tmp==0
	tmp = 1;
end
vAE = vAE/tmp;
vAA = vAA/tmp;
vAP = vAP/tmp;
vAI = vAI/tmp;

tmp = vPE+vPA+vPP+vPI;
if tmp==0
	tmp = 1;
end
vPE = vPE/tmp;
vPA = vPA/tmp;
vPP = vPP/tmp;
vPI = vPI/tmp;

tmp = vIE+vIA+vIP+vII;
if tmp==0
	tmp = 1;
end
vIE = vIE/tmp;
vIA = vIA/tmp;
vIP = vIP/tmp;
vII = vII/tmp;

pij = [vAE vAA vAP vAI; vPE vPA vPP vPI; vIE vIA vIP vII];

end