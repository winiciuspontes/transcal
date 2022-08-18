import sympy  as sp
import numpy as np 
import pandas as pd 
import math as mt 
import matplotlib.pyplot as plt
import transcal as tc
import scipy as spy

#CONSTANTES UTILIZADAS:
k = 143.5
m_ponto_3 = 0.243
m_ponto_comb = 0.00398 * 0.19 #Tirar dúvida sobre o m_comb (nao sabemos qual é o valor)
m_ponto_zp = m_ponto_3  * 0.19
m_ponto_zs = m_ponto_3 * 0.0961
# m_ponto_fenda_zd = 0.247 * 0.19  #Não estou certo se o valor é esse, retirei da pagina 75 do TG do Gasturb
T3 = 498.42
T4 = 1100
P1 = 100 * 1000 #Tabela  GasTurb
P2 = 99  * 1000 #Tabela  GasTurb
P3 = 354.420 * 1000 
delta_P3_4 = (354.420 - 333.155) * 1000 
delta = 0.8 # Adotei para o calculo o Orificio de adminissao de canto vivo, pag 70 do TG
theta = 73 * (10**6) 
d_int = 0.024 # Valor obtido da tabela 4
m_ponto_h_zp = 0.254 
razao_delta_P3_4_q_ref = 20
razao_delta_P3_4_q_p3 = 0.06
print("Valor de m_ponto_3", m_ponto_3)
print("Valor de m_ponto_zp", m_ponto_zp)
print("Valor de m_ponto_zs", m_ponto_zs)

#Eq.30, 31, 32, 33
phi_zp = 1.249
phi_zs = 0.8
phi_zd = 0.237
phi_global = 0.237
phi_pobre = 0.45
phi_rico = 2.553
phi_estequiometrica = 0.06818

#VALORES DO PHI DE ACORDO COM AS EQUAÇÕES DO TG
phi_pobre = tc.phi_pobre_function(T3)
phi_rico = tc.phi_rico_function(T3)
# m_ponto_zs = (phi_global/phi_zs) * m_ponto_3

#Calculando o delta T_phi com phi = 1 (Eq.39)
delta_T_phi = 2185 - (0.5 * T3)
delta_T_zp = 2830 - (800 * phi_zp) #Aproximado para quando T3=400K  - Não sabemos que de onde vem esse valor de 2830
delta_T_zs = 1600 - (0.5 * T3)

#Calcunado A_ref, d_ref, A_ft
b = tc.b_function(phi_zp) #Eq. antes da 30
d_ref = 0.03863 # Valor obtido da tabela 4)
# a_ref = tc.area_tranferencia_aerodinamica(k, m_ponto_3, T3, P3, razao_delta_P3_4_q_ref, razao_delta_P3_4_q_p3)
a_ref = 0.0076
a_ft = 0.004256 #Eq. 34 - No trabalho do Caio ele multiplicou por 10 ^ 6

d_ref = tc.altura_referencia(a_ref, d_int)
d_ft = tc.altura_tubo_chama(a_ft, d_int, d_ref)

l_zp = (3/4) * d_ft #Enunciado antes da Eq.35 - Pegamos o l_zp_max, para pegar o minimo trocamos o 3/4 por 2/3
l_zs = (1/2) * d_ft #Enunciado antes da Eq.35
l_zd = (3/2) * d_ft #Enunciado antes da Eq.35
l_zr = (1/2) * d_ft #Enunciado antes da Eq.35
lcc = tc.comprimento_camara_combustao(l_zp,l_zs,l_zd) #Eq.35 - Nosso resultado bate com a Tabela 11

quanti_relativa_ar = tc.quant_ar_zona_secundaria(phi_global, phi_rico, phi_zs)

#Calculando a Eq.37
m_ponto_arref = tc.porcentagem_ar_resfriamento(T3, m_ponto_3) #Eq.37 - Equacao nao confirmada(não achamos o ela no Tg)
print("Valor de m_ponto_arref", m_ponto_arref)



m_ponto_zd = tc.m_ponto_zd_funcao(m_ponto_zp, m_ponto_zs, m_ponto_arref, m_ponto_3) # m_ponto_arref desconsiderado para não dar negativo
m_ponto_zd = 0.1749 # m_ponto_arref desconsiderado para não dar negativo - 0.04360957259559913
print("Valor de m_ponto_zd", m_ponto_zd)


#Eq. 40 - calculando o eta
eta_zr = tc.eta_zr(T3,P3)
print("Valor de eta_zr", eta_zr)

t_max_zr = tc.t_max_zr(T3, eta_zr, delta_T_phi)
print("Valor de t_max_zr", t_max_zr)

t_media_zr = tc.t_media_zr(T3, t_max_zr)
print("Valor de t_media_zr", t_media_zr)

#Eq.43 - Calculando o eta_zp:

eta_zp = tc.eta_zp(T3, P3)
print("Valor de eta_zp", eta_zp)
#Eq.42 - Calculando a Tsaida_zp - Caio trocou a ordem novamente:
t_saida_zp = tc.t_saida_zp(T3, eta_zp, delta_T_zp)
print("Valor de t_saida_zp", t_saida_zp)
#Eq.45 - Calculando o Tsaida_zs rico
eta_zs_rico = tc.eta_zs(phi_zs)
#Eq. 44 - Calculando a Tsaida_zs rico - Caio trocou a ordem das eq novamente:
t_saida_zs_rico = tc.t_saida_zs(T3, eta_zs_rico, delta_T_zs)
print("Valor de t_saida_zs_rico", t_saida_zs_rico)
v_zs = tc.v_zs(a_ft, l_zs)
print("Valor de v_zs", v_zs)
# Eq. 48 - Calculando o psi_maiusculo:
#Vamos calcular o psi_t3 porque o psi_300 é um caso especifico. Perg para o professor
psi_t3 = tc.phi_t3(phi_zs, m_ponto_comb, v_zs, P3, T3) 
# psi_300 = tc.phi_t300(m_ponto_comb, P3, phi_zs,  v_zs) 
#Eq. 47 - Calculando o D*:
d_asterisco = tc.d_ast(delta_P3_4,P3)
#Eq.46 - Calculando o eta_zs pobre(caio trocou a ordem novamente):
eta_zs_pobre = tc.eta_zs_pobre(psi_t3, phi_zs, d_asterisco)
#Eq. 44 para Tsaida_zs pobre - Caio trocou a ordem das eq novamente:
t_saida_zs_pobre= tc.t_saida_zs(T3, eta_zs_pobre, delta_T_zs) 
#Eq. veio da tabela 12 
t_saida_zd = 1134
print(t_saida_zd)



# Calculando o tamanho:
x1 = np.linspace(0, l_zr, 25)
x2 = np.linspace(l_zr, l_zp, 25)
x3 = np.linspace(l_zp, (l_zp + l_zs), 25)
x4 = np.linspace((l_zp + l_zs), lcc, 25)


x_total = np.linspace(0, lcc, 100)

tg1_lista = []
tg2_lista = []
tg3_lista = []
tg4_lista = []
tg_lista = []

# Fazendo o looping para conseguir plotar o gráficos das temperaturas
for c in range(len(x1)):
    Tg1 = t_media_zr
    tg_lista.append(Tg1)
    tg1_lista.append(Tg1)

for d in range(len(x2)):
    Tg2 = tc.tg_2(t_media_zr, t_saida_zp, l_zp,l_zr, x2[d])
    tg2_lista.append(Tg2)

for e in range(len(x3)):
    Tg3 = tc.tg_3(t_saida_zp, t_saida_zs_rico, l_zs,l_zp, x3[e])
    tg3_lista.append(Tg3)

for f in range(len(x4)):
    Tg4 = tc.tg_4(t_saida_zs_rico, t_saida_zd, l_zd, x4[f], l_zp, l_zs)
    tg4_lista.append(Tg4)

tg_total = tg1_lista + tg2_lista + tg3_lista + tg4_lista
x_total_sem_resf =tc.x_total_sem_resf
tg_total_sem_resf=tc.tg_total_sem_resf

plt.title('Temperatura Gases e Temperatura parede sem resfriamento')
plt.plot(x_total_sem_resf,tg_total_sem_resf,'y')
plt.plot(x_total, tg_total)
plt.ylabel('Temperatura (K)')
plt.xlabel('Comprimento da camara (mm)')
plt.grid()
plt.show()
plt.ylabel('Tg')
plt.xlabel('x [mm]')

mg1_ponto = []
mg2_ponto = []
mg3_ponto = []
mg4_ponto = []

man_1_lista = []
man_2_lista = []
man_3_lista = []
man_4_lista = []


for indice, valor in enumerate(x1):
    mg1_ponto_calculo = tc.mg1_ponto(m_ponto_zp)
    mg1_ponto.append(mg1_ponto_calculo)
    man_1_calculo = m_ponto_3 - mg1_ponto_calculo
    man_1_lista.append(man_1_calculo)


for indice, valor in enumerate(x2):
    mg2_ponto_calculo = tc.mg2_ponto(mg1_ponto[indice], m_ponto_zp, valor, l_zr,l_zp)
    mg2_ponto.append(mg2_ponto_calculo)
    man_2_calculo = m_ponto_3 - mg2_ponto_calculo
    man_2_lista.append(man_2_calculo)

for indice, valor in enumerate(x3):
    mg3_ponto_calculo = tc.mg3_ponto(mg2_ponto[indice], m_ponto_zs, valor, l_zp, l_zs)
    mg3_ponto.append(mg3_ponto_calculo)
    man_3_calculo = m_ponto_3 - mg3_ponto_calculo
    man_3_lista.append(man_3_calculo)


for indice, valor in enumerate(x4):
    mg4_ponto_calculo = tc.mg4_ponto(mg3_ponto[indice], m_ponto_zd, valor, l_zp, l_zs, l_zd)
    mg4_ponto.append(mg4_ponto_calculo)
    man_4_calculo = m_ponto_3 - mg4_ponto_calculo
    man_4_lista.append(man_4_calculo)


mg_total = mg1_ponto + mg2_ponto + mg3_ponto + mg4_ponto 
man_total = man_1_lista + man_2_lista + man_3_lista + man_4_lista 

plt.title('Fluxo de massas')
plt.ylabel('Vazão mássica')
plt.xlabel('Comprimento da camara (mm)')
plt.grid()
plt.plot(x_total, mg_total)
plt.show()


#########################################################################################################################


#Posição das Fendas Escolhidas
posicao_fenda_1 = 29*1e-3
x_t = np.linspace(1e-4, lcc, 100)
id_1 = np.argmin(np.abs(x_t-posicao_fenda_1))
m_ponto_gases_fenda_1 = mg_total[id_1] 
m_ponto_anular_fenda_1 = man_total[id_1]
t_gases_fenda_1 = tg_total[id_1]

#Posição da Fenda

# Área da Fenda 1
s_1 = 0.001
t_1 = 0.002
a_fenda_1 = 2*np.pi*s_1*(d_ref + d_ft)
A_an = a_ref-a_ft 
m_fenda_1 = m_ponto_anular_fenda_1*(a_fenda_1/A_an)

index = 0
for i in x_t:
    if i> posicao_fenda_1:
        mg_total[index] = mg_total[index] + m_fenda_1
        man_total[index] = man_total[index] - m_fenda_1
    index+=1

#Fenda 2
posicao_fenda_2= 50*1e-3
id_2 = np.argmin(np.abs(x_t-posicao_fenda_2))
m_ponto_gases_fenda_2 = mg_total[id_2]
m_ponto_anular_fenda_2 = man_total[id_2]
t_gases_fenda_2  = tg_total[id_2]

# Área da Fenda 2
s_2 = 0.001
t_2 = 0.002 
A_fenda_2 = 2*np.pi*s_2*(d_ref + d_ref)
m_fenda_2 = m_ponto_anular_fenda_2*(A_fenda_2/A_an)


index = 0
for i in x_t:
    if i> posicao_fenda_2:
        mg_total[index] = mg_total[index] + m_fenda_2
        man_total[index] = man_total[index] - m_fenda_2
    index+=1

# Eficiência do resfriamento Fenda 1
pu_an_1 = m_ponto_anular_fenda_1/A_an 
pu_g_1 = m_ponto_gases_fenda_1/a_ft 
visc_ar = (0.03863 + 0.00749*T3 - (5.8564*(10**(-6))*(T3**2)) + (2.7769*(10**(-9))*(T3**3)) - (4.600774*(10**(-13))*(T3**4)))*(10**(-5)) 
visc_g_1 = (0.03863 + 0.00749*t_gases_fenda_1-5.8564*(10**-6)*(t_gases_fenda_1**2) + 2.7769*(10**-9)*(t_gases_fenda_1**3) - 4.600774*(10**-13)*(t_gases_fenda_1**4))*(10**-5) #viscosidade dinamica do gás no interior do tubo de chama
razao_m_1 = pu_an_1/pu_g_1
print(f'Razao m fenda 1: {razao_m_1}')

# Eficiência do resfriamento Fenda 2
pu_an_2 = m_ponto_anular_fenda_2/A_an 
pu_g_2 = m_ponto_gases_fenda_2/a_ft 
visc_g_2 =(0.03863 + 0.00749*t_gases_fenda_2-5.8564*(10**-6)*(t_gases_fenda_2**2) + 2.7769*(10**-9)*(t_gases_fenda_2**3) - 4.600774*(10**-13)*(t_gases_fenda_2**4))*(10**-5) #viscosidade dinamica do gás no interior do tubo de chama
razao_m_2 = pu_an_2/pu_g_2
#razao_m_2 = 0.52
print(f'Razao m fenda 2: {razao_m_2}')



#Temperatura dos gases proximo da parede

# Cálculo da Temperatura na Parede ao Longo do Tubo de Chama 
index = 0
vetor_temp = []
graph_eff = []

for posX in x_t:
	m_ponto_gas = mg_total[index] 
	m_ponto_an = man_total[index]
	t_gases = tg_total[index]

	if posX < posicao_fenda_1:
		visc_g =(0.03863 + 0.00749*t_gases-5.8564*(10**-6)*(t_gases**2) + 2.7769*(10**-9)*(t_gases**3) - 4.600774*(10**-13)*(t_gases**4))*(10**-5) 
		posXEff = posX
		nr = 0
	elif posX >= posicao_fenda_1 and posX < posicao_fenda_2:
		visc_g = visc_g_1 
		posXEff = posX
		posRel = posX - posicao_fenda_1 + 0.001
		razao_m = razao_m_1

		if razao_m>0.5 and razao_m<=1.3:
			nr = 1.10*(razao_m**0.65)*((visc_ar/visc_g)**0.15)*((posRel/s_1)**(-0.2))*((t_1/s_1)**(-0.2)) 
		elif razao_m>1.3 and razao_m<=4:
			nr = 1.28*((visc_ar/visc_g)**0.15)*((posRel/s_1)**(-0.2))*((t_1/s_1)**(-0.2))
		else:
			nr = 0
			print("Favor verificar intervalo.")
	else:
		visc_g = visc_g_2 
		posXEff = posX
		razao_m = razao_m_2
		posRel = posX - posicao_fenda_2 + 0.001

		if razao_m>0.5 and razao_m<=1.3:
			nr = 1.10*(razao_m**0.65)*((visc_ar/visc_g)**0.15)*((posRel/s_2)**(-0.2))*((t_2/s_2)**(-0.2)) 
		elif razao_m>1.3 and razao_m<=4:  #Editei aqui
			nr = 1.28*((visc_ar/visc_g)**0.15)*((posRel/s_2)**(-0.2))*((t_2/s_2)**(-0.2))
		else:
			nr = 0
			print("Favor verificar intervalo.")

	graph_eff.append(nr)
	Tg_w = t_gases - nr * (t_gases - T3)

	#TRANSFERENCIA DE CALOR NA PAREDE
	# RADIAÇÃO
	#T_w1 Temperatura na superficie interna
	#T_w2 Temperatura na superficie externa
	T_w1, T_w2 = sp.symbols('T_w1,T_w2')

	k_w = 26 #[W/m.K] condutividade térmica
	t_w = 0.0005 
	sigma = 5.67*10**(-8) #Stefan-Boltzmann
	e_w = 0.4
	q = 0.00398/0.243
	l_b_int = d_ft 
	l_b_ext = 1.2 * d_ft 
	e_g_int = 1 - np.exp(-0.290*P3*(q*l_b_int)**0.5*t_gases**(-1.5)) 
	e_g_ext = 1 - np.exp(-0.290*P3*(q*l_b_ext)**0.5*t_gases**(-1.5))
	Z = 0.4 #para alumínio - página 66 do TG

	R_1 = 0.5*sigma*(1+e_w)*e_g_int*(t_gases**1.5)*((t_gases**2.5)-(T_w1**2.5)) 
	R_2 = Z*sigma*(T_w2**4 - T3**4)       

	#K1-2 :
	K_1_2 = (k_w/t_w)*(T_w1 - T_w2) 

	#CONVECÇÃO
	k_g = 5.92657*10**(-4)+9.80957*(10**(-5))*t_gases-4.89398*(10**(-8))*(t_gases**2)+1.5011410*(10**(-11))*(t_gases**3)
	Re_x = (m_ponto_gas*posXEff)/(a_ft*visc_g)  

	if posX <= posicao_fenda_1:
		D_L = 4*a_ft/d_ft
		C_1 = 0.017*(k_g/(D_L**0.2))*((m_ponto_gas/(a_ft*visc_g))**0.8)*(t_gases-T_w1)
	elif posX >= posicao_fenda_1 and posX < posicao_fenda_2:
		if razao_m > 0.5 and razao_m <= 1.3:
			C_1 = 0.069*(k_g/posXEff)*(Re_x**0.7)*(Tg_w-T_w1)
		elif razao_m > 1.3 and razao_m < 4.0:
			C_1 = 0.10*(k_g/posXEff)*Re_x**0.8*(posXEff/s_1)**(-0.36)*(Tg_w-T_w1)
	else:
		if razao_m > 0.5 and razao_m <= 1.3:
			C_1 = 0.069*(k_g/posXEff)*(Re_x**0.7)*(Tg_w-T_w1)
		elif razao_m > 1.3 and razao_m < 4.0:
			C_1 = 0.10*(k_g/posXEff)*Re_x**0.8*(posXEff/s_2)**(-0.36)*(Tg_w-T_w1)
		
	D_an = d_ref - d_ft
	k_a = 5.92657*(10**(-4))+9.80957*(10**(-5))*T3-4.89398*(10**(-8))*(T3**2)+1.5011410*(10**(-11))*(T3**3) 
	C_2 = 0.020*(k_a/(D_an**0.2))*((m_ponto_an/(A_an*visc_ar))**0.8)*(T_w2 - T3) 

	eq1 = sp.Eq(R_1+C_1-R_2-C_2,0)
	eq2 = sp.Eq(R_1+C_1-K_1_2,0)
	solver = sp.nsolve((eq1,eq2),(T_w1, T_w2), (1500,800))
	tfinal = max(solver[0],solver[1])
	vetor_temp.append(tfinal)
	print(rf"Interação {index} T_w1 = {round(solver[0],2):,.2f}, T_w2 = {round(solver[1],2):,.2f}")
	index+=1

lcc = tc.comprimento_camara_combustao(l_zp,l_zs,l_zd)



plt.figure(4,figsize=(12, 7), dpi=80)
plt.plot(x_t*(10**3),vetor_temp,'b',label = 'Temperatura máxima na superfície')
plt.hlines(800,0,65,linestyles='--', label='Limite do Material', color='g')
plt.title('Temperatura da parede ao longo do tubo de chama')
plt.ylabel('Temperatura (K)')
plt.xlabel('Comprimento da camara (mm)')
plt.ylim(400,1000)
plt.xlim(0,60)
plt.grid()
plt.show()

plt.plot(x_t*(10**3),graph_eff,'b',label = 'Eficiência')
plt.vlines(posicao_fenda_1*1e3,0,1,colors='k', linestyles='--')
plt.vlines(posicao_fenda_2*1e3,0,1,colors='y', linestyles='--')
plt.title('Eficiência do resfriamento')
plt.ylabel('Eficiência')
plt.xlabel('Comprimento da camara (mm)')
plt.grid()
plt.show()








