from sqlalchemy import false
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
m_ponto_zs = m_ponto_3 * 0.19
m_ponto_fenda_zd = 0.247 * 0.19  #Não estou certo se o valor é esse, retirei da pagina 75 do TG do Gasturb
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


#CALCULO DOS PHIS DE ACORDO COM O CAIO:

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
m_ponto_zs = (phi_global/phi_zs) * m_ponto_3
print("Valor de m_ponto_zs", m_ponto_zs)
#Calculando o delta T_phi com phi = 1 (Eq.39)
delta_T_phi = 2185 - (0.5 * T3)
delta_T_zp = 2830 - (800 * phi_zp) #Aproximado para quando T3=400K  - Não sabemos que de onde vem esse valor de 2830
delta_T_zs = 1600 - (0.5 * T3)


b = tc.b_function(phi_zp) #Eq. antes da 30
d_ref = 0.03863 # Valor obtido da tabela 4)

#Acima da Eq. 30 
b = 245*(1.39 + np.log(phi_zp))

# a_ref_quimico = tc.area_tranferencia_quimica(P3, d_ref, T3, b, m_ponto_3, theta)
# a_ref_aerodinamica = tc.area_tranferencia_aerodinamica(k, m_ponto_3, T3, P3, razao_delta_P3_4_q_ref, razao_delta_P3_4_q_p3)
# print("Valor de a_ref_aerodinamica", a_ref_aerodinamica)
a_ref_aerodinamica = 0.00759996 # A_ref do Caio

# d_ref = tc.altura_referencia(a_ref_aerodinamica, d_int) # Eq. Tabela 6 
d_ref = 0.03863 # Eq. Tabela 6 
print("Valor de d_ref", d_ref)
# a_ft = tc.area_ref_trans(a_ref_aerodinamica) #Eq. 34 - No trabalho do Caio ele multiplicou por 10 ^ 6
a_ft = 0.004256 #Eq. 34 - No trabalho do Caio ele multiplicou por 10 ^ 6

print("Valor de a_ft", a_ft)
# d_ft = tc.altura_tubo_chama(a_ft, d_int, d_ref) #Eq. Tabela 6 (para o caso anular)
d_ft = 0.02163
print("Valor de d_ft", d_ft)

#Calculando o Lcc com os parametros necessarios:
l_zp = (3/4) * d_ft #Enunciado antes da Eq.35 - Pegamos o l_zp_max, para pegar o minimo trocamos o 3/4 por 2/3
l_zs = (1/2) * d_ft #Enunciado antes da Eq.35
l_zd = (3/2) * d_ft #Enunciado antes da Eq.35
l_zr = (1/2) * d_ft #Enunciado antes da Eq.35
lcc = tc.comprimento_camara_combustao(l_zp,l_zs,l_zd) #Eq.35 - Nosso resultado bate com a Tabela 11

#Eq. 36 (phi_global_rico/phi_zs) é convencionada da literatura como sendo 0.8 
razao_phi_global_rico_phi_zs = 0.8 

#Calculando a Eq.37
m_ponto_arref = tc.porcentagem_ar_resfriamento(T3, m_ponto_3) #Eq.37 - Equacao nao confirmada(não achamos o ela no Tg)
print("Valor de m_ponto_arref", m_ponto_arref)

#Calculando a Eq.38
m_ponto_zd = tc.m_ponto_zd_funcao(m_ponto_zp, m_ponto_zs, 0, m_ponto_3) # m_ponto_arref desconsiderado para não dar negativo
print("Valor de m_ponto_zd", m_ponto_zd)


#Eq. 40 - calculando o eta
eta_zr = tc.eta_zr(T3,P3)
print("Valor de eta_zr", eta_zr)
#Eq.39 - As eq estavam fora de ordem
t_max_zr = tc.t_max_zr(T3, eta_zr, delta_T_phi)
print("Valor de t_max_zr", t_max_zr)
#Eq.41 - Calculando a T_med - dando dif de 10 graus :
t_media_zr = tc.t_media_zr(T3, t_max_zr)
print("Valor de t_media_zr", t_media_zr)
# delta_T_zp = 3000 - t_media_zr

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


#Eq. 49 - Calculando o V_zs:
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

#Tabela 7


tg_lista = []

x1 = np.linspace(0, l_zr, 50)
x2 = np.linspace(l_zr, l_zp, 50)
x3 = np.linspace(l_zp, (l_zp + l_zs), 50)
x4 = np.linspace((l_zp + l_zs), lcc, 50)


x_total = np.linspace(0, lcc, 200)

tg1_lista = []
tg2_lista = []
tg3_lista = []
tg4_lista = []
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
#Gráfico aparentemente certo - bateu com o TG do caio (grafico 4)

# plt.figure(figsize = (15,10))
# plt.plot(x1*1000, tg1_lista, label='Zona de Recirc.')
# plt.plot(x2*1000, tg2_lista, label='Zona Primária')
# plt.plot(x3*1000, tg3_lista, label='Zona Secundária')
# plt.plot(x4*1000, tg4_lista, label='Zona Secundária')
# plt.show()

s = 0.415/100
area_fenda = tc.area_fenda(d_ref, d_ft, s) # A altura da fenda é a soma das alturas de referencia e do tubo de chama 

#Eq. 51 (confirmar com o grupo se vamos colocar a eq 51 como sendo um array igual o de cima) - feito dentro do for

a_an = tc.A_an_int(d_int, d_ref, d_ft)

print("Valor de area_fenda", area_fenda)
print("Valor de a_an", a_an)


mg_ponto_zr_lista = []
mg_ponto_zp_lista = []
mg_ponto_zs_lista = []
mg_ponto_zd_lista = []

man_1_lista = []
man_2_lista = []
man_3_lista = []
man_4_lista = []

# m_ponto_g1 =  tc.mg1_ponto(m_ponto_zp)
# m_ponto_g2 =  tc.mg2_ponto(m_ponto_g1, m_ponto_zp, 0.014, l_zr,l_zp)
# m_ponto_g3 =  tc.mg3_ponto(m_ponto_g2, m_ponto_zs, 0.033, l_zp, l_zs)
# m_ponto_g4 =  tc.mg4_ponto(m_ponto_g3, m_ponto_zd, 0.043,l_zp, l_zs, l_zd)

m_ponto_fenda1_lista = []
m_ponto_fenda2_lista = []
m_ponto_fenda3_lista = []
m_ponto_fenda4_lista = []

razao_rho_an_mi_an1_lista = []
razao_rho_an_mi_an2_lista = []
razao_rho_an_mi_an3_lista = []
razao_rho_an_mi_an4_lista = []

razao_rho_g_mi_g_lista1 = []
razao_rho_g_mi_g_lista2 = []
razao_rho_g_mi_g_lista3 = []
razao_rho_g_mi_g_lista4 = []


# razao_rho_g_mi_g = tc.razao_rho_g_mi_an(m_ponto_g1, a_ft)
# razao_rho_g_mi_g = tc.razao_rho_g_mi_an(m_ponto_g1, a_ft)
# razao_rho_g_mi_g = tc.razao_rho_g_mi_an(m_ponto_g1, a_ft)
# razao_rho_g_mi_g = tc.razao_rho_g_mi_an(m_ponto_g1, a_ft)


for indice,valor in enumerate (x1):
    mg_ponto_zr = tc.mg1_ponto(m_ponto_zp)
    man_1_calculo = m_ponto_3 - mg_ponto_zr
    mg_ponto_zr_lista.append(mg_ponto_zr)
    man_1_lista.append(man_1_calculo)
    m_ponto_fenda1 = tc.m_ponto_fenda_funcao(area_fenda, man_1_calculo, a_an)
    m_ponto_fenda1_lista.append(m_ponto_fenda1)
    razao_rho_an_mi_an1_calculo = tc.razao_rho_an_mi_an(m_ponto_fenda1, area_fenda)
    razao_rho_an_mi_an1_lista.append(razao_rho_an_mi_an1_calculo)
    razao_rho_g_mi_g1 = tc.razao_rho_g_mi_an(mg_ponto_zr, a_ft)
    razao_rho_g_mi_g_lista1.append(razao_rho_g_mi_g1)

for indice, valor in enumerate (x2):
    mg_ponto_zp = tc.mg2_ponto(mg_ponto_zr_lista[indice], m_ponto_zp, valor, l_zr,l_zp)
    man_2_calculo = m_ponto_3 - mg_ponto_zp
    mg_ponto_zp_lista.append(mg_ponto_zp)
    man_2_lista.append(man_2_calculo)
    m_ponto_fenda2 = tc.m_ponto_fenda_funcao(area_fenda, man_2_calculo, a_an)
    m_ponto_fenda2_lista.append(m_ponto_fenda2)
    razao_rho_an_mi_an2_calculo = tc.razao_rho_an_mi_an(m_ponto_fenda2, area_fenda)
    razao_rho_an_mi_an2_lista.append(razao_rho_an_mi_an2_calculo)
    razao_rho_g_mi_g2 = tc.razao_rho_g_mi_an(mg_ponto_zp, a_ft)
    razao_rho_g_mi_g_lista2.append(razao_rho_g_mi_g2)

for indice, valor in enumerate (x3):
    mg_ponto_zs = tc.mg3_ponto(mg_ponto_zp_lista[indice], m_ponto_zs, valor, l_zp, l_zs)
    man_3_calculo = m_ponto_3 - mg_ponto_zs
    mg_ponto_zs_lista.append(mg_ponto_zs)
    man_3_lista.append(man_3_calculo)
    m_ponto_fenda3 = tc.m_ponto_fenda_funcao(area_fenda, man_3_calculo, a_an)
    m_ponto_fenda3_lista.append(m_ponto_fenda3)
    razao_rho_an_mi_an3_calculo = tc.razao_rho_an_mi_an(m_ponto_fenda3, area_fenda)
    razao_rho_an_mi_an3_lista.append(razao_rho_an_mi_an3_calculo)
    razao_rho_g_mi_g3 = tc.razao_rho_g_mi_an(mg_ponto_zs, a_ft)
    razao_rho_g_mi_g_lista3.append(razao_rho_g_mi_g3)

m_ponto_zd = 0.12329510485804025
for indice,valor in enumerate (x4):
    mg_ponto_zd = tc.mg4_ponto(mg_ponto_zs_lista[indice], m_ponto_zd, valor, l_zp, l_zs, l_zd)
    man_4_calculo = m_ponto_3 - mg_ponto_zd
    mg_ponto_zd_lista.append(mg_ponto_zd)
    man_4_lista.append(man_4_calculo)
    m_ponto_fenda4 = tc.m_ponto_fenda_funcao(area_fenda, man_4_calculo, a_an)
    m_ponto_fenda4_lista.append(m_ponto_fenda4)
    razao_rho_an_mi_an4_calculo = tc.razao_rho_an_mi_an(m_ponto_fenda4, area_fenda)
    razao_rho_an_mi_an4_lista.append(razao_rho_an_mi_an4_calculo)
    razao_rho_g_mi_g4 = tc.razao_rho_g_mi_an(mg_ponto_zd, a_ft)
    razao_rho_g_mi_g_lista4.append(razao_rho_g_mi_g4)

mg_total = mg_ponto_zr_lista + mg_ponto_zp_lista + mg_ponto_zs_lista + mg_ponto_zd_lista
m_an_total = man_1_lista + man_2_lista + man_3_lista + man_4_lista
m_ponto_fenda_total = m_ponto_fenda1_lista + m_ponto_fenda2_lista + m_ponto_fenda3_lista + m_ponto_fenda4_lista 
razao_rho_an_mi_an_total = razao_rho_an_mi_an1_lista + razao_rho_an_mi_an2_lista +razao_rho_an_mi_an2_lista + razao_rho_an_mi_an4_lista 
razao_rho_g_mi_g_total = razao_rho_g_mi_g_lista1 + razao_rho_g_mi_g_lista2 + razao_rho_g_mi_g_lista3 + razao_rho_g_mi_g_lista4



mg_total = mg_ponto_zr_lista + mg_ponto_zp_lista + mg_ponto_zs_lista + mg_ponto_zd_lista
m_an_total = man_1_lista + man_2_lista + man_3_lista + man_4_lista
m_ponto_fenda_total = m_ponto_fenda1_lista + m_ponto_fenda2_lista + m_ponto_fenda3_lista + m_ponto_fenda4_lista 
razao_rho_an_mi_an_total = razao_rho_an_mi_an1_lista + razao_rho_an_mi_an2_lista +razao_rho_an_mi_an2_lista + razao_rho_an_mi_an4_lista 
razao_rho_g_mi_g_total = razao_rho_g_mi_g_lista1 + razao_rho_g_mi_g_lista2 + razao_rho_g_mi_g_lista3 + razao_rho_g_mi_g_lista4

mizinho = []
for c in range(len(razao_rho_an_mi_an_total)):
    razao = razao_rho_an_mi_an_total[c]/ (razao_rho_g_mi_g_total[c]*10)
    mizinho.append(razao)

# print(mg_total)
# print('#####################################################################################################################################################')
# print(m_an_total)
# print('#####################################################################################################################################################')
# print(m_ponto_fenda_total)
# print('#####################################################################################################################################################')
# print(razao_rho_an_mi_an_total)
# print('#####################################################################################################################################################')
# print(razao_rho_g_mi_g_total)
# print('#####################################################################################################################################################')
# print(mizinho)
# plt.plot(x_total,mg_total )
# plt.show()
# plt.figure(figsize = (15,10))
# plt.plot(x1*1000, mg_ponto_zr_lista)
# plt.plot(x2*1000, mg_ponto_zp_lista)
# plt.plot(x3*1000, mg_ponto_zs_lista)
# plt.plot(x4*1000, mg_ponto_zd_lista)

# plt.legend()
# plt.show()


#Eq. 50




#Eq. 51


# #Eq. 52


# #Eq. 53
# razao_rho_an_mi_an1 = tc.razao_rho_g_mi_an(m_ponto_g1, a_ft)
# razao_rho_an_mi_an2 = tc.razao_rho_g_mi_an(m_ponto_g2, a_ft)
# razao_rho_an_mi_an3 = tc.razao_rho_g_mi_an(m_ponto_g3, a_ft)
# razao_rho_an_mi_an4 = tc.razao_rho_g_mi_an(m_ponto_g4, a_ft)



#Eq. 55
mi_ar = tc.mi_ar_funcao(T3)
mi_g= tc.mi_g_funcao(T3)

# #Eq. 54 
# t = ((1/2) * 10**-3)
# # n_r = tc.eta_r(0.5, mi_ar, mi_g_funcao, t, s)
# n_r_sem_resfriamento = tc.eta_r_sem_resfriamento(0.6, mi_ar, mi_g)
# n_r_com_refriamento = tc.eta_r_com_resfriamento(0.6, mi_ar, mi_g,t,s)
# print(n_r_sem_resfriamento)


# tgw1_lista = []
# for c in range(len(x1)):
#     t_g_w_1 = tc.t_g_w(tg1_lista[c],n_r_sem_resfriamento,T3)
#     tgw1_lista.append(t_g_w_1)

# tgw2_lista = []
# for c in range(len(x1)):
#     t_g_w_2 = tc.t_g_w(tg2_lista[c],n_r_sem_resfriamento,T3)
#     tgw2_lista.append(t_g_w_2)


# tgw3_lista = []
# for c in range(len(x1)):
#     t_g_w_3 = tc.t_g_w(tg3_lista[c],n_r_sem_resfriamento,T3)
#     tgw3_lista.append(t_g_w_3)

# tgw4_lista = []
# for c in range(len(x1)):
#     t_g_w_4 = tc.t_g_w(tg4_lista[c],n_r_sem_resfriamento,T3)
#     tgw4_lista.append(t_g_w_4)



# #Eq. 74
# m_ponto_h_zd = tc.m_ponto_h_zd_funcao(m_ponto_3,m_ponto_zp,m_ponto_zs,m_ponto_fenda_zd)



#Eq. 75 -78 Areas das orificios respectivos de cada zona (se trata de uma iteração)
#Aqui temos que calcular A_orificios  de cada zona e isso é resultado de uma iteração entre as funcoes.
#Pensei em fazer um while para cada zona, basicamente ao fim dessas iterações o valor arbritado no inicio c_d_h_inicial dve ser proximo ao c_d_h_final
delta_P_1_2 = P2 - P1
c_d_h_inicial = 2 # Valor Abritario para inicio da iteração para encontrar um Cdh que convirja


# while round(c_d_h_inicial)!=round(c_d_h_final):
#     A_h = tc.area_orificio(m_ponto_zp,T3,P3,c_d_h_inicial,delta_P_1_2)
#     alpha = a_an / A_h #Renan: Acho que a relação está invertida 
#     #Abaixo considerei m_ponto_fenda como m_an1
#     m_ponto_an = m_ponto_3 - man_1_lista[0] - m_ponto_h_zp
#     beta = m_ponto_zp / m_ponto_an #Renan: precisamos variar o zp,zs,zd (cada zona i...)? ou só o de entrada zp?
#     mi = beta / alpha
#     k = tc.k_funcao(delta_P_1_2,mi,beta)
#     c_d_h_final = tc.c_d_h(k,delta_P_1_2,beta)



# #Eq. 78 - Coeficiente de descarga
# #Essa equacao utiliza os mesmos parametros da equacao (77) a fim de formar o gráfico 3
# delta=0.8 # ou 0.6 de acordo com o canto do orificio, convexo ou vivo

# c_d_h=tc.cdh(k, 0.8, beta)


# #Eq. 79 - relacao area e diametro dos orificios
# #Dados da tabela 4, precisamos fazer a relação para cada fileira de furos internos e externos i=1,2,3...

# d_hi=tc.relacao_area_diam()

#Equacao 80 ja comentada é apenas a soma dos valores da tabela 4 (numero de orificios)

#Equacao 86 - Area anular externa
# os valores de Dint, Dref e Dft sao os mesmos obtidos e utilizados na tabela 6
#substitui o Din por Dint que está na fórmula do Caio porque não existe o subscrito in de acordo com a lista de subscritos

#Equacao 87 - Area anular interna
#Substituir os mesmos valores da tabela 6 utilizados acima

#Equacao 88 - Utiliza de parametros as equacoes 87 e 86 mas não consegui identificar o que seria o A_h_i (Area_fileira orificio) em valores



# ####### Cálculo da Temperatura na Parede ao Longo do Tubo de Chama
#Posição das Fendas Escolhidas
#Fenda 1
posFenda1= 21*1e-3

id_1 = np.argmin(np.abs(x_total-posFenda1))
vazaoMassicaGases_1 = m_ponto_fenda_total[id_1] 
vazaoMassicaAnular_1 = m_an_total[id_1]
temperaturaGases_1 = tg_total[id_1]

#Posição da Fenda

# Área da Fenda 1
s_1 = 0.001  # altura da fenda [m]
t_1 = 0.002   #espessura da fenda [m]
A_fenda_1 = 2*np.pi*s_1*(d_ref + d_ft)
A_an = a_ref_aerodinamica - a_ft 
# Vazão Mássica de Ar que Entra em cada Fenda
m_fenda_1 = vazaoMassicaAnular_1*(A_fenda_1/A_an)

index = 0
for i in x_total:
    if i> posFenda1:
        [index] = mg[index] + m_fenda_1
        m_an[index] = m_an[index] - m_fenda_1
    index+=1

#Fenda 2
posFenda2= 29*1e-3
id_2 = np.argmin(np.abs(x_t-posFenda2))
vazaoMassicaGases_2 = mg[id_2]
vazaoMassicaAnular_2 = m_an[id_2]
temperaturaGases_2  = Tg[id_2]

# Área da Fenda 2
s_2 = 0.001  # altura da fenda [m]
t_2 = 0.002    #espessura da fenda [m]
A_fenda_2 = 2*np.pi*s_2*(D_ref + D_ft)
# Vazão Mássica de Ar que Entra em cada Fenda
m_fenda_2 = vazaoMassicaAnular_2*(A_fenda_2/A_an)

index = 0
for i in x_t:
    if i> posFenda2:
        mg[index] = mg[index] + m_fenda_2
        m_an[index] = m_an[index] - m_fenda_2
    index+=1

# Eficiência do resfriamento Fenda 1
pu_an_1 = vazaoMassicaAnular_1/A_an # produto entre a densidade e a velocidade do ar que atravessa a área anular
pu_g_1 = vazaoMassicaGases_1/A_ft # produto entre a densidade e velocidade do gás no tubo de chama (na região onde está posicionada a fenda)
visc_ar = (0.03863 + 0.00749*T_3 - (5.8564*(10**(-6))*(T_3**2)) + (2.7769*(10**(-9))*(T_3**3)) - (4.600774*(10**(-13))*(T_3**4)))*(10**(-5)) #viscosidade dinamica do ar
visc_g_1 = (0.03863 + 0.00749*temperaturaGases_1-5.8564*(10**-6)*(temperaturaGases_1**2) + 2.7769*(10**-9)*(temperaturaGases_1**3) - 4.600774*(10**-13)*(temperaturaGases_1**4))*(10**-5) #viscosidade dinamica do gás no interior do tubo de chama
razao_m_1 = pu_an_1/pu_g_1
print(f'Razao m fenda 1: {razao_m_1}')

# Eficiência do resfriamento Fenda 2
pu_an_2 = vazaoMassicaAnular_2/A_an # produto entre a densidade e a velocidade do ar que atravessa a área anular
pu_g_2 = vazaoMassicaGases_2/A_ft # produto entre a densidade e velocidade do gás no tubo de chama (na região onde está posicionada a fenda)
visc_g_2 =(0.03863 + 0.00749*temperaturaGases_2-5.8564*(10**-6)*(temperaturaGases_2**2) + 2.7769*(10**-9)*(temperaturaGases_2**3) - 4.600774*(10**-13)*(temperaturaGases_2**4))*(10**-5) #viscosidade dinamica do gás no interior do tubo de chama
razao_m_2 = pu_an_2/pu_g_2
#razao_m_2 = 0.52
print(f'Razao m fenda 2: {razao_m_2}')

#Temperatura dos gases proximo da parede

####### Cálculo da Temperatura na Parede ao Longo do Tubo de Chama





index = 0
vetorTemp = []
graphEff = []


for indice, posX in enumerate(x_total):
	mdotgas = mg_total[index] 
	mdotan = m_an_total[index]
	tgases = tg_total[index]

	if posX < posFenda1:
		visc_g =(0.03863 + 0.00749* tgases-5.8564*(10**-6)*(tgases**2) + 2.7769*(10**-9)*(tgases**3) - 4.600774*(10**-13)*(tgases**4))*(10**-5) 
		posXEff = posX
		nr = 0
	elif posX >= posFenda1 and posX < posFenda2:
		visc_g = mi_g  	#viscosidade dinamica do gás no interior do tubo de chama
		posXEff = posX
		posRel = posX - posFenda1 + 0.001
		razao_m = mizinho

		if razao_m > 0.5 and razao_m<=1.3:
			nr = 1.10*(razao_m**0.65)*((mi_ar/visc_g)**0.15)*((posRel/s_1)**(-0.2))*((t_1/s_1)**(-0.2)) 
		elif razao_m>1.3 and razao_m<=4:  #Editei aqui
			nr = 1.28*((mi_ar/visc_g)**0.15)*((posRel/s_1)**(-0.2))*((t_1/s_1)**(-0.2))
		else:
			nr = 0
			print("Favor verificar intervalo.")
	else:
		visc_g = mi_g  	#viscosidade dinamica do gás no interior do tubo de chama
		posXEff = posX
		razao_m = razao_m_2
		posRel = posX - posFenda2 + 0.001

		if razao_m>0.5 and razao_m<=1.3:
			nr = 1.10*(razao_m**0.65)*((mi_ar/visc_g)**0.15)*((posRel/s_2)**(-0.2))*((t_2/s_2)**(-0.2)) 
		elif razao_m>1.3 and razao_m<=4:  #Editei aqui
			nr = 1.28*((mi_ar/visc_g)**0.15)*((posRel/s_2)**(-0.2))*((t_2/s_2)**(-0.2))
		else:
			nr = 0
			print("Favor verificar intervalo.")

	graphEff.append(nr)
	Tg_w = tgases - nr*(tgases-T3)

	#TRANSFERENCIA DE CALOR NA PAREDE
	#Resolução será feita utilizando a biblioteca Sympy do Python
	# RADIAÇÃO
	#T_w1 Temperatura na superficie interna tubo de chama
	#T_w2 Temperatura na superficie externa tubo de chama
	T_w1, T_w2 = sp.symbols('T_w1,T_w2')

	k_w = 26 #[W/m.K] condutividade térmica - página 79
	t_w = 0.0005 #liner wall thickness m - retirado do GasTurb pag 325 VALIDAR ESSE VALOR (atualizado do desenho do TG)
	sigma = 5.67*10**(-8)   #constante de Stefan-Boltzmann
	e_w = 0.4 #emissividade do material - página 79 
	q = 0.00398/0.243 #VERIFICAR COM O PROFESSOR????
	l_b_int = d_ft   #comprimento característico do gás - interno
	l_b_ext = 1.2*d_ft #comprimento característico do gás - externo
	e_g_int = 1 - np.exp(-0.290*P3*(q*l_b_int)**0.5*tgases**(-1.5)) #emissividade do gás 
	e_g_ext = 1 - np.exp(-0.290*P3*(q*l_b_ext)**0.5*tgases**(-1.5)) #emissividade do gás 
	Z = 0.4 #para alumínio - página 66 do TG

	#fluxo de calor por radiação do gás (devo usar e_g_int ou e_g_ext)
	R_1 = 0.5*sigma*(1+e_w)*e_g_int*(tgases**1.5)*((tgases**2.5)-(T_w1**2.5)) 
	R_2 = Z*sigma*(T_w2**4 - T3**4)  #fluxo de calor por radiação da parede do tubo de chama         

	#O fluxo de calor através da parede K1-2 é dado por:
	K_1_2 = (k_w/t_w)*(T_w1 - T_w2) #Fluxo de calor através da parede

	#CONVECÇÃO
	#condutividade térmica do gás, função de sua temperatura
	k_g = 5.92657*10**(-4)+9.80957*(10**(-5))*tgases-4.89398*(10**(-8))*(tgases**2)+1.5011410*(10**(-11))*(tgases**3)
	#Número de Reynolds na posição x
	Re_x = (mdotgas*posXEff)/(a_ft*visc_g)  

	if posX <= posFenda1:
		D_L = 4*a_ft/d_ref
		C_1 = 0.017*(k_g/(D_L**0.2))*((mdotgas/(a_ft*visc_g))**0.8)*(tgases-T_w1)
	elif posX >= posFenda1 and posX < posFenda2:
		if razao_m > 0.5 and razao_m <= 1.3:
			C_1 = 0.069*(k_g/posXEff)*(Re_x**0.7)*(Tg_w-T_w1)
		elif razao_m > 1.3 and razao_m < 4.0:
			#Essa equação estava errada no trabalho do rapaz.
			C_1 = 0.10*(k_g/posXEff)*Re_x**0.8*(posXEff/s_1)**(-0.36)*(Tg_w-T_w1)
	else:
		if razao_m > 0.5 and razao_m <= 1.3:
			C_1 = 0.069*(k_g/posXEff)*(Re_x**0.7)*(Tg_w-T_w1)
		elif razao_m > 1.3 and razao_m < 4.0:
			#Essa equação estava errada no trabalho do rapaz.
			C_1 = 0.10*(k_g/posXEff)*Re_x**0.8*(posXEff/s_2)**(-0.36)*(Tg_w-T_w1)
		
	D_an = d_ref - d_ft
	k_a = 5.92657*(10**(-4))+9.80957*(10**(-5))*T3-4.89398*(10**(-8))*(T3**2)+1.5011410*(10**(-11))*(T3**3) 
	#fluxo de calor por convecção da parede externa do tubo de chama
	C_2 = 0.020*(k_a/(D_an**0.2))*((mdotan/(a_an*mi_ar))**0.8)*(T_w2 - T3) 

	eq1 = spy.Eq(R_1+C_1-R_2-C_2,0)
	eq2 = spy.Eq(R_1+C_1-K_1_2,0)
	sol = spy.nsolve((eq1,eq2),(T_w1, T_w2),(1500,800))
	tfinal = min(sol[0],sol[1])
	vetorTemp.append(tfinal)
	print(rf"For iteration {index} the solution is T_w1 = {round(sol[0],2):,.2f}, T_w2 = {round(sol[1],2):,.2f}")
	index+=1
