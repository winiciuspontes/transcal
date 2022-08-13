from sqlalchemy import false
import sympy  as sp
import numpy as np 
import pandas as pd 
import math as mt 
import matplotlib.pyplot as plt
import transcal as tc

#CONSTANTES UTILIZADAS:
k = 143.5
m_ponto_3 = 0.243
m_ponto_comb = 0.252 #Tirar dúvida sobre o m_comb (nao sabemos qual é o valor)
m_ponto_zp = m_ponto_3 # Não temos certeza sobre esse valor, precisamos tirar dúvida com o professor 
T3 = 498.42 
P3 = 354.420 * 1000 
delta_P_3_4 = (354.420 - 333.155) * 1000 
theta = 73 * (10**6)
d_int = 0.024 # Valor obtido da tabela 4



#CALCULO DOS PHIS:
phi_zp = 1.249
phi_zs = 0.8
phi_zd = 0.237
phi_global = 0.237
phi_pobre = 0.45
phi_rico = 2.553


#Calculando o delta T_phi com phi = 1 (Eq.39)
delta_T_phi = 2185 - (0.5 * T3)
delta_T_zp = 2830 - (800 * phi_zp) #Aproximado para quando T3=400K  - Não sabemos que de onde vem esse valor de 2830
delta_T_zs = 1600 - (0.5 * T3)


b = tc.b_function(phi_zp) #Eq. antes da 30
d_ref = 0.03863 # Valor obtido da tabela 4

# a_ref_aerodinamico = tc.area_tranferencia_aerodinamica(k, m_ponto_3, T3, P3)
# d_ref = tc.altura_referencia(a_ref_aerodinamico, d_int)


a_ref = 0.00759996 # A_ref do Caio

d_ref = tc.altura_referencia(a_ref, d_int) # Eq. Tabela 6 

a_ft = tc.area_ref_trans(a_ref) #Eq. 34 - No trabalho do Caio ele multiplicou por 10 ^ 6

d_ft = tc.altura_tubo_chama(a_ft, d_int, d_ref) #Eq. Tabela 6 (para o caso anular)


#Calculando o Lcc com os parametros necessarios:
l_zp = (3/4) * d_ft #Enunciado antes da Eq.35 - Pegamos o l_zp_max, para pegar o minimo trocamos o 3/4 por 2/3
l_zs = (1/2) * d_ft #Enunciado antes da Eq.35
l_zd = (3/2) * d_ft #Enunciado antes da Eq.35
l_zr = (1/2) * d_ft #Enunciado antes da Eq.35
lcc = tc.comprimento_camara_combustao(l_zp,l_zs,l_zd) #Eq.35 - Nosso resultado bate com a Tabela 11

#Eq. 36 (phi_global_rico/phi_zs) é convencionada da literatura como sendo 0.8 

#Calculando a Eq.37
m_ponto_arref = tc.porcentagem_ar_resfriamento(T3, m_ponto_3) #Eq.37 - Equacao nao confirmada(não achamos o ela no Tg)

#Calculando a Eq.38 - Não achamos o valor de m_zp e m_zs. Conseguiriamos calcular o m_zd.

#Eq. 40 - calculando o eta
eta_zr = tc.eta_zr(T3,P3)
#Eq.39 - As eq estavam fora de ordem
t_max_zr = tc.t_max_zr(T3, eta_zr, delta_T_phi)

#Eq.41 - Calculando a T_med - dando dif de 10 graus :
t_media_zr = tc.t_media_zr(T3, t_max_zr)


#Eq.43 - Calculando o eta_zp:

eta_zp = tc.eta_zp(T3, P3)

#Eq.42 - Calculando a Tsaida_zp - Caio trocou a ordem novamente:
t_saida_zp = tc.t_saida_zp(T3, eta_zp, delta_T_zp)

#Eq.45 - Calculando o Tsaida_zs rico
eta_zs_rico = tc.eta_zs(phi_zs)

#Eq. 44 - Calculando a Tsaida_zs rico - Caio trocou a ordem das eq novamente:
t_saida_zs_rico = tc.t_saida_zs(T3, eta_zs_rico, delta_T_zs)  


#Eq. 49 - Calculando o V_zs:
v_zs = tc.v_zs(a_ft, l_zs)


# Eq. 48 - Calculando o psi_maiusculo:
#Vamos calcular o psi_t3 porque o psi_300 é um caso especifico. Perg para o professor
psi_t3 = tc.phi_t3(phi_zs, m_ponto_comb, v_zs, P3, T3) 
# psi_300 = tc.phi_t300(m_ponto_comb, P3, phi_zs,  v_zs) 


#Eq. 47 - Calculando o D*:
d_asterisco = tc.d_ast(delta_P_3_4,P3)


#Eq.46 - Calculando o eta_zs pobre(caio trocou a ordem novamente):
eta_zs_pobre = tc.eta_zs_pobre(psi_t3, phi_zs, d_asterisco)

#Eq. 44 para Tsaida_zs pobre - Caio trocou a ordem das eq novamente:
t_saida_zs_pobre= tc.t_saida_zs(T3, eta_zs_pobre, delta_T_zs) 

#Eq. veio da tabela 12 
t_saida_zd = 1134

#Tabela 7

# Tg1 = t_media_zr 
# print(Tg1)
# Tg2 = tc.tg(t_media_zr, t_saida_zp, 0, 0, l_zp, l_zr, 0, 0, 0.01622, 2)
# print(Tg2)
# Tg3 = tc.tg(0, t_saida_zp, t_saida_zs_rico, 0, l_zp,0, l_zs, 0, 0.02704, 3)
# print(Tg3)
# Tg4 = tc.tg(0, 0, t_saida_zs_rico,t_saida_zd, l_zp , 0, l_zs, l_zd, 0.05945, 4)
# print(Tg4)

tg_lista = []

x1 = np.linspace(0, l_zr, 30)
x2 = np.linspace(l_zr, 0.01622, 30)
x3 = np.linspace(0.01622, 0.02704, 30)
x4 = np.linspace(0.02704, 0.05945, 30)
x_total = np.linspace(0, 0.05945, 120)


#Fazendo o looping para conseguir plotar o gráficos das temperaturas
for c in range(len(x1)):
    Tg1 = t_media_zr
    tg_lista.append(Tg1)

for d in range(len(x2)):
    Tg2 = tc.tg(t_media_zr, t_saida_zp, 0, 0, l_zp, l_zr, 0, 0, x2[d], 2)
    tg_lista.append(Tg2)

for e in range(len(x3)):
    Tg3 = tc.tg(0, t_saida_zp, t_saida_zs_rico, 0, l_zp,0, l_zs, 0, x3[e], 3)
    tg_lista.append(Tg3)

for f in range(len(x4)):
    Tg4 = tc.tg(0, 0, t_saida_zs_rico,t_saida_zd, l_zp , 0, l_zs, l_zd, x4[f], 4)
    tg_lista.append(Tg4)

# plt.scatter(x_total,tg_lista)
# plt.show()


m_ponto_g1 = tc.mg_ponto(m_ponto_zp, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
m_ponto_g2 = tc.mg_ponto(m_ponto_zp, m_ponto_g1, 0, 0, 0, 0, l_zr, l_zp, 0, 0, 0.01622,  2)
m_ponto_g3 = tc.mg_ponto(0, 0, m_ponto_3, 0, m_ponto_g2 , 0, 0, l_zp, l_zs, 0, 0.02704,  3)
m_ponto_g4 = tc.mg_ponto(0, 0, 0, m_ponto_g3, 0 , m_ponto_3, 0, l_zp, l_zs, l_zd, 0.05945,  4)

print(m_ponto_g1, m_ponto_g2, m_ponto_g3, m_ponto_g4)


#Eq. 50
area_fenda = tc.area_fenda(d_ref, d_ft) # A altura da fenda é a soma das alturas de referencia e do tubo de chama 


#Eq.51 (confirmar com o grupo se vamos colocar a eq 51 como sendo um array igual o de cima)
m_an1 = m_ponto_3 - m_ponto_g1
m_an2 = m_ponto_3 - m_ponto_g2
m_an3 = m_ponto_3 - m_ponto_g3
m_an4 = m_ponto_3 - m_ponto_g4

a_an = a_ft/0.65
m_ponto_fenda1 = tc.m_ponto_fenda_funcao(area_fenda, m_an1, a_an)
m_ponto_fenda2 = tc.m_ponto_fenda_funcao(area_fenda, m_an2, a_an)
m_ponto_fenda3 = tc.m_ponto_fenda_funcao(area_fenda, m_an3, a_an)
m_ponto_fenda4 = tc.m_ponto_fenda_funcao(area_fenda, m_an4, a_an)
