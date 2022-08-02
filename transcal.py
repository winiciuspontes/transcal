
#importanto as bibliotecas necessárias
#from tkinter.tix import Tree
import sympy  as sp
import numpy as np 
import pandas as pd 
import math as mt 
import matplotlib as plt
from sympy import beta





#Todos os simbólicos usados
x, y, z, phi, T, T_max, T1, T2, T4,T_mr,delta_T_phi,delta_T_zp,delta_T_zs, delta,P,  P1, P2, P3, A_ref, A_h_eff, phi_ZP, theta, q_ref, A_h_reff, delta_P_3_4, R, m_ponto_3, T3, k, D_ref, phi_global, m_ponto_zp, phi_pobre, phi_rico,A_ft, D_ref, L_cc, L_zp, L_zs, L_zd,L_sz,L_zr, L_dz,  phi_Zs, m_ponto_arref, m_ponto_zd, eta_zr, eta_zp, eta_zs, p_3,T_med_zr, T_max_zr, delta_p, psi_t3, m_ponto_comb,V_zs, teste, T_saida_zp, T_saida_zs, T_saida_zd, mg_ponto_zr, mg_ponto_zp, mg_ponto_zs, m_ponto_zs, rho_an,u_an,rho_g, u_g, m_ponto_fenda, A_fenda, m_ponto_g, mi_ar, mi_g, T_g,  eta_r,D_ft, m_ponto_fenda_zp, m_ponto_fenda_zp, m_ponto_fenda_zs, m_ponto_fenda_zd,mi, beta,N_h_i_ext,N_h_i_int  = sp.symbols(
    [
    'x','y','z', 'phi', 'T', 'T_max', 'T1', 'T2', 'T4', 'T_mr', 'delta_T_phi', 'delta_T_zp','delta_T_zs', 'delta','P', 'P1',
    'P2', 'P3', 'A_ref', 'A_h_eff', 'phi_ZP', 'theta', 'q_ref', 'A_h_reff', 
    'delta_P_3_4','R', 'm_ponto_3','T3', 'k', 'D_ref', 'phi_global', 'm_ponto_zp', 'phi_pobre', 
    'phi_rico','A_ft','D_ref','L_cc','L_zp','L_zs','L_zd', 'L_sz', 'L_zr','L_dz',  
    'phi_Zs','m_ponto_arref', 'm_ponto_zd','eta_zr','eta_zp','eta_zs', 'p_3', 'T_med_zr', 
    'T_max_zr','delta_p', 'psi_t3', 'm_ponto_comb', 'V_zs', 'teste', 'T_saida_zp', 
    'T_saida_zs', 'T_saida_zd', 'mg_ponto_zr', 'mg_ponto_zp', 'mg_ponto_zs', 'm_ponto_zs', 
    'rho_an','u_an','rho_g', 'u_g', 'm_ponto_fenda', 'A_fenda','m_ponto_g', 'mi_ar', 'mi_g', 'T_g',  'eta_r', 'D_ft', 'm_ponto_fenda_zp', 'm_ponto_fenda_zd', 'm_ponto_fenda_zs', 'm_ponto_fenda_zd', 'mi', 'beta','N_h_i_ext','N_h_i_int'
    ])



#Equacao 5
def razao_ar_comb_esteq(a, PMar, PMcomb):
    A_C_esteq = a*(1 + 3.7274 + 0.0458) * (PMar/PMcomb)
    return A_C_esteq


#Equacao 8 - Define se a mistura é rica ou pobre
def razao_equivalencia(A_C_esteq, A_C_real):
    phi = A_C_esteq/A_C_real
    return ('rica' if phi > 1 else 'pobre')


#Equacao 19 - Perda de pressao quente

def perda_pressao_quente(rho, u, T4, T3):
    delta_P_quente = 0.5 * rho * (u**2) * ((T4/T3) - 1 )
    return delta_P_quente


#Equacao 20 - perda de pressao fria 

def perda_pressao_fria(P_difusor, P_turbo_chama):
    delta_P_fria = sum(P_difusor, P_turbo_chama)
    return delta_P_fria


#Equacao 21 - perda de pressao do difusor 
#Sendo 𝑃1 e 𝑃2 as pressões totais nas seções transversais logo ao início e logo ao final do difusor, respectivamente
delta_P_difusor = P1 - P2  


#Equacao 22 - perda de pressao do tubo de chama 

def P_turboChama(A_ref, A_h_eff, q_ref):   

    delta_P_turbo_chama = ((A_ref/A_h_eff) ** 2) * q_ref

    return delta_P_turbo_chama

#Equacao 23 - dimensionar os orifícios do tubo de chama 

def dimensao_orificio_tubo_chama(A_ref, delta_p_3_4, q_ref, delta_P_difusor):
    dimensao_orificioTuboChama = A_ref/(((delta_p_3_4/q_ref) - (delta_P_difusor/q_ref))**0.5)
    return dimensao_orificioTuboChama


#Equacao 24 - dimensionar os orifícios do tubo de chama
#C_d e A_h serão arrays(listas)
def A_h_eff_function(C_d, A_h):
    somatorio_A_h_reff = 0
    for c in range(len(C_d)):
        somatorio_A_h_reff =  somatorio_A_h_reff + C_d[c] * A_h[c]  
        
    return somatorio_A_h_reff


#Equacao 26 - qualidade transversal de Temperatura 

def qualidade_tranversal_temp(T_max, T4, T3):
    qualidade_trans_temp = (T_max - T4)/(T4 - T3)
    return qualidade_trans_temp


#Equacao 27 - Perfil radial de Temperatura 

def perfil_radial_temp(T_mr, T4, T3):
    perfil_rad_temp = (T_mr - T4)/(T4 - T3)
    return perfil_rad_temp

#Equacao 28 - Obtenção da área de referência pela aerodinâmica 

def area_tranferencia_aerodinamica(k, m_ponto_3, T3, P3):
    camara_analunar = True
    if(camara_analunar):
        razao_delta_P_3_4_q_ref = 20/0.06 # Valor retirado da Tabela 5 
        area_transf =  ( (k * (((m_ponto_3 * (T3**0.5) )/(P3))**2  * (razao_delta_P_3_4_q_ref)) ) ** 0.5 ) 
    return area_transf


#Equacao 29- Obtenção da área de referência pela química 

def area_tranferencia_quimica(P3, D_ref, T3, b, m_ponto_3, theta):
    A_ref =  (theta * m_ponto_3) / ((P3**1.75)  * (D_ref ** 0.75) * np.exp(T3/b))
    return A_ref





#Equacao 30 
def phi_zp_function(phi_global, m_ponto_zp, m_ponto_3):
    phi_ZP = phi_global/(m_ponto_zp/m_ponto_3)
    return phi_ZP


#Equacao 31
def phi_global_function(m_comb_op, m_ponto_3_op, m_comb, m_ponto_3):
    phi_global = (m_comb_op/m_ponto_3_op)/ (m_comb/m_ponto_3)
    return phi_global


#Equacao do b, antes da 30

def b_function(phi_ZP):
    b_calculo = 170 * (2 - np.log(phi_ZP))
    return b_calculo

#Equacao 32 

def phi_pobre_function(T3):
    phi_pobre = 0.70547 - 0.00046 * T3
    return phi_pobre

#Equacao 33

def phi_rico_function(T3):
    phi_rico = 1.46695 + 0.00172 * T3
    return phi_rico

#Equacao 34

def area_ref_trans(A_ref):
    A_ft = 0.56*A_ref
    return A_ft

#Tabela 6 

def altura_referencia(A_ref,D_int):
    D_ref = ((np.sqrt( ((4 * A_ref)/np.pi) + (D_int**2))) - D_int)/2
    return D_ref

def altura_tubo_chama(A_ft,D_int,D_ref):
    D_ft = (A_ft/(np.pi()*(D_int+D_ref)))
    return D_ft


#Equacao 35 - determinação da área de combustao

def comprimento_camara_combustao():
  L_cc = L_zp + L_zs + L_zd
  return L_cc

#Equacao 36 - No entanto, talvez tenha valor fixo de 0.8
def quant_ar_zona_secundaria(phi_global,rico,phi_Zs):
  quanti_relativa_ar = (phi_global+rico)/phi_Zs
  return quanti_relativa_ar


#Equacao 37 

def porcentagem_ar_resfriamento(m_ponto_arref, T3, m_ponto_3):
    m_ponto_arref = (0.1*T3 - 30) * m_ponto_3
    return m_ponto_arref


#Equacao 38 

def vazao_ar_zona_diluicao(m_ponto_zp, m_ponto_zs, m_ponto_arref, m_ponto_3):
    m_ponto_zd = (1 - (sum(m_ponto_zp, m_ponto_zs, m_ponto_arref))/m_ponto_3 )
    return m_ponto_zd


#Equacao 39 
def t_max_zr(T3, eta_zr, delta_T_phi):
    t_max_zr = T3 + eta_zr*delta_T_phi
    return t_max_zr

#Equacao 40 ele não estava aceitando os valores como simbolico 
def eta_zr(T3,p_3):
    variavel_dentro_tanh = 1.5475 * (10**-3) * (T3 + 108 * np.log(p_3) - 1863)
    return (0.92 + (0.12 * mt.tanh(variavel_dentro_tanh)))


#Equacao 41

def t_media_zr(T3, T_max_zr):
    t_media_calculo = (1/3 * T3) + (2/3 * T_max_zr)
    return t_media_calculo

#Equacao 42 

def t_saida_zp(T3, eta_zs, delta_T_zs):
    T_saida_zs = T3 + (eta_zs * delta_T_zs) 
    return T_saida_zs

#Equacao 43 tem o mesmo problema da equacao 40 ]
def eta_zp(T3,p_3):
    variavel_dentro_tanh = 1.5475 * (10**-3) * (T3 + 108 * np.log(p_3) - 1863)
    return (0.92 + (0.12 * mt.tanh(variavel_dentro_tanh)))
#Trocar pelo o que é requerido


#Equacao 44 

def t_saida_zs(T3, eta_zs, delta_T_zs):
    T_saida_zs = T3 + (eta_zs * delta_T_zs) 
    return T_saida_zs

#Equacao 45 
def eta_zs_rico(phi_Zs):
    return 1/phi_Zs


#Equacao 46

def eta_zs_pobre(psi_300,phi_Zs,d_ast):
    0.911*mt.log(psi_300)+8.02*phi_Zs-1.097+d_ast
    return 1/()

#Equacao 47 

def d_ast(delta_p, P3):

    d_ast_calculo = 0.736 - 0.0173 * P3/ delta_p
    return d_ast_calculo

#Equacao 48 

def phi_t3(phi_Zs, m_ponto_comb, V_zs, P):
    psi_t3 = (10 ** (-3.054 * (phi_Zs **-1.25))) * (T3 ** (1.2327 * (phi_Zs ** -1.205))) * m_ponto_comb/ V_zs * (P**(2*phi_Zs))
    return psi_t3


#Equacao 49 

def V_zs(A_ft, L_sz):
    V_zs = A_ft * L_sz
    return V_zs


#Tabela 8 - Definicao do Tg - conferir se nao esqueceu nenhum termo

def tg(T_med_zr, T_saida_zp,T_saida_zs,T_saida_zd, L_zp, L_zr, L_zs, L_dz, x, n_regiao):
        
    if(n_regiao == 1):
        Tg = T_med_zr

    elif(n_regiao == 2):
        Tg = T_med_zr + ((T_saida_zp + T_med_zr)/(L_zp - L_zr) * (x - L_zr))

    elif(n_regiao == 3):
        Tg = T_saida_zp + ((T_saida_zs - T_saida_zp)/L_zs * (x - L_zp))

    elif(n_regiao == 4):
        Tg = T_saida_zp + ((T_saida_zd - T_saida_zs)/ L_dz * (x - L_zp - L_zs))
    
    return Tg


#Tabela 8 - Definicao da vazao massica dos gases localmente

def mg_ponto(m_ponto_zp, mg_ponto_zr, m_ponto_zs, m_ponto_zd, L_zr, L_zp, L_zs, L_zd, n_regiao_massa):



    if(n_regiao_massa == 1):
        mg_ponto = 3/4 * m_ponto_zp

    elif(n_regiao_massa == 2):
        mg_ponto = mg_ponto_zr + ((m_ponto_zp - mg_ponto_zr) * (x - L_zr)/(L_zp - L_zr) )

    elif(n_regiao_massa == 3):
        mg_ponto = mg_ponto_zp + ((m_ponto_zs - mg_ponto_zp) * (x - L_zp)/L_zs)

    elif(n_regiao_massa == 4):
        mg_ponto = mg_ponto_zs + ((m_ponto_zd - mg_ponto_zs) * (x - (L_zp + L_zs))/L_zd)
    
    return mg_ponto


#Equacao 50

def area_fenda(D_ref, D_ft, s):
    area_fenda_calculo = (2 * np.pi * s) * (D_ref + D_ft)
    return area_fenda_calculo

#Equacao 51

def m_ponto_fenda_funcao(area_fenda, m_an, A_an):
    vazao_massica_fenda_calculo = m_an * (area_fenda/A_an)
    return vazao_massica_fenda_calculo


#Da equacao 52 ate 54 a relação é simples e não compensa fazer uma função 


#Equacao 55

def mi_ar_funcao(T3):
    mi_ar_calculo = (0.03863 + (0.00749 * 𝑇3) - ((5.8564 * 10**-6) * T3 ** 2) + ((2.7769 * 10**-9) * T3**3) - ((4.600774 * 10 **-13) * T3**4 )) * 10 **-5
    return mi_ar_calculo


#Equacao 56

def mi_g_funcao(T_g):
    mi_g_calculo = (0.03863 + (0.00749 * T_g) - ((5.8564 * 10**-6) * T_g ** 2) + ((2.7769 * 10**-9) * T_g**3) - ((4.600774 * 10 **-13) * T_g**4 )) * 10 **-5
    return mi_g_calculo


#Da equacao 57 ate 60 eu não entendi muito bem 


#Equacao 61

def r1_funcao(sigma, epsilon_w, epsilon_g,  T_g, T_w1):
    R1_calculo = 0.5 * sigma * (1 + epsilon_w) * epsilon_g * (T_g**1.5) * ((T_g ** 2.5) - (T_w1 ** 2.5))
    return R1_calculo


#Equacao 63 

def r2_funcao(Z, sigma, T_w2, T3):
    r2_calculo = Z * sigma * ((T_w2 ** 4) - (T3**4))
    return r2_calculo

#Equacao 64

def c1_funcao_0_5(k_g, x, Re, T_wg, T_w1):
    c1_calculo = 0.069 * (k_g/x) * (Re ** 0.7 ) * (T_wg - T_w1)
    return c1_calculo

#Equacao 65

def c1_funcao_1_3(k_g, x, Re, T_wg, T_w1, s):
    c1_calculo = 0.010 * (k_g/x)* (Re ** 0.8 ) * ((x/s) **-0.36) * (T_wg - T_w1)
    return c1_calculo

#Equacao 66
def kg_funcao(T_g):
    k_g_calculo = (5.92657 * 10**4) +( 9.80957 * 10-5 ** T_g) - ((4.89398 * (10**-8) * (T_g**2))) + ((1.501141010 **-11) * (T_g**3))
    return k_g_calculo


#Equacao 67 - Mas ele usa uma outra relacao também que é a eq 68

def reynolds(rho, u, x, mi):
    Re = (rho * u * x)/mi 
    return Re

#Equacao 69 

def c2_funcao(Ka, D_an, m_ponto_an, A_an, mi_a, T_w2, T3):
    c2_calculo = 0.02 * (Ka/(D_an ** 0.2) * ((m_ponto_an)/(A_an*mi_a) ** 0.8)) * (T_w2 - T3)
    return c2_calculo

#Da equacao 69 ate 74 as equacoes são simples (geralmente de subtração) então não vale a pena 
#fazer uma funcao 


#Equacao 75 

def area_orificio(m_ponto_hi, T3, P3, C_dh, delta_P_h):
    a_hi = ((143.5 * m_ponto_hi**2 * T3)/(P3**2 * C_dh**2 * (delta_P_h/P3))) ** 0.5
    return a_hi

#Equacao 76 é uma somatoria de valores, portanto, acho melhor fazer sem função 



#Equacao 77 

def k_funcao(delta, mi, beta):
    k_calculo = 1 + delta**2 * (((2* mi**2) + (4*mi**2 + ((mi**2)/(delta**2) * (4*beta - beta**2)) )) ** 0.5)
    return k_calculo


#Equacao 78 

def c_d_h(K, delta, beta):
    c_d_h_calculo = (K - 1)/(delta * ((4*K**2) - ((K*(2 - beta)**2) ) ** 0.5) )
    return c_d_h_calculo

#Equacao 79

def d_hi(A_hi, N_hi):
    d_hi_calculo = 2 * np.sqrt(A_hi/(np.pi * N_hi) )
    return d_hi_calculo


#Equacao 80 - é uma relacao de soma, portanto, não vale a pena criar uma função para isso
#Bruno: Estou implemetando essa func de soma, pois ja vou fazer a verificação da validade de Nh_i_ext e Nh_i_int, assim na rotina precisamos apenas invocar essa func
def N_h_i(N_h_i_int,N_h_i_ext,D_int, D_ref, D_ft,d_hi):
    razao_circuferencia = (np.pi*((D_int+D_ref+D_ft))/(2*d_hi))
    if verifica_N_ext(N_h_i_ext,razao_circuferencia)==True:
        print("N_ext Ok!")
    else:
        print("N_ext maior que a razao de circuferencia")
    if verifica_N_int(N_h_i_int,razao_circuferencia)==True:
        print("N_int Ok!")
    else:
        print("N_int maior que a razao de circuferencia")
        


#Equacao 81
def verifica_N_ext(N_h_i_ext,razao_circuferencia):   
    if N_h_i_ext < razao_circuferencia:
        return True
    else:
        return False
    

#Equacao 82
def verifica_N_int(N_h_i_int,razao_circuferencia):
    if N_h_i_int < razao_circuferencia:
        return True
    else:
        return False
    
# Equacao 83 - é uma relacao de soma, portanto, não vale a pena criar uma função para isso
    
# Equacao 84
def A_h_i_ext(A_an_ext,A_ext,A_int,A_h_i):
    return ((A_an_ext)/(A_ext+A_int))*A_h_i

# Equacao 85
def A_h_i_ext(A_an_ext,A_an_int,A_h_i):
    return ((A_an_ext)/(A_an_ext+A_an_int))*A_h_i

# Equacao 86
def  A_an_ext(D_in,D_ref,D_ft):
    return (np.pi/4)*(pow((D_in+(2*D_ref)),2)-pow((D_in+D_ref+D_ft),2))

# Equacao 87
def  A_an_int(D_in,D_ref,D_ft):
    return (np.pi/4)*(pow((D_in+D_ref-D_ft),2)-pow(D_in,2))

#Terminar até a equacação 87 (Ver com o grupo se vão querer fazer com funções)

# if __name__ == "__main__":
#     # Iniciar metodologia
#     print("Iniciando...")
   
