from sqlalchemy import false
import sympy  as sp
import numpy as np 
import pandas as pd 
import math as mt 
import matplotlib as plt
import transcal as tc


k = 143.5
massa_ponto_3 = 0.243    
T3 = 498.42 
P3 = 354.420 * 1000 
delta_P_3_4 = (354.420 - 333.155) * 1000 
theta = 73 * (10**6)
d_int = 0.024 # Valor obtido da tabela 4
b = tc.b_function(1.249)
d_ref = 0.039 # Valor obtido da tabela 4

# a_ref_aerodinamico = tc.area_tranferencia_aerodinamica(k, massa_ponto_3, T3, P3)
# d_ref = tc.altura_referencia(a_ref_aerodinamico, d_int)


a_ref = 0.00759996 # A_ref do Caio
d_ref = tc.altura_referencia(a_ref, d_int)
print(d_ref)


