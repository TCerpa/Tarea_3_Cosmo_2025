import numpy as np
import camb
import matplotlib.pyplot as plt

#Set cosmology
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.66, ombh2=0.02242, omch2=0.11933)
#Configuracion parametros de espectro primordial n_s , A_s con valores por defecto de CAMB
pars.InitPower.set_params()
#Calcula el espectro de potencia P(k)
pars.set_matter_power(kmax=10)
#Obtenemos los resultados de la cosmología usada, el espectro primordial y el espectro de potencia en un objeto CAMBresults
results = camb.get_results(pars)

#Guardamos los resultados: las funciones transferencia calculadas por CAMB. T_a(k,z), con a un sufijo para las distintas especies.
f_transfer = results.get_matter_transfer_data()
#Extraemos los valores de k/h usados en el calculo de T_a(k,z). 
kh = f_transfer.transfer_data[0, :, 0]#esta notacion extrae la segunda dimension de un array en 3D. Segun CAMB estos son los valores k/h.
#print(kh)
#Despejamos los valores de k.
k = kh * results.Params.h
#Tomamos la funcion transferencia total de materia para todos los k con z=0
transfer = f_transfer.transfer_data[6, :, 0]  
#Genera espectro de potencia primordial como una función de k
primordial_PK = results.Params.scalar_power(k)
#Hacemos el calculo explicito del matter power spectrum según lo definido por CAMB
matter_power = primordial_PK * transfer**2 * k**4 / (k**3 / (2 * np.pi**2))

plt.loglog(kh, matter_power)
plt.xlabel(r'$k[h Mpc^{-1}]$')
plt.title('Matter power spectrum')

plt.show()
