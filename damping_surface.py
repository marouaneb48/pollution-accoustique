from cmath import *

## l'idée est d'écrire toutes les petites fonctions dans une grande fonction 
##   pour que damping surface soit utilisé comme un module dans d'autres application
## la fonction plot_alpha prend comme paramètre le nom du paramètre.\
##les caractéristiques du matériau sont stockées à l'intérieur plot_alpha 3

## pour voir les résultats, éxécutez le fichier main
## !! Gardez le fichier damping_surface dans le meme dossier que le fichier main pour que le code soit éxécuté 



def plot_alpha(material):
    A=1       
    B=1


    if material == 'BIRCHLT':
        # Birch LT
        phi = 0.529  # porosity
        gamma_p = 7.0 / 5.0
        sigma = 151429.0  # resitivity
        rho_0 = 1.2
        alpha_h = 1.37  # tortuosity
        c_0 = 340.0
    elif material == 'MELAMINE':
        # Melamine foam
        phi = 0.99  # porosity
        gamma_p = 7.0 / 5.0
        sigma = 14000.0  # resitivity
        rho_0 = 1.2
        alpha_h = 1.02  # tortuosity
        c_0 = 340.0
    elif material == 'POLYURETHANE':
        # Polyurethane foam
        phi = 0.98  # porosity
        gamma_p = 7.0 / 5.0
        sigma = 45000.0  # resitivity
        rho_0 = 1.2
        alpha_h = 2.01  # tortuosity
        c_0 = 340.0
    elif material == 'SUTHERLAND':
        # Sutherland's law
        a = 0.555           # certains de ces paramètres ne servent pas pour la suite.
        b = 120.0
        T0 = 524.07
        TC = 20
        TR = 9 / 5 * (273.15 + TC)
        mu0 = 1.827 * 1e-05
        mu = mu0 * (a * T0 + b) / (a * TR + b) * (TR / T0) ** (3 / 2)
        k0 = 1.6853 * 1e-09
        sigma = mu / k0

        phi = 0.7
        gamma_p = 7.0 / 5.0
        rho_0 = 1.2
        alpha_h = 1.15
        c_0 = 340.0


    # parameters of the geometry
    L = 0.01

    # parameters of the mesh
    resolution = 12  # := number of elements along L

 

    eta_0=1
    ksi_0=1/c_0**2
    eta_1=phi/alpha_h
    ksi_1=phi*gamma_p/((c_0**2))
    a=sigma*(phi**2) *gamma_p /((c_0**2)*rho_0*alpha_h)


    def g(k):
        if k==0:
            return 1
        else:
            return 0
        

    def f(x,k,ksi_0,eta_0,omega,L):
        lambda0=lambda_0(k,ksi_0,eta_0,omega)

        return (lambda0*eta_0 -x)*exp(-lambda0*L)\
            +(lambda0*eta_0+x)*exp(lambda0*L)

    def lambda_0(k,ksi_0,eta_0,omega):
        if k**2>=ksi_0*(omega**2)/eta_0:
            return sqrt(k**2 - ksi_0*omega**2/eta_0)
        elif k**2<= ksi_0*omega**2/eta_0:
            return complex(0,sqrt(ksi_0*(omega**2)/eta_0 -k**2))

    def lambda_1(k,ksi_1,eta_1, omega,a):
        return  complex(1/sqrt(2) *sqrt(k**2-ksi_1*omega**2/eta_1\
            +sqrt((k**2-ksi_1*omega**2/eta_1)**2+\
            (a*omega/eta_1)**2)), -1/sqrt(2) \
                *sqrt(ksi_1*omega**2/eta_1 -k**2 +\
                    sqrt((k**2-ksi_1*omega**2/eta_1)**2+\
                (a*omega/eta_1)**2)))

    def chi(k,alpha,ksi_0,eta_0,omega,ksi_1,eta_1,a):  
        lambda0=lambda_0(k,ksi_0,eta_0,omega)       ## pour réduire les expressions des lambda, on les calcules on les stocke
        lambda1=lambda_1(k,ksi_1,eta_1, omega,a)

        x=g(k)*((lambda0*eta_0-lambda1*eta_1)/f(lambda1*eta_1,k,ksi_0,eta_0,omega,L) -\
            (lambda0*eta_0-alpha)/f(alpha,k,ksi_0,eta_0,omega,L))
        return x

    def gamma(k,alpha,ksi_0,eta_0,omega,ksi_1,eta_1,a):
        lambda0=lambda_0(k,ksi_0,eta_0,omega)
        lambda1=lambda_1(k,ksi_1,eta_1, omega,a)
        y=g(k)*((lambda0*eta_0+lambda1*eta_1)/f(lambda1*eta_1,k,ksi_0,eta_0,omega,L) -\
            (lambda0*eta_0+alpha)/f(alpha,k,ksi_0,eta_0,omega,L))
        return y

    def e(k,alpha,ksi_0, eta_0,omega,ksi_1,eta_1,a,L):
        
        lambda0=lambda_0(k,ksi_0,eta_0,omega)
        chi0=chi(k,alpha,ksi_0,eta_0,omega,ksi_1,eta_1,a)
        gamma0=gamma(k,alpha,ksi_0,eta_0,omega,ksi_1,eta_1,a)

        if k**2 >=ksi_0*(omega**2)/eta_0:    

            return (A+B*(abs(k)**2))*((1/(2*lambda0) *((abs(chi0)**2)*(1-exp(-2*lambda0*L)))+abs(gamma0)**2*(exp(2*lambda0*L)-1))\
                +2*L*chi0*(gamma0.conjugate()).real)+ B*(lambda0/2)*((abs(chi0)**2)*(1-exp(-2*lambda0*L))+\
                    (abs(gamma0)**2)*(exp(2*lambda0*L)-1))-2*B*(lambda0**2)*L*(chi0*(gamma0.conjugate())).real
        
        elif k**2<ksi_0*(omega**2)/eta_0:
            return (A+B*abs(k)**2)*(L*((abs(chi0)**2)+(abs(gamma0)**2))+complex(0,(1/lambda0)*(chi0*(gamma0.conjugate())*(1-exp(-2*lambda0*L))).imag))+\
                B*L*(abs(lambda0)**2)*((abs(chi0)**2)+(abs(gamma0)**2))+\
                    complex(0,B*lambda0*(chi0*(gamma0.conjugate())*(1-exp(-2*lambda0*L))).imag) 

    #la fonction erreur prend comme paramètre alphac qui est un couple qu'on transformera plus tard en complexe 
    def erreur(alphac,ksi_0, eta_0,omega,ksi_1,eta_1,a,L):    
        somme=0
        liste=[i*pi/L for i in range(-resolution,resolution+1)]
        for k in liste :
            somme+=e(k,complex(alphac[0],alphac[1]),ksi_0, eta_0,omega,ksi_1,eta_1,a,L)       
        return somme 



    def reel_erreur(alphac,omega):  ## ici on prend la partie réelle de l'erreur 
        
        return erreur(alphac,ksi_0, eta_0,omega,ksi_1,eta_1,a,L).real








    from scipy.optimize import minimize

    import matplotlib.pyplot as plt


    omega_list=[k*pi for k in range(2,30000,100)]   #on prend un pas de 100
    alphac_list_reel=[]
    alphac_list_imag=[]
    erreur_reel_list=[]


    for w in omega_list:
        z=minimize(reel_erreur,[30,-35],args=w,method='CG',options={'gtol': 1e-4})   #on minimise l'erreur et on prend l'argument du min
        alphac_list_reel.append(z.x[0])  ## on stocke la partie réelle dans la liste alphac_list_reel
        
        alphac_list_imag.append(z.x[1])   ## on stocke la partie réelle dans la liste alphac_list_imag
        erreur_reel_list.append(reel_erreur([z.x[0],z.x[1]],w))  ## on stocke l'erreur dans la liste erreur_reel_list
        
     ## et on trace les courbes 
    plt.figure()
        
    plt.plot(omega_list,alphac_list_reel)
    plt.xlabel('omega')
    plt.ylabel('partie réelle de alpha')
    plt.title(material)
    plt.figure()
    plt.plot(omega_list,alphac_list_imag,'tab:orange')
    plt.xlabel('omega')
    plt.ylabel('partie imaginaire de alpha')
    plt.title(material)
    plt.figure()
    plt.plot(omega_list,erreur_reel_list)
    plt.xlabel('omega')
    plt.ylabel('erreur')
    plt.title(material)
    plt.show()

    





