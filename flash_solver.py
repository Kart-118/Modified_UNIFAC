import moduni as md 
import numpy as np 
import scipy.optimize as scopt 

def Vp(comp :md.Compound,T):
    '''
    Docstring for Psat
    
    :param comp: Array of compounds 
    :type comp: md.Compound
    :param T: Temperature in K
    '''
    vp = []

    for element in comp: 
        val = element.Vp(T)
        vp.append(val) 

    vp = np.array(vp) 

    return vp 

def P_bubble(comp :md.Compound, T , X,params,bool = False):
    '''
    Docstring for P_bubble
    
    :param comp: Array of Compounds 
    :type comp: md.Compound
    :param T: Temperature in K
    :param X: Array of liquid phase mole fractions 
    :param params: UNIFAC parameters generated explicitly using md.generate_params
    ''' 
    X = np.array(X) 

    Psat = Vp(comp,T) 
    gammas = md.gammas_moduni(T,X,params) 

    P = X*gammas*Psat 
    P = np.sum(P) 

    if bool:
        Y = X*gammas*Psat/P 
        return P,Y
    
    else:
        return P

def P_dew(comp: md.Compound,T,Y,params,bool = False): 

    Y = np.array(Y) 

    Psat = Vp(comp,T) 
    gammas = md.gammas_moduni(T,Y,params) 
    val = Y/gammas/Psat
    val = np.sum(val) 
    P = 1/val 

    while True:
        X = P*Y/gammas/Psat
        gammas = md.gammas_moduni(T,X,params)
        val = Y/gammas/Psat
        val = np.sum(val) 
        Pnew = 1/val  

        err = np.abs(Pnew - P)/Pnew

        if err<1e-06:
            break 

        else:
            P = Pnew 
        
    if bool:
        return P,X 
    else:
        return P

def VLE_flash_PT(comp :md.Compound,P,T,Z,params): 
    P = P*1e+05 
    Z = np.array(Z) 

    Pb = P_bubble(comp,T,Z,params) 
    Pd = P_dew(comp,T,Z,params) 

    if P>=Pb:
        X = Z 
        Y = 0*Z 
        alpha = 0 

        return alpha,X,Y 
    
    elif P<=Pd:
        X = 0*Z 
        Y = Z 
        alpha = 1 

        return alpha,X,Y 
    
    else:
        # Rachford Rice algorithm 
        # Not ready for non-condensables as of now 

        alpha = (Pb-P)/(Pb-Pd) 
        gammas = md.gammas_moduni(T,Z,params) 
        Psat = Vp(comp,T) 

        while True:
            K = gammas*Psat/P 

            def f(alpha):
                R = (K-1)*Z/(1+(K-1)*alpha) 
                return np.sum(R) 
            
            alpha_new = scopt.newton(f,alpha,maxiter=1000) 

            err = np.abs(alpha_new-alpha)/alpha_new 

            if err<1e-06:
                break 

            else:
                alpha = alpha_new 
                X = Z/(1+(K-1)*alpha) 
                gammas = md.gammas_moduni(T,X,params) 

        Y = K*X 

        return alpha,X,Y 
    

# prop = md.Compound('Propane') 
# but = md.Compound('Butane') 
# pent = md.Compound('Pentane') 
# hexx = md.Compound('Hexane') 

# wat = md.Compound('Water') 
# eth = md.Compound('Ethanol') 
# prop.add_groups(1,2) 
# prop.add_groups(2,1) 

# but.add_groups(1,2)
# but.add_groups(2,2) 

# pent.add_groups(1,2) 
# pent.add_groups(2,3) 

# hexx.add_groups(1,2) 
# hexx.add_groups(2,4) 

# comp = [prop,but,pent,hexx] 
# params = md.generate_params(comp) 

# Z = [0.3,0.1,0.15,0.45]
# print(VLE_flash_PT(comp,2,50+273.15,Z,params))
# wat.add_groups(16,1) 
# eth.add_groups(1,1) 
# eth.add_groups(2,1) 
# eth.add_groups(14,1) 

# comp = [eth,wat] 
# params = md.generate_params(comp) 

# print(P_bubble(comp,80+273.15,[0.5,0.5],params,bool = True)) 

# print(P_dew(comp,80+273.15,[0.6576469, 0.3423531],params,bool = True))  

