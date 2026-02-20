import xlwings as xw 
import numpy as np 

wbk = xw.Book('MODUNI.xlsx') 
data = xw.Book('Data.xlsx')

gc = wbk.sheets('Groups') 
wa = wbk.sheets('A')  
wb = wbk.sheets('B')  
wc = wbk.sheets('C')  
vp = data.sheets('Data') 

comp_names = vp.range('C3:C347').value 
R = 8314

def Psat(T,C): 
    C1,C2,C3,C4,C5 = C 
    ans = C1+C2/T+C3*np.log(T)+C4*T**C5 
    return np.exp(ans) 

class Compound: 

    def __init__(self,name):
        self.name = name 
        self.groups = []  
        found = 0 

        for elem in comp_names: 
            if elem == self.name: 
                ind = comp_names.index(elem) 
                self.vpconsts = vp.range('D'+str(ind+3)+':H'+str(ind+3)).value 
                more_data = vp.range('I'+str(ind+3)+':N'+str(ind+3)).value  
                mw,Tc,Pc,Vc,Zc,w = more_data 
                self.mw = mw 
                self.Tc = Tc 
                self.Pc = Pc*10 
                self.w = w
                found = 1 
                break 
        
        if found ==0: 
            self.vpconsts = np.zeros(5) 
            self.mw = 0 
            self.Tc = 0 
            self.Pc = 0 
            self.w = 0 


    def add_groups(self,sub_grp_no,quant): 
        row = sub_grp_no + 2 
        main_grp_no = int(gc.range(row,4).value)
        self.groups.append([sub_grp_no,quant,main_grp_no])  

    def Vp(self,T): 
        consts = self.vpconsts 
        val = Psat(T,consts) 
        return val

def generate_params(compounds): 

    sub_groups = [] 
    main_groups =[] 
    Q = [] 
    R = []

    for comp in compounds: 
        info = comp.groups 

        for grp in info: 
            [sg,quant,mg] = grp 
            match = 0 

            for element in sub_groups: 
                
                if sg==element: 
                    match = 1 
                    break # Confirm later
            
            if match ==0: 
                sub_groups.append(sg) 
                main_groups.append(mg) 
                val_q = gc.range(2+int(sg),7).value 
                val_r = gc.range(2+int(sg),6).value 
                Q.append(val_q) 
                R.append(val_r) 

    a=[] 
    b=[] 
    c=[]

    for i in main_groups: 
        ai =[] 
        bi =[] 
        ci =[]

        row = 2+int(i) 

        for j in main_groups: 
            col = 2+int(j) 
            vala = wa.range(row,col).value 
            valb = wb.range(row,col).value 
            valc = wc.range(row,col).value 
            ai.append(vala)
            bi.append(valb)
            ci.append(valc)

        a.append(ai) 
        b.append(bi) 
        c.append(ci) 

    nu =[] 

    for comp in compounds: 
        info = comp.groups 
        nu_comp = [] 

        for element in sub_groups: 
            match = 0 

            for grp in info: 
                [sg,quant,mg] = grp 

                if sg==element: 
                    nu_comp.append(quant)
                    match = 1  
                    break 

            if match==0: 
                nu_comp.append(0) 

        nu.append(nu_comp) 

    a = np.array(a) 
    b = np.array(b) 
    c = np.array(c) 
    R = np.array(R) 
    Q = np.array(Q) 
    nu = np.array(nu) 
    return R,Q,nu,a,b,c

def combinatorial(x,nu,R,Q): 
    x= np.array(x) 
    q = nu.dot(Q) 
    r = nu.dot(R) 

    V = r/np.sum(x*r) 
    F = q/np.sum(x*q) 
    V_p = (r**(3/4))/np.sum((r**(3/4))*x)

    ans = 1-V_p+np.log(V_p)-5*q*(1-V/F+np.log(V/F)) 
    return ans 

def gamma_mixture(T,x,nu,Q,a,b,c): 
    x = np.array(x) 
    mat = x.dot(nu) 
    X = mat/np.sum(mat) 
    theta = (Q*X)/np.sum(Q*X) 
    psi = np.exp(-(a/T + b + c*T)) 

    t1 = theta.dot(psi) 
    t2 = psi.dot(theta/t1) 

    ans = Q*(1-np.log(t1)-t2) 
    return ans 

def gamma_pure(T,x,nu,Q,a,b,c): 
    gk = [] 
    n = len(x)
    i=0
    for element in nu: 
       x_c = np.zeros(n) 
       x_c[i] = 1
       gk.append(gamma_mixture(T,x_c,nu,Q,a,b,c)) 
       i+=1        
    gk = np.array(gk)

    return gk 

def residual(T,x,nu,Q,a,b,c): 
    gamma_p = gamma_pure(T,x,nu,Q,a,b,c) 
    gamma_mix = gamma_mixture(T,x,nu,Q,a,b,c) 
    gammas_r = [] 
    i=0
    for element in nu: 
        gam_p = gamma_p[i] 
        diff = gamma_mix-gam_p
        gammas_r.append(np.sum(element*diff)) 
        i+=1

    gammas_r = np.array(gammas_r)
    return gammas_r


def gammas_moduni(T,x,params): 
    x = np.array(x) 
    R,Q,nu,a,b,c = params 

    ln_gc = combinatorial(x,nu,R,Q) 
    ln_gr = residual(T,x,nu,Q,a,b,c) 

    ln_g = ln_gc + ln_gr 
    ans = np.exp(ln_g)

    return ans

def list_subgrps(): 
    list_comp = gc.range('B2:C127').value 
    return list_comp 

# Uses Modified Raoult's Law
def find_Tbubble_binary(P,x,comp,params,Tg=350,alpha = 20): 
    X = np.array([x,1-x]) 
    comp1,comp2 = comp 
    
    while True: 
        gammas = gammas_moduni(Tg,X,params) 
        Psat= np.array([comp1.Vp(Tg),comp2.Vp(Tg)])*1e-05
        rhs = X*gammas*Psat 
        rhs = np.sum(rhs) 

        if np.abs(rhs-P)/P < 1e-04: 
            break 

        else: 
            delT = alpha*(P-rhs)  
            Tg+=delT   

    return Tg  

# Uses Modified Raoult's Law
def generate_binary_VLE_isobaric(P,comp,params,points = 100,T_g = 350,alpha_val =20): 
    T = []
    X= np.linspace(0,1,points + 1) 
    gamma1=[] 
    
    for element in X: 
        val = find_Tbubble_binary(P,element,comp,params,Tg=T_g,alpha=alpha_val) 
        T.append(val)
        gammas = gammas_moduni(val,[element,1-element],params)
        gam1,_= gammas
        gamma1.append(gam1) 
        
    T = np.array(T)
    gamma1 = np.array(gamma1) 
    comp1 = comp[0]
    Y = (X*gamma1*comp1.Vp(T)*1e-05)/P 

    return T,X,Y 

# Uses Modified Raoult's Law
def generate_binary_VLE_isothermal(T,comp,params,points = 100): 
    X = np.linspace(0,1,points +1) 
    P = [] 
    gamma1=[]
    comp1,comp2 = comp

    for element in X: 
        x = np.array([element,1-element])
        gammas = gammas_moduni(T,x,params) 
        gam1,_ = gammas
        Psat= np.array([comp1.Vp(T),comp2.Vp(T)])*1e-05 
        val = x*gammas*Psat 
        P.append(np.sum(val)) 
        gamma1.append(gam1)

    P = np.array(P) 
    gamma1 = np.array(gamma1) 
    Y = (X*gamma1*comp1.Vp(T)*1e-05)/P 

    return P,X,Y

# Programming PR-EOS for High Pressure VLE (Gamma - Phi Approach) 

def eos_constants(T,P,comp): 
    P = P*1e+05  
    Tc = comp.Tc 
    w = comp.w 
    Pc = (comp.Pc)*1e+05 
    Tr = T/Tc 
    m = 0.37464 + 1.5422*w - 0.26992*w**2 
    alpha = (1+m*(1-np.sqrt(Tr)))**2 
    a = 0.45724*alpha*(R*Tc)**2/Pc 
    b = 0.07779*R*Tc/Pc 
    da_dT = (-a*m)/np.sqrt(alpha*T*Tc) 

    return a,b,da_dT

def Z_vap(T,P,a,b): 
    P = P*1e+05 
    A = a*P/(R*T)**2 
    B = b*P/(R*T) 

    p0 = (B**3 + B**2 -A*B) 
    p1 = (A-2*B-3*B**2) 
    p2 = (B-1) 
    p3 = 1 

    roots_Z = np.roots([p3,p2,p1,p0]) 
    Z = roots_Z[roots_Z.imag==0] 
    z_vap = np.max(Z) 

    return z_vap 

def fug_vap(T,P,Z,a,b,da_dT): 
    P = P*1e+05 
    A = a*P/(R*T)**2 
    B = b*P/(R*T) 
    sq = np.sqrt(2) 

    hr_rt = (Z-1) + (1/R)*(da_dT - a/T)*(1/(2*b*sq))*np.log((Z + B*(1+sq))/(Z +B*(1-sq))) 
    sr_r = np.log(Z-B) + (da_dT/(2*R*b*sq))*np.log((Z + B*(1+sq))/(Z +B*(1-sq))) 

    gr_rt = hr_rt - sr_r 
    f_P = np.exp(gr_rt) 
    fug = f_P*P 
    fug = fug*1e-05 

    return fug 

def fug_vap_comps(T,P,compounds): 
    fug = [] 

    for comp in compounds: 
        a,b,da_dT = eos_constants(T,P,comp) 
        Z = Z_vap(T,P,a,b) 
        f = fug_vap(T,P,Z,a,b,da_dT) 
        fug.append(f) 

    fug = np.array(fug) 

    return fug 

# wat = Compound('Water') 
# eth = Compound('Ethanol') 

# wat.add_groups(16,1) 
# eth.add_groups(1,1) 
# eth.add_groups(2,1) 
# eth.add_groups(14,1) 

# comp = [eth,wat] 
# params = generate_params(comp) 

# print(fug_vap_comps(352,1.013,comp))  