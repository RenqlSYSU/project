import numpy as np
import gc #garbage collector
ga = 9.80665 # Gravitational acceleration
a  = 6378388 # the radius of earth, m
Cpd= 1005.7 # J/(K*kg)
Rd = 287.05 # J/(K*kg)
epsl = 0.622 # Mv/Md

def calc_uv2vr_cfd(u,v,lat,lon):
    # u,v 4 deminsion; lat,lon 1 dimension
    # dv/dx-du/dy
    lat1 = lat*np.pi/180.0
    lon1 = lon*np.pi/180.0
    coslat = np.broadcast_to(np.cos(lat1).reshape(
        1,len(lat1),1), u.shape)
    term = 1.0/(a*coslat)*(center_diff(v,lon1,2)-
        center_diff(u*coslat,lat1,1))
    return term*100000

def calc_n2(t,z,level,varname,nmask):
    if varname == 'N2':
        term = calc_pot_temp(t, np.array(level), 1)
        term1 = upward_diff_z(term,z)*ga/term
    if varname == 'dthdz':
        lev4d = np.broadcast_to(np.array(level).reshape(
            len(level), 1, 1), t.shape)
        term1 = np.power(1000.0/lev4d,0.286)*upward_diff_z(
                t,z) + t*upward_diff_z(
                np.power(1000.0/lev4d,0.286),z)
        #term = calc_pot_temp(t, np.array(level), 1)
        #term1 = upward_diff_z(term,z)
    if varname == 'dthdz_moist':
        term1 = calc_dthdz_moist(t,z,np.array(level))
    
    nx = sum(x<0 for x in term1.flatten())
    print('negative grid is %d'%nx)
    if nmask:
        term1 = np.ma.array(term1,mask=(term1<0))
    return term1

def calc_eff_egr(t,z,u,f0,level,w):
    term = calc_eff_dthdz(t,z,level,w,False)
    t = calc_pot_temp(t, np.array(level), 1)
    brunt = np.sqrt(term*ga/t) 
    del w,t
    gc.collect()
    # negative is nan

    brunt = np.ma.array(brunt,mask=(brunt<1e-10))
    f00 = np.broadcast_to(f0.reshape(len(f0), 1), u.shape)
    term = 0.3098*np.abs(upward_diff_z(u,z))*f00/brunt
    return term 

def calc_eff_dthdz(t,z,level,w,nmask):
    term = calc_pot_temp(t, np.array(level), 1)
    dthdz = upward_diff_z(term,z)
    dthdzm = calc_dthdz_moist(t,z,np.array(level))
    term1 = dthdz + np.where(w>0,0,-1)*dthdzm
    #term1 = dtdz + np.where(w>0,0,w)*dthdzm
    
    nx = sum(x<0 for x in term1.flatten())
    print('negative grid is %d'%nx)
    if nmask:
        term1 = np.ma.array(term1,mask=(term1<0))
    return term1

def calc_dthdz_moist(t,z,level):
    lev4d = np.broadcast_to(level.reshape(len(level), 1, 1), t.shape)
    Lv = 2.5e6-2323.0*(t-273.15) # Latent Heat of Vaporization of Water
    es = 6.1078*np.exp((17.13*t-273.15)/(t-38.0)) # saturation vapour pressure
    rs = epsl*es/(lev4d*100.0-es) # Saturation mixing ratio mv/md
    dtdz = -(ga/Cpd)*(1+Lv*rs/(Rd*t))/(1+epsl*Lv*Lv*rs/(Cpd*Rd*t*t)) # dtdz
    dthdz = np.power(1000.0/lev4d,0.286)*dtdz + t*upward_diff_z(
            np.power(1000.0/lev4d,0.286),z)
    return dthdz

def calc_egr(t,z,u,f0,level):
    t = calc_pot_temp(t, np.array(level), 1)
    brunt = np.sqrt(upward_diff_z(t,z)*ga/t) 
    # negative is nan

    brunt = np.ma.array(brunt,mask=(brunt<1e-10))
    f00 = np.broadcast_to(f0.reshape(len(f0), 1), u.shape)
    term = 0.3098*np.abs(upward_diff_z(u,z))*f00/brunt
    #dudz = np.sqrt(np.power(upward_diff_z(u, z),2)+
    #    np.power(upward_diff_z(v, z),2))
    #term = 0.3098*dudz*f00/brunt
    return term

def calc_teqv(t,q):
    Lv   = 2.5104e6  # [J/kg]=[m2/s2]  Latent Heat of Vaporization of Water
    return (t+(Lv/Cpd)*q)  # equivalent temperature 

def calc_pot_temp(t,p,dim):
    # dim indicating which dimension of t 
    # is similar to p, unit of p is hPa
    term = np.moveaxis(t,dim,-1)*np.power(1000.0/p,0.286)
    return np.moveaxis(term,-1,dim)

def upward_diff_z(var,z):
    term1 = var.copy()
    term1[:,0:-1,:,:] = np.divide(var[:,0:-1,:,:]-var[:,1:,:,:],
        z[:,0:-1,:,:]-z[:,1:,:,:])
    term1[:,-1,:,:] = term1[:,-2,:,:] 
    return term1

def center_diff_z(var,z):
    term1 = var.copy()
    term1[:,1:-1,:,:] = np.divide(var[:,0:-2,:,:]-var[:,2:,:,:],
        z[:,0:-2,:,:]-z[:,2:,:,:])
    term1[:,0,:,:] = np.divide(var[:,0,:,:]-var[:,1,:,:],
        z[:,0,:,:]-z[:,1,:,:])
    term1[:,-1,:,:] = np.divide(var[:,-2,:,:]-var[:,-1,:,:],
        z[:,-2,:,:]-z[:,-1,:,:])
    return term1

def center_diff(var,z,dim):
    varn = np.moveaxis(var,dim,-1)
    term = varn.copy()
    if len(term.shape)==4:
        term[:,:,:,1:-1] = np.divide(varn[:,:,:,0:-2]-varn[:,:,:,2:],
            z[0:-2]-z[2:])
        term[:,:,:,0] = np.divide(varn[:,:,:,0]-varn[:,:,:,1],
            z[0]-z[1])
        term[:,:,:,-1] = np.divide(varn[:,:,:,-2]-varn[:,:,:,-1],
            z[-2]-z[-1])
    if len(term.shape)==3:
        term[:,:,1:-1] = np.divide(varn[:,:,0:-2]-varn[:,:,2:],
            z[0:-2]-z[2:])
        term[:,:,0] = np.divide(varn[:,:,0]-varn[:,:,1],
            z[0]-z[1])
        term[:,:,-1] = np.divide(varn[:,:,-2]-varn[:,:,-1],
            z[-2]-z[-1])
    return np.moveaxis(term,-1,dim)

