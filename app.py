import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

st.set_page_config(page_title="Simulador cinemático", layout="wide")
st.title("Simulador cinemático: paquete desde un dron")
st.caption("Modelo del parcial. Eje y positivo hacia abajo; el impacto ocurre cuando y = H.")

def bisect(f, lo, hi, tol=1e-10):
    flo=f(lo)
    for _ in range(200):
        mid=(lo+hi)/2; fm=f(mid)
        if abs(fm)<tol or hi-lo<tol: return mid
        if flo*fm<=0: hi=mid
        else: lo,flo=mid,fm
    return (lo+hi)/2

def simulate(p, dt=0.005):
    H,vD,vw,vP,ta=(p[k] for k in ("H","vD","vw","vP","ta"))
    tx1,ty1,tx2,ty2=(p[k] for k in ("tx1","ty1","tx2","ty2"))
    VT1,VT2=p["VT1"],p["VT2"]
    y1=lambda t: VT1*(t-ty1*(1-np.exp(-t/ty1)))
    if y1(ta)>=H: tf=bisect(lambda t:y1(t)-H,0,ta)
    else:
        ya=y1(ta); vya=VT1*(1-np.exp(-ta/ty1))
        y2=lambda t: ya+VT2*(t-ta)+(vya-VT2)*ty2*(1-np.exp(-(t-ta)/ty2))
        hi=ta+max(1,H/max(VT2,.1))
        while y2(hi)<H: hi=ta+2*(hi-ta)
        tf=bisect(lambda t:y2(t)-H,ta,hi)
    t=np.append(np.arange(0,tf,dt),tf)
    x=np.zeros_like(t); y=np.zeros_like(t); vx=np.zeros_like(t); vy=np.zeros_like(t); ax=np.zeros_like(t); ay=np.zeros_like(t)
    b=t<ta; a=~b
    vx[b]=vw+(vD-vw)*np.exp(-t[b]/tx1); vy[b]=VT1*(1-np.exp(-t[b]/ty1))
    ax[b]=-(vD-vw)/tx1*np.exp(-t[b]/tx1); ay[b]=VT1/ty1*np.exp(-t[b]/ty1)
    vxa=vw+(vD-vw)*np.exp(-ta/tx1); vya=VT1*(1-np.exp(-ta/ty1))
    xa=vw*ta+(vD-vw)*tx1*(1-np.exp(-ta/tx1)); ya=y1(ta); s=t[a]-ta
    vx[a]=vw+(vxa-vw)*np.exp(-s/tx2); vy[a]=VT2+(vya-VT2)*np.exp(-s/ty2)
    ax[a]=-(vxa-vw)/tx2*np.exp(-s/tx2); ay[a]=-(vya-VT2)/ty2*np.exp(-s/ty2)
    x[b]=vw*t[b]+(vD-vw)*tx1*(1-np.exp(-t[b]/tx1)); y[b]=y1(t[b])
    x[a]=xa+vw*s+(vxa-vw)*tx2*(1-np.exp(-s/tx2))
    y[a]=ya+VT2*s+(vya-VT2)*ty2*(1-np.exp(-s/ty2))
    return dict(t=t,x=x,y=y,vx=vx,vy=vy,ax=ax,ay=ay,tf=tf)

with st.sidebar:
    st.header("Parámetros")
    H=st.slider("Altura del dron H (m)",20.0,150.0,80.0,1.0)
    vw=st.slider("Velocidad del viento vw (m/s)",-5.0,20.0,6.0,0.5)
    vP=st.slider("Velocidad de plataforma vP (m/s)",0.0,25.0,10.0,0.5)
    ta=st.slider("Instante de apertura ta (s)",1.0,4.5,3.0,0.1)
    vD=st.slider("Velocidad del dron vD (m/s)",0.0,35.0,18.0,0.5)
    VT1=st.slider("Velocidad terminal antes VT1 (m/s)",10.0,70.0,45.0,1.0)
    VT2=st.slider("Velocidad terminal después VT2 (m/s)",1.0,15.0,5.0,0.5)

p={"H":H,"vD":vD,"vw":vw,"vP":vP,"ta":ta,"tx1":5.0,"ty1":4.0,"tx2":2.0,"ty2":1.2,"VT1":VT1,"VT2":VT2}
s=simulate(p); t=s["t"]; tf=s["tf"]; xP0=s["x"][-1]-vP*tf; xp=xP0+vP*t
vx,vy=s["vx"][-1],s["vy"][-1]; relx=vx-vP; okv=vy<=6; okx=abs(relx)<=4
c1,c2,c3=st.columns(3)
c1.metric("Tiempo de impacto",f"{tf:.3f} s")
c2.metric("Posición horizontal al impacto",f"{s['x'][-1]:.3f} m")
c3.metric("Posición inicial requerida xP0",f"{xP0:.3f} m")
fig,axs=plt.subplots(2,3,figsize=(16,8))
axs[0,0].plot(s["x"],s["y"],label="Paquete"); axs[0,0].axhline(H,ls="--",label="y=H"); axs[0,0].scatter([s["x"][-1]],[H]); axs[0,0].invert_yaxis(); axs[0,0].set(title="Trayectoria",xlabel="x (m)",ylabel="y (m; positivo abajo)")
axs[0,1].plot(t,s["x"],label="Paquete"); axs[0,1].plot(t,xp,"--",label="Centro plataforma"); axs[0,1].set(title="Posición horizontal",xlabel="t (s)",ylabel="x (m)")
axs[0,2].plot(t,s["y"]); axs[0,2].axhline(H,ls="--"); axs[0,2].set(title="Posición vertical",xlabel="t (s)",ylabel="y (m)")
axs[1,0].plot(t,s["vx"],label="vx"); axs[1,0].plot(t,s["vy"],label="vy"); axs[1,0].set(title="Velocidades",xlabel="t (s)",ylabel="m/s")
axs[1,1].plot(t,s["ax"],label="ax"); axs[1,1].plot(t,s["ay"],label="ay"); axs[1,1].set(title="Aceleraciones",xlabel="t (s)",ylabel="m/s²")
axs[1,2].axis("off")
for ax in axs.flat[:5]:
    ax.axvline(ta,ls=":",color="gray",label="Apertura"); ax.axvline(tf,ls="-.",color="black",label="Impacto"); ax.grid(alpha=.3); ax.legend(fontsize=8)
axs[1,2].text(0,1,f"RESULTADOS EN EL IMPACTO\\n\\nvx = {vx:.3f} m/s\\nvy = {vy:.3f} m/s\\nax = {s['ax'][-1]:.3f} m/s²\\nay = {s['ay'][-1]:.3f} m/s²\\n\\nVelocidad relativa horizontal = {relx:.3f} m/s\\n\\nSeguridad vertical (vy ≤ 6): {'CUMPLE' if okv else 'NO CUMPLE'}\\nSeguridad horizontal (|vx-vP| ≤ 4): {'CUMPLE' if okx else 'NO CUMPLE'}\\nSeguridad total: {'CUMPLE' if okv and okx else 'NO CUMPLE'}",va="top",family="monospace",fontsize=10)
fig.suptitle("Simulador cinemático del parcial",fontsize=15); fig.tight_layout()
st.pyplot(fig)
st.info("El programa ajusta xP0 para que el centro de la plataforma coincida con el paquete en el instante de impacto. Las constantes temporales se mantienen en los valores del enunciado.")
