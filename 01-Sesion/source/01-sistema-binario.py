from vpython import *

scene.caption = ""
scene.forward = vector(0,-.3,-1)

G = 6.7e-11 # Constante Gravitacional

#Primer Objeto
giant = sphere(pos=vector(-1e11,0,0), radius=2e10, color=color.red, 
                make_trail=True, trail_type='points', interval=10, retain=50)
# Su masa
giant.mass = 2e30
giant.p = vector(0, 0, -1e4) * giant.mass

#Segundo Objeto
dwarf = sphere(pos=vector(1.5e11,0,0), radius=1e10, color=color.yellow,
                make_trail=True, interval=10, retain=50)
# Su masa
dwarf.mass = 1e30
dwarf.p = -giant.p

# Paso de tiempo
dt = 1e5
while True:
    rate(200)
    r = dwarf.pos - giant.pos # Distancia entre ambos objetos
    F = G * giant.mass * dwarf.mass * r.hat / mag2(r) # Fuerza gravitacional    
    # ¿Que unidades fisicas tiene giant.p?     
    giant.p = giant.p + F*dt
    # ¿Que unidades fisicas tiene dwarf.p?    
    dwarf.p = dwarf.p - F*dt
    # ¿A que ecuacion de cinematica se les parece estas ecuaciones?     
    giant.pos = giant.pos + (giant.p/giant.mass) * dt 
    dwarf.pos = dwarf.pos + (dwarf.p/dwarf.mass) * dt