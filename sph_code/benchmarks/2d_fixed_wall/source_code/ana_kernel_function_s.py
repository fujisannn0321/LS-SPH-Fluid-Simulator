#=========================#
#  module                 # 
#=========================#
import math

#=========================#
#  cubic spline           # 
#=========================#
def cubic_spline_d0W(r,h):
    q = r/h
    if (0.0 <= q <= 1.0):
        W = (2.0-q)**3.0 - 4.0*(1.0-q)**3.0
    elif(1.0 < q <= 2.0):
        W = (2.0-q)**3.0
    else:
        W = 0.0
    W = 5.0 / (14.0 * math.pi * h**2.0) * W
    return(W)

def cubic_spline_d1W(r,h):
    q = r/h
    if (0.0 <= q <= 1.0):
        dW = (2.0-q)**2.0 - 4.0*(1.0-q)**2.0
    elif(1.0 < q <= 2.0):
        dW = (2.0-q)**2.0
    else:
        dW = 0.0
    dW = 5.0 / (14.0 * math.pi * h**2.0) * dW * (-3.0/h)
    return(dW)

def cubic_spline_d2W(r,h):
    q = r/h
    if (0.0 <= q <= 1.0):
        ddW = (2.0-q) - 4.0*(1.0-q)
    elif(1.0 < q <= 2.0):
        ddW = (2.0-q)
    else:
        ddW = 0.0
    ddW = 5.0 / (14.0 * math.pi * h**2.0) * ddW * (6.0/h**2.0)
    return(ddW)

#=========================#
#  quintic spline         # 
#=========================#
def quintic_spline_d0W(r,h):
    q = r/h
    if (0.0 <= q <= 1.0):
        W = (3.0-q)**5.0 - 6.0*(2.0-q)**5.0 + 15.0*(1.0-q)**5.0
    elif(1.0 < q <= 2.0):
        W = (3.0-q)**5.0 - 6.0*(2.0-q)**5.0
    elif(2.0 < q <= 3.0):
        W = (3.0-q)**5.0
    else:
        W = 0.0
    W = 7.0 / (478.0 * math.pi * h**2.0) * W
    return(W)

def quintic_spline_d1W(r,h):
    q = r/h
    if (0.0 <= q <= 1.0):
        dW = (3.0-q)**4.0 - 6.0*(2.0-q)**4.0 + 15.0*(1.0-q)**4.0
    elif(1.0 < q <= 2.0):
        dW = (3.0-q)**4.0 - 6.0*(2.0-q)**4.0
    elif(2.0 < q <= 3.0):
        dW = (3.0-q)**4.0
    else:
        dW = 0.0
    dW = 7.0 / (478.0 * math.pi * h**2.0) * dW *(-5.0/h)
    return(dW)

def quintic_spline_d2W(r,h):
    q = r/h
    if (0.0 <= q <= 1.0):
        ddW = (3.0-q)**3.0 - 6.0*(2.0-q)**3.0 + 15.0*(1.0-q)**3.0
    elif(1.0 < q <= 2.0):
        ddW = (3.0-q)**3.0 - 6.0*(2.0-q)**3.0
    elif(2.0 < q <= 3.0):
        ddW = (3.0-q)**3.0
    else:
        ddW = 0.0
    ddW = 7.0 / (478.0 * math.pi * h**2.0) * ddW *(20.0/h**2.0)
    return(ddW)

#=========================#
#  Wendland C2            # 
#=========================#
def wendland_C2_d0W(r,h):
    q = r/h
    if (0.0 <= q <= 2.0):
        W = (1.0 - 0.5*q)**4.0 * (1.0 + 2.0*q)
    else:
        W = 0.0
    W = 7.0 / (4.0 * math.pi * h**2.0) * W
    return(W)

def wendland_C2_d1W(r,h):
    q = r/h
    if (0.0 <= q <= 2.0):
        dW = -2.0*(1.0 - 0.5*q)**3.0 * (1.0 + 2.0*q) \
            +2.0*(1.0 - 0.5*q)**4.0
    else:
        dW = 0.0
    dW = 7.0 / (4.0 * math.pi * h**2.0) * dW / h
    return(dW)

def wendland_C2_d2W(r,h):
    q = r/h
    if (0.0 <= q <= 2.0):
        ddW = +3.0*(1.0 - 0.5*q)**2.0 * (1.0 + 2.0*q) \
              -4.0*(1.0 - 0.5*q)**3.0 \
              -4.0*(1.0 - 0.5*q)**3.0
    else:
        ddW = 0.0
    ddW = 7.0 / (4.0 * math.pi * h**2.0) * ddW / h**2.0
    return(ddW)

#=========================#
#  Wendland C4            # 
#=========================#
def wendland_C4_d0W(r,h):
    q = r/h
    if (0.0 <= q <= 2.0):
        W = (1.0 - 0.5*q)**6.0 * (1.0 + 3.0*q + 35.0*q**2.0/12.0)
    else:
        W = 0.0
    W = 9.0 / (4.0 * math.pi * h**2.0) * W
    return(W)

def wendland_C4_d1W(r,h):
    q = r/h
    if (0.0 <= q <= 2.0):
        dW = -3.0*(1.0 - 0.5*q)**5.0 * (1.0 + 3.0*q + 35.0*q**2.0/12.0) \
             +1.0*(1.0 - 0.5*q)**6.0 * (3.0 + 35.0*q/6.0)
    else:
        dW = 0.0
    dW = 9.0 / (4.0 * math.pi * h**2.0) * dW / h
    return(dW)

def wendland_C4_d2W(r,h):
    q = r/h
    if (0.0 <= q <= 2.0):
        ddW = +15.0/2.0*(1.0 - 0.5*q)**4.0 * (1.0 + 3.0*q + 35.0*q**2.0/12.0) \
              -3.0     *(1.0 - 0.5*q)**5.0 * (3.0 + 35.0*q/6.0) \
              -3.0     *(1.0 - 0.5*q)**5.0 * (3.0 + 35.0*q/6.0) \
              +1.0     *(1.0 - 0.5*q)**6.0 * (35.0/6.0)
    else:
        ddW = 0.0
    ddW = 9.0 / (4.0 * math.pi * h**2.0) * ddW / h**2.0
    return(ddW)

#=========================#
#  Wendland C6            # 
#=========================#
def wendland_C6_d0W(r,h):
    q = r/h
    if (0.0 <= q <= 2.0):
        W = (1.0 - 0.5*q)**8.0 * (1.0 + 4.0*q + 25.0*q**2.0/4.0 + 4.0*q**3.0)
    else:
        W = 0.0
    W = 39.0 / (14.0 * math.pi * h**2.0) * W
    return(W)

def wendland_C6_d1W(r,h):
    q = r/h
    if (0.0 <= q <= 2.0):
        dW = -4.0*(1.0 - 0.5*q)**7.0 * (1.0 + 4.0*q + 25.0*q**2.0/4.0 + 4.0*q**3.0) \
             +1.0*(1.0 - 0.5*q)**8.0 * (4.0 + 25.0*q/2.0 + 12.0*q**2.0)
    else:
        dW = 0.0
    dW = 39.0 / (14.0 * math.pi * h**2.0) * dW / h
    return(dW)

def wendland_C6_d2W(r,h):
    q = r/h
    if (0.0 <= q <= 2.0):
        ddW = +14.0*(1.0 - 0.5*q)**6.0 * (1.0 + 4.0*q + 25.0*q**2.0/4.0 + 4.0*q**3.0) \
              -4.0 *(1.0 - 0.5*q)**7.0 * (4.0 + 25.0*q/2.0 + 12.0*q**2.0) \
              -4.0 *(1.0 - 0.5*q)**7.0 * (4.0 + 25.0*q/2.0 + 12.0*q**2.0) \
              +1.0 *(1.0 - 0.5*q)**8.0 * (25.0/2.0 + 24.0*q)
    else:
        ddW = 0.0
    ddW = 39.0 / (14.0 * math.pi * h**2.0) * ddW / h**2.0
    return(ddW)

# END #