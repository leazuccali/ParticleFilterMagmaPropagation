###
from __future__ import division

from trajectoire_pas_i import *
import math

eps = 1e-14  #


# Match input variable order as closely as possible
# def calc_mogi(x,y,xoff=0,yoff=0,d=3e3,dV=1e6,nu=0.25,output='cyl'):
# def forward(E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN):
def forward(x, y, depth, length, width, dip, opening, strike,  slip = 0, rake = 0, nu=0.25):
    '''
    Calculate surface displacements for Okada85 dislocation model
   Args:
        x: x-coordinate(s) of the observation point(s), east direction (m)
        y: y-coordinate(s) of the observation point(s), north direction (m)
        depth: depth of the fault reference point  (m)
        length: length of the fault plane along strike (m)
        width: width of the fault plane along dip  (default = 1000m)
        dip: dip angle of the fault plane, (rad)
        opening: tensile opening (dislocation) across the fault plane (m)
        strike: strike angle of the fault plane, in degrees (converted to radians internally)
        slip: slip magnitude along the fault plane (used with rake to compute strike-slip and dip-slip components)
        rake: rake angle of the slip, in degrees (converted to radians internally)
        nu: Poisson's ratio of the medium
    Returns:
        ue: east-component surface displacement at each observation point
        un: north-component surface displacement at each observation point
        uz: vertical surface displacement at each observation point

    '''

    e = x
    n = y
    #print(strike)
    strike = np.deg2rad(strike)

    #dip = np.deg2rad(dip)
    rake = np.deg2rad(rake)
    L = length
    W = width

    U1 = np.cos(rake) * slip
    U2 = np.sin(rake) * slip
    U3 = opening


    d = depth + np.sin(dip) * W / 2
    ec = e + np.cos(strike) * np.cos(dip) * W / 2
    nc = n - np.sin(strike) * np.cos(dip) * W / 2
    x = np.cos(strike) * nc + np.sin(strike) * ec + L / 2
    y = np.sin(strike) * nc - np.cos(strike) * ec + np.cos(dip) * W
    p = y * np.cos(dip) + d * np.sin(dip)
    q = y * np.sin(dip) - d * np.cos(dip)

    ux = - U1 / (2 * np.pi) * chinnery(ux_ss, x, p, L, W, q, dip, nu) - \
         U2 / (2 * np.pi) * chinnery(ux_ds, x, p, L, W, q, dip, nu) + \
         U3 / (2 * np.pi) * chinnery(ux_tf, x, p, L, W, q, dip, nu)
    uy =  -U1 / (2* np.pi) * chinnery(uy_ss, x, p, L, W, q, dip, nu) - \
          U2 / (2*np.pi) * chinnery(uy_ds, x, p, L, W, q, dip, nu) + \
          U3 / (2*np.pi) * chinnery(uy_tf, x, p, L, W, q, dip, nu)

    uz = - U1 / (2 * np.pi) * chinnery(uz_ss, x, p, L, W, q, dip, nu) - \
         U2 / (2 * np.pi) * chinnery(uz_ds, x, p, L, W, q, dip, nu) + \
         U3 / (2 * np.pi) * chinnery(uz_tf, x, p, L, W, q, dip, nu)

    ue = np.sin(strike) * ux - np.cos(strike) * uy
    un = np.cos(strike) * ux + np.sin(strike) * uy

    return ue, un, uz

#essai_froward = forward()

'''
% Notes for I... and K... subfunctions:
%
%	1. original formulas use Lame's parameters as mu/(mu+lambda) which
%	   depends only on the Poisson's ratio = 1 - 2*nu
%	2. tests for cos(dip) == 0 are made with "cos(dip) > eps"
%	   because cos(90*np.pi/180) is not zero but = 6.1232e-17 (!)
%	   NOTE: don't use cosd and sind because of incompatibility
%	   with Matlab v6 and earlier...
'''


def chinnery(f, x, p, L, W, q, dip, nu):
    ''' % Chinnery's notation [equation (24) p. 1143]'''
    u = (f(x, p, q, dip, nu) -
         f(x, p - W, q, dip, nu) -
         f(x - L, p, q, dip, nu) +
         f(x - L, p - W, q, dip, nu))
    return u


'''
% Displacement subfunctions
% strike-slip displacement subfunctions [equation (25) p. 1144]
'''


def ux_ss(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = xi * q / (R * (R + eta)) + \
        I1(xi, eta, q, dip, nu, R) * np.sin(dip)
    k = (q != 0)
    # u[k] = u[k] + np.arctan2( xi[k] * (eta[k]) , (q[k] * (R[k])))
    u[k] = u[k] + np.arctan((xi[k] * eta[k]) / (q[k] * R[k]))
    return u


def uy_ss(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = (eta * np.cos(dip) + q * np.sin(dip)) * q / (R * (R + eta)) + \
        q * np.cos(dip) / (R + eta) + \
        I2(eta, q, dip, nu, R) * np.sin(dip)
    return u


def uz_ss(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    db = eta * np.sin(dip) - q * np.cos(dip)
    u = (eta * np.sin(dip) - q * np.cos(dip)) * q / (R * (R + eta)) + \
        q * np.sin(dip) / (R + eta) + \
        I4(db, eta, q, dip, nu, R) * np.sin(dip)
    return u


def ux_ds(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = q / R - \
        I3(eta, q, dip, nu, R) * np.sin(dip) * np.cos(dip)
    return u


def uy_ds(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = ((eta * np.cos(dip) + q * np.sin(dip)) * q / (R * (R + xi)) -
         I1(xi, eta, q, dip, nu, R) * np.sin(dip) * np.cos(dip))
    k = (q != 0)
    u[k] = u[k] + np.cos(dip) * np.arctan((xi[k] * eta[k]) / (q[k] * R[k]))
    return u


def uz_ds(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    db = eta * np.sin(dip) - q * np.cos(dip)
    u = (db * q / (R * (R + xi)) -
         I5(xi, eta, q, dip, nu, R, db) * np.sin(dip) * np.cos(dip))
    k = (q != 0)
    # u[k] = u[k] + np.sin(dip) * np.arctan2(xi[k] * eta[k] , q[k] * R[k])
    u[k] = u[k] + np.sin(dip) * np.arctan((xi[k] * eta[k]) / (q[k] * R[k]))
    return u


def ux_tf(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = q ** 2 / (R * (R + eta)) - \
        I3(eta, q, dip, nu, R) * (np.sin(dip) ** 2)
    return u


def uy_tf(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = - (eta * np.sin(dip) - q * np.cos(dip)) * q / (R * (R + xi)) - \
        np.sin(dip) * xi * q / (R * (R + eta)) - \
        I1(xi, eta, q, dip, nu, R) * (np.sin(dip) ** 2)
    k = (q != 0)
    # u[k] = u[k] + np.sin(dip) * np.arctan2(xi[k] * eta[k] , q[k] * R[k])
    u[k] = u[k] + np.sin(dip) * np.arctan((xi[k] * eta[k]) / (q[k] * R[k]))
    return u


def uz_tf(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    db = eta * np.sin(dip) - q * np.cos(dip)
    u = (eta * np.cos(dip) + q * np.sin(dip)) * q / (R * (R + xi)) + \
        np.cos(dip) * xi * q / (R * (R + eta)) - \
        I5(xi, eta, q, dip, nu, R, db) * np.sin(dip) ** 2
    k = (q != 0)  # not at depth=0?
    u[k] = u[k] - np.cos(dip) * np.arctan((xi[k] * eta[k]) / (q[k] * R[k]))
    return u


def I1(xi, eta, q, dip, nu, R):
    db = eta * np.sin(dip) - q * np.cos(dip)
    if np.cos(dip) > eps:
        I = (1 - 2 * nu) * (- xi / (np.cos(dip) * (R + db))) - \
            np.sin(dip) / np.cos(dip) * \
            I5(xi, eta, q, dip, nu, R, db)
    else:
        I = -(1 - 2 * nu) / 2 * xi * q / (R + db) ** 2
    return I


def I2(eta, q, dip, nu, R):
    I = (1 - 2 * nu) * (-np.log(R + eta)) - \
        I3(eta, q, dip, nu, R)
    return I


def I3(eta, q, dip, nu, R):
    yb = eta * np.cos(dip) + q * np.sin(dip)
    #print(yb)
    db = eta * np.sin(dip) - q * np.cos(dip)
    if np.cos(dip) > eps:
        I = (1 - 2 * nu) * (yb / (np.cos(dip) * (R + db)) - np.log(R + eta)) + \
            np.sin(dip) / np.cos(dip) * \
            I4(db, eta, q, dip, nu, R)
    else:
        I = (1 - 2 * nu) / 2 * (eta / (R + db) + yb * q / (R + db) ** 2 - np.log(R + eta))
    return I


def I4(db, eta, q, dip, nu, R):
    if np.cos(dip) > eps:
        I = (1 - 2 * nu) * 1.0 / np.cos(dip) * \
            (np.log(R + db) - np.sin(dip) * np.log(R + eta))
    else:
        I = - (1 - 2 * nu) * q / (R + db)
    return I


def I5(xi, eta, q, dip, nu, R, db):
    X = np.sqrt(xi ** 2 + q ** 2)
    if np.cos(dip) > eps:
        I = (1 - 2 * nu) * 2 / np.cos(dip) * \
            np.arctan((eta * (X + q * np.cos(dip)) + X * (R + X) * np.sin(dip)) /
                      (xi * (R + X) * np.cos(dip)))
        I[xi == 0] = 0
    else:
        I = -(1 - 2 * nu) * xi * np.sin(dip) / (R + db)
    return I






def une_particule_deplacement_pas_i(X_haut, Z_haut, X_bas, Z_bas, longueur, ouverture, distance_parcourue_i, R, G, mu, E, vol, P_load, pas_okada, xmin, xmax, type_observations):
    '''

    Args:
        X_haut, Z_haut : coordonnées (X,Z) du front du dike
        X_bas, Z_bas : coordonnées (X,Z) du bas du dike
        longueur, ouverture : longueur et ouverture théorique du dike
        distance_parcourue_i : distance parcourue par le dike depuis le point de départ
        R, G, mu, E, vol : paramètres de la particule
        P_load : valeur charge ou décharge
        xmin, xmax : étendue axe des abscisses
        pas_okada : discrétisation de la grille
        type_observations :  'regulier','normal','random'

    Returns:
        ux : vecteur déplacement ux associé au dike
        uz: vecteur déplacement uz associé au dike
        grille_x : liste des points observés à la surface
        vec_x_c, vec_z_c : point central entre le haut et le bas du dike
        longueur_okada, ouverture_okada : longueur et ouverture pour okada
        vec_dip : valeur du dip
        strike : 180 si X_bas >=0; 0 sinon

    '''



    longueur_okada = np.sqrt( (X_haut - X_bas)**2 + (Z_haut - Z_bas)**2 )


    vol_eff = (((1 - 0.25**2)*G*abs(P_load)) * distance_parcourue_i**3 ) / E

    ouverture_okada = (math.pi * vol_eff) / (2 * longueur_okada**2)


    if type_observations == 'regulier':
    # ## SOIT
        grille_x = np.arange(xmin, xmax, pas_okada)
        grille_y = np.arange(0, pas_okada, pas_okada)
        mX, mY = np.meshgrid(grille_x, grille_y)
    else:
    ### SOIT
        grille_x = pas_okada
        liste_pas_okada = pas_okada.tolist()
        len_liste_pas_okada = len(liste_pas_okada)
        mX = np.array([[liste_pas_okada]])
        liste_zeros = [0] * len_liste_pas_okada
        mY = np.array([[liste_zeros]])


    #print("mX", mX)
    #print("mY", mY)

    ### Calculs du centre entre le haut et la bas de la trajectoire en chaque temps
    vec_x_c = (X_haut + X_bas) / 2
    # print('vecxc', vec_x_c)
    vec_z_c = -(Z_haut + Z_bas) / 2
    # print('veczc', vec_z_c)
    vec_z_c_abs = abs(vec_z_c)

    # nb_points = len(vec_X_haut)

    # Calcul du strike
    if X_bas >= 0:
        strike = 180
    else:
        strike = 0

    ##Calcul du dip pour chaque temps
    val_tan = (Z_haut - Z_bas) / (X_haut - X_bas)
    vec_dip = np.arctan(abs(val_tan))

    # Calcul des vecteurs de déplacement ue, un et uz
    ux = []
    uy = []
    uz = []

    ux, uy, uz = forward((mX - vec_x_c), mY, vec_z_c_abs, longueur_okada, 1000, vec_dip, ouverture_okada, strike)



    return ux, uz, grille_x, mY, vec_x_c, vec_z_c, longueur_okada, ouverture_okada, vec_dip, strike
