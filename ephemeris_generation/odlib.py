#David Peng

from math import *
import numpy as np
import copy


#Returns scalar dot product of two vectors of any size 
def dot(vector1, vector2):
    dotProduct = 0
    for i in range(0, len(vector1)):
        dotProduct = dotProduct + vector1[i] * vector2[i]
        #adds product of component of each vector to the sum
    return dotProduct


#Returns vector cross product of two 3x3 vectors
def cross(v1, v2):
    mat1 = v1[1] * v2[2] - v1[2] * v2[1] #finds cross product of x component
    mat2 = v1[0] * v2[2] - v1[2] * v2[0] #finds cross product of y component
    mat3 = v1[0] * v2[1] - v1[1] * v2[0] #finds cross product of z component
    return [mat1, -mat2, mat3]


#Returns triple product of three vectors
def tri_product(v1, v2, v3):
    return dot(v1, cross(v2, v3))
    #finds dot product of v1 and the cross product of v2 and v3


#finds mean
def mean(x_list):
     x_sum = 0
     for i in x_list:
          x_sum += i
     x_avg = x_sum / len(x_list)
     return x_avg


#finds standard deviation
def stdev(x_list):
     x_sum = 0
     x_avg = mean(x_list)
     for i in x_list:
          x_sum += (i - x_avg)**2
     sigma = sqrt(x_sum / len(x_list))
     return sigma


#decimalizes RA
#takes RA in hours, minutes, seconds, and
#returns angle in decimalized degrees
def HMStoDeg (h, m, s):
    total_h = h + m/60 + s/3600
    return total_h * 360 / 24


#displays RA in hours, minutes, seconds to the reader,
#given decimalized value in degrees
def RAdecimalToHMS(dec):
     theHour = dec * 24 / 360
     intHour = int(theHour)
     theMinute = (theHour - intHour) * 60
     intMinute = int(theMinute)
     theSecond = (theMinute - intMinute) * 60
     return [intHour, abs(intMinute), abs(theSecond)]


#decimilaizes declination given degrees, minutes, seconds of arc
def DMStoDeg (deg, arcmin, arcsec):
    if (deg < 0): 
        mag_deg = -1 * deg
    else:
        mag_deg = deg
    dec = mag_deg + arcmin / 60. + arcsec / 3600.
    if (deg > 0):
        return dec
    else: # deg < 0
        return -1 * dec



#displays declination in degrees, minutes, seconds of arc,
#given decimalized value in degrees
def DECdecimalToDMS(dec):
     intDegree = int(dec)
     theMinute = (dec - intDegree) * 60
     intMinute = int(theMinute)
     theSecond = (theMinute - intMinute) * 60
     return[intDegree, abs(intMinute), abs(theSecond)]


#finds magnitude 
def mag(vector):
     s = 0
     for i in vector:
          s += i**2
     return sqrt(s)


#Find angles given sin and cos values
#quadrant
def find_angle(sin,cos):
##     return atan2(sin,cos)
     theta = 0
     if (sin > 0 and cos > 0): #quad 1
          theta = asin(sin)
     elif (sin > 0 and cos < 0): #quad 2
          theta = pi - asin(sin)
     elif (sin < 0 and cos < 0): #quad 3
          theta = pi + asin(-sin)
     elif sin < 0 and cos > 0: #quad 4
          theta = asin(sin)
     elif sin == 1 and cos == 0: #pi/2
          theta = pi/2
     elif sin == -1 and cos == 0: # 3/2 pi
          theta = 3/2 * pi
     elif sin == 0 and cos == 1:
          theta = 0
     else: # sin == 0 and cos == -1
          theta = pi
     return theta


#Given decimal radian values of two sides and the opposite angle
#finds the other side and angles 
def spherical_triangle(aa, bb, bbig_C, isDec):
     if (isDec):
          a = aa * pi / 180
          b = bb * pi / 180
          big_C = bbig_C * pi / 180
     else:
          a = aa
          b = bb
          big_C = bbig_C

     cos_c = cos(a)*cos(b) + sin(a)*sin(b)*cos(big_C)
     c = acos(cos_c)
     cos_bigB = (cos(b) - cos(a)*cos(c))/(sin(a)*sin(c))
     cos_bigA = (cos(a) - cos(b)*cos(c))/(sin(b)*sin(c))
     B = acos(cos_bigB)
     A = acos(cos_bigA)

     if (isDec):
          c *= 180/pi
          A *= 180/pi
          B *= 180/pi

     return [c, A, B]
     
     

#Finds distance between two longitudinal/latitudinal coordinates
def calc_distance(long1,lat1,long2,lat2, isEarth):
     #changing R drastically changes the final result...
     
     b = (90-lat1)*pi/180
     c = (90-lat2)*pi/180

     A = (long1 - long2) * pi/180
     if (A < 0):
          A = -A
     cos_g = cos(b)*cos(c) + sin(b)*sin(c)*cos(A)
     g = acos(cos_g)

     if (isEarth):
          R = 6000
          a = R * g
     else:
          a = g

     return a


#Multiplies matrices
#
def mat33x31(m1,m2):
     m = []
     for r in range(len(m1)):
          cur_r = 0
          for c in range(len(m1[r])):
               cur_r += m1[r][c] * m2[c]
          m.append(cur_r)
     return m


#Performs series of coordinate rotations (3)
def coordinate_rot(angles,vector):
     a = [i*pi/180 for i in angles]

     #about y     
     mat1 = [ [cos(a[0]), 0, -sin(a[0])], [0,1,0], [sin(a[0]), 0, cos(a[0])]]
     vi = mat33x31(mat1,vector)

     #about x
     mat2 = [ [1, 0, 0], [0, cos(a[1]), -sin(a[1])], [0, sin(a[1]), cos(a[1])]]
     vii = mat33x31(mat2, vi)

     #about z
     mat3 = [ [cos(a[2]), -sin(a[2]), 0], [sin(a[2]), cos(a[2]), 0], [0, 0, 1]]
     viii = mat33x31(mat3, vii)

     return viii


#performs single coordinate rotation
#angle: in degrees
#vector 3x1
#axis: 'x', y or z

def coor_rot(angle, vector, axis):
    a = angle * pi / 180 # convert to radians
    val = 0
    mat = np.zeros(3)
    if axis == 'y':
        mat = [  [cos(a), 0, -sin(a)],
                 [0,1,0],
                 [sin(a), 0, cos(a)]]
    elif axis == 'x':
        mat = [ [1, 0, 0],
                [0, cos(a), sin(a)],
                [0, -sin(a), cos(a)]]
    elif axis == 'z':
        mat = [ [cos(a), sin(a), 0],
                [-sin(a), cos(a), 0],
                [0, 0, 1]]
    val = mat33x31(mat, vector)
    return val


#Finds percent error in % (ex perfect accuracy would return 100)
def percent_error(correct, wrong):
     dif = correct - wrong
     percent = abs(dif/correct * 100)  
     return percent


#finds determinant of 2x2
def det2x2(m):
    det_m = m[0][0]*m[1][1]-m[0][1]*m[1][0]
    return det_m


#finds determinant of 3x3
def det3x3(m):
    a = m[0][0] * det2x2([[m[1][1],m[1][2]],[m[2][1],m[2][2]]])
    b = -m[0][1]*det2x2([[m[1][0],m[1][2]],[m[2][0],m[2][2]]])
    c = m[0][2]*det2x2([[m[1][0],m[1][1]],[m[2][0],m[2][1]]])
    det_m = a + b + c
    return det_m


#splices column. makes copy of matrix to be safe
def spliceColumn(o_m3x3, m3x1, j):
    m3x3 = copy.deepcopy(o_m3x3)
    m3x3[0][j] = m3x1[0]
    m3x3[1][j] = m3x1[1]
    m3x3[2][j] = m3x1[2]
    return m3x3


#performs cramer given a 3x3 matrix and a 3x1 matrix of constants
def cramer(m3x3, constant):
    D = det3x3(m3x3)
    Ax = spliceColumn(m3x3, constant, 0)
    By = spliceColumn(m3x3, constant, 1)
    Cz = spliceColumn(m3x3, constant, 2)
    x = det3x3(Ax) / D
    y = det3x3(By) / D
    z = det3x3(Cz) / D
    return [x, y, z]

    
#Finds r vector
def find_r(x, y, z):
     r_vec = [x, y, z]
     return r_vec


#Finds r dot vector
def find_r_dot(pos1, pos2):
     k = 0.0172029847

     
     dx = pos2[0] - pos1[0]
     dy = pos2[1] - pos1[1]
     dz = pos2[2] - pos1[2]
     
     delT = k * (mag(pos2) - mag(pos1))

     sx = dx/delT
     sy = dy/delT
     sz = dz/delT

     r_dot_vec = [sx, sy, sz]
     return r_dot_vec


#Finds h vector and magnitude
def find_h(r_vec, r_dot_vec):
     h_vec = cross(r_vec, r_dot_vec)
     return h_vec


#Finds v^2
def v_squared(r_dot):
     vv = dot(r_dot, r_dot)
     return vv


#OD Element a
#Finds a
#r is list with vector and last element magnitude
def find_a(r, vv):
##     print('2/r',2/r,'v^2', vv)
     a_vec = 1 / (2 / r - vv)
##     print('a',a_vec)
     return a_vec


#Finds e
def find_e(h, a):
     ee = sqrt(1 - h**2/a)
     return ee


#Finds i in rad
def find_i(h_v):
     i = acos(h_v[2] / mag(h_v))
     return i


#Find omega in rad
def find_omega(h_v, i):
     h = mag(h_v)
     cos_o = -h_v[1] / (h * sin(i))
     sin_o = h_v[0] / (h * sin(i))
     omega = find_angle(sin_o, cos_o)
     while omega < 0:
         omega += 2 * pi
##     print('omega',omega)
     return omega


#Finds U, to find W
def find_U(omega, r, i, pos):
     x = pos[0]
     y = pos[1]
     
     cos_u = (pos[0] * cos(omega) + pos[1] * sin(omega)) / r
     sin_u = pos[2] / (r * sin(i))
     U = find_angle(sin_u, cos_u)
     return U


#Finds V, to find W
def find_V(a, ee, h_v, r_v, r_dot_v):
     r = mag(r_v)
     h = mag(h_v)
     cos_v = 1 / ee * (a * (1 - ee**2) / r - 1)
     sin_v = 1 / ee * (a * (1 - ee**2) / h * dot(r_v, r_dot_v) / r)
     V = find_angle(sin_v, cos_v)
     return V

     
#Finds W, in rad
def find_W(U, V):
     if U < V:
          U += 2 * pi
     W = U - V
##     print('w',W)
     return W


#Finds E
def find_E(ee, r, a, V):
     E = acos((1 / ee) * (1 - r/a))
     if E > 2 * pi:
         E = E - 2 * pi
     return E


#Finds M
def find_M(E, ee):
     M = E - ee * sin(E)
     return M


#Finds rhohat using RA and DEC
def find_rhohat1(ra_deg, dec_deg):
     ra = radians(ra_deg)
     dec = radians(dec_deg)
     rho_x = cos(ra) * cos(dec)
     rho_y = sin(ra) * cos(dec)
     rho_z = sin(dec)
     return (rho_x, rho_y, rho_z)


#Finds rhohat using vectors from horizons
#rhohat / rhohat_v
def find_rhohat2(pos):
     rho = mag(pos)
     rho_v = np.array(pos)
     rhohat = rho_v / rho
     
     return rhohat


# Finds sqrt sum of components squared
# l is a list
def sqrt_sum(l):
     s = 0
     for i in len(l):
          s += i**2
     s = sqrt(s)
     return s


#t2 in JD
#m2 in radians
def find_T(t2, m2, a):
    k = 0.0172029847
    n = sqrt(1/a**3)
    T = t2 - m2 / (n * k)
    return T



