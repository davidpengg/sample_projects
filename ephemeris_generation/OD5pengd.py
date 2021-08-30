#OD 5
#David Peng

from odlib import *

#ra, dec in decimalized form
def run_func(ra_input, dec_input, is_print):
     #units
     #distance: au
     #time: gaussian

     #constants
     k = .0172020947
     c = 173.1446
     E = 23.4367487509

     a2 = - 1

     tiny = 10**-15
     max_iteration = 2000

     #2002 KM6
     #Sommers Bausch
     #2458668.7734709373
     #
     
     R_1 = np.array([-1.070812736389900E-1,9.275913635721775E-1,\
                     4.020656126092508E-1])
     R_2 = np.array([-2.073234317681465E-1, 9.133043423994438E-1,\
                     3.958762175519032E-1])
     R_3 = np.array([-4.456361424861880E-1, 8.380253101520228E-1,\
                     3.632408231491512E-1])

     #6-28-19 first
     t1 = 2458662.7769969213
     ra1 = HMStoDeg(17, 39, 44.56)
     dec1 = DMStoDeg(-27, 25, 35.7)

     #7-4-19 first 
     t2 = 2458668.7734709373
     ra2 = ra_input
     dec2 = dec_input

     #7-19-19 first
     t3 = 2458683.7118790857
     ra3 = HMStoDeg(17, 39, 41.28)
     dec3 = DMStoDeg(-21, 48, 20.5)

     #Actual values: taken from ECLIPTIC plane
     #actual r_2
     r_2_right = np.array([1.590863996905743E-01,-1.560550274185269E+00,\
                           -2.362122173499574E-02])

     #positions immediately before and after middle time on 7-4-19
     pos1 = np.array([1.590309637393197E-01,-1.560556339103218E+00,\
                       -2.363054912479409E-02])
     pos2 = np.array([1.591418354799099E-01,-1.560544207827113E+00,\
                       -2.361189432587143E-02])

     #actual r_2_dot
     r_2_dot_right = (pos2 - pos1) / (k * (2458668.776943158\
                                           -2458668.769998714))

     #Sets apparent time to constants
     t1_ap = t1
     t2_ap = t2
     t3_ap = t3


     #Calculates rho hat using ra, dec
     rhohat_1_calc = find_rhohat1(ra1, dec1)
     rhohat_2_calc = find_rhohat1(ra2, dec2)
     rhohat_3_calc = find_rhohat1(ra3, dec3)

     rhohat_1 = np.array(rhohat_1_calc)
     rhohat_2 = np.array(rhohat_2_calc)
     rhohat_3 = np.array(rhohat_3_calc)

     rhohat_calc = [rhohat_1_calc, rhohat_2_calc, rhohat_3_calc]


     #Start of OD4 specific code
     tau0 = k * (t3 - t1)
     tau1 = k * (t1 - t2)
     tau3 = k * (t3 - t2)


     #Calculate Dij
     D11 = dot(cross(R_1, rhohat_2), rhohat_3)
     D12 = dot(cross(R_2, rhohat_2), rhohat_3)
     D13 = dot(cross(R_3, rhohat_2), rhohat_3)

     D21 = dot(cross(rhohat_1, R_1), rhohat_3)
     D22 = dot(cross(rhohat_1, R_2), rhohat_3)
     D23 = dot(cross(rhohat_1, R_3), rhohat_3)

     D31 = dot(rhohat_1, cross(rhohat_2, R_1))
     D32 = dot(rhohat_1, cross(rhohat_2, R_2))
     D33 = dot(rhohat_1, cross(rhohat_2, R_3))
     D0 = dot(rhohat_1, cross(rhohat_2, rhohat_3))


     #First approx for 
     a1 = tau3 / tau0
     a3 = - tau1 / tau0


     #finds f component, given tau_i
     def find_f(tau, r_v, r_dot_v):
          r = mag(r_v)
          O4 = tau**4 / (24 * r**3) \
               * (3 * ((dot(r_dot_v, r_dot_v) / r**2) - (1 / r**3)) \
                  - 15 * (dot(r_v, r_dot_v) / r**2)**2 \
                  + 1 / r**3)
               
          val = 1 - (tau**2 / (2 * r**3)) \
                + ((dot(r_v, r_dot_v) * tau**3) / (2 * r**5)) + O4

          return val


     #finds g component, given tau_i
     def find_g(tau, r_v, r_dot_v):
          r = mag(r_v)
          val = tau - (tau**3 / (6 * r**3))\
                + (dot(r_v, r_dot_v) * tau**4 / (4 * r**5))

          return val
          

     #iterative
     iteration = 0 # initialize iteration counter
     dif_r = dif_r_dot = [100.,100.,100.]


     #sets current r and rdot approximations to 'old' values
     r_2_old = np.array([1000, 1000, 1000])
     r_2_dot_old = np.array([1000, 1000, 1000])


     #performs time correction
     def time_corr(t, rho):
          new_t = t - rho / c
          return new_t


     while iteration <= max_iteration: # iterates while under max iteration

          rho_1 = (a1 * D11 + a2 * D12 + a3 * D13 ) / (a1 * D0)
          rho_2 = (a1 * D21 + a2 * D22 + a3 * D23 ) / (a2 * D0)
          rho_3 = (a1 * D31 + a2 * D32 + a3 * D33 ) / (a3 * D0)

          #time correction
          t1 = time_corr(t1_ap, rho_1)
          t2 = time_corr(t2_ap, rho_2)
          t3 = time_corr(t3_ap, rho_3)

          #recalculates tau
          tau0 = k * (t3 - t1)
          tau1 = k * (t1 - t2)
          tau3 = k * (t3 - t2)
          
          r_1 = rho_1 * rhohat_1 - R_1
          r_2 = rho_2 * rhohat_2 - R_2
          r_3 = rho_3 * rhohat_3 - R_3

          if iteration == 0:
               r_2_dot = (r_3 - r_1) / tau0

          g1 = find_g(tau1, r_2_old, r_2_dot_old)
          g3 = find_g(tau3, r_2_old, r_2_dot_old)

          f1 = find_f(tau1, r_2_old, r_2_dot_old)
          f3 = find_f(tau3, r_2_old, r_2_dot_old)
          
          a1 = g3 / (f1 * g3 - f3 * g1)
          a3 = g1 / (f3 * g1 - f1 * g3)

          b1 = f3 / (f3 * g1 - f1 * g3)
          b3 = f1 / (f1 * g3 - f3 * g1)

          r_2 = a1 * r_1 + a3 * r_3
          r_2_dot = b1 * r_1 + b3 * r_3

          if iteration != 0:
               dif_r = [abs(r_2[i] - r_2_old[i]) for i in range(3)]
               dif_r_dot = [abs(r_2_dot[i] - r_2_dot_old[i]) for i in range(3)]

          #Determines if difference is tiny enough
          if all([dif_r[i] < tiny and dif_r_dot[i] < tiny for i in range(3)]):
               break

          r_2_old = r_2
          r_2_dot_old = r_2_dot
          
          iteration += 1


     r_2 = coor_rot(E, r_2, 'x') #rotate to ecliptic
     r_2_dot = coor_rot(E, r_2_dot, 'x')

     r_v = r_2
     r_dot_v = r_2_dot
     
     r = mag(r_v)
     r_dot = mag(r_dot_v)

     h_v = find_h(r_v, r_dot_v) #h is vector
     h = mag(h_v)

     vv = v_squared(r_dot_v) #vv is v^2, a constant

     a = find_a(r, vv)
     print('a',a)
     ee = find_e(h,a)

     i = find_i(h_v)
             
     omega = find_omega(h_v, i)

     U = find_U(omega, r, i, r_v)
     V = find_V(a, ee, h_v, r_v, r_dot_v)
     w = find_W(U, V)

     E = find_E(ee, r, a, V)

     M = find_M(E, ee)
     M += 2 * pi

     def print_it(val, real):
          print('{:<42}'.format(val), real)

     if is_print == True:
          
          ro = 11
          print_it('Expected r', r_2_right)
          print_it('Calculated r', [round(i, ro) for i in r_2])
          print_it('Percent difference for each component (%)',\
                   [round(percent_error(r_2_right[i], r_2[i]),ro)\
                    for i in range(3)])
          print_it('Expected r dot', r_2_dot_right)
          print_it('Calculated r dot', [round(i, ro) for i in r_2_dot])
          print_it('Percent difference for each component (%)',\
                   [round(percent_error(r_2_dot_right[i], r_2_dot[i]),ro)\
                    for i in range(3)])
          print()

          real_a = 2.640264641359547E+00
          real_e = 4.059093870368781E-01
          real_i = 9.527029779017093E+00
          real_o = 2.809685857380750E+02
          real_w = 3.567203322934803E+02
          real_M = 3.592508065642508E+02
          
          print('Semimajor Axis')
          print('Expected a:', real_a)
          print('Calculated a:', a)
          print('Percent error:', percent_error(real_a, a), '%')

          print()

          #EC
          print('Eccentricity')
          print('Expected e:', real_e)
          print('Calculated e:', ee)
          print('Percent error:', percent_error(real_e, ee), '%')

          print()

          #IN
          print('Inclination')
          print('Expected i:', real_i)
          print('Calculated i:', i * 180 / pi)
          print('Percent error:', percent_error(real_i, i * 180 / pi), '%')

          print()

          #OM
          print('Longitude of Ascending Node')
          print('Expected omega:', real_o)
          print('Calculated omega:', omega * 180 / pi)
          print('Percent error:', percent_error(real_o, \
                                                omega * 180 / pi), '%')

          print()

          #W
          print('Argument of Perifocus')
          print('Expected w:', real_w)
          print('Calculated w:', w * 180 / pi)
          print('Percent error:', percent_error(real_w, \
                                                w * 180 / pi), '%')

          print()

##          #n/a
##          print('Eccentric Anomaly')
##          print('Expected E: ?')
##          print('Calculated E:', E * 180 / pi)
##          print('Percent error: ?') #percent_error(real_E, \
##                                    #            E * 180 / pi), '%')

          print()

          #L
          print('Mean Anomaly')
          print('Expected M:', real_M)
          print('Calculated M:', M * 180 / pi)
          print('Percent error:', percent_error(real_M, \
                                                M * 180 / pi), '%')
     return [a, ee, i, omega, w, E, M]

run_func(HMStoDeg(17, 38, 19.53), DMStoDeg(-25, 44, 27.6), True)
