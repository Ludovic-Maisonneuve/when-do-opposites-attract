from generation import *

f_P = [0.5, 0.5]
f_M = [0.01, 0.99]  # Pmut, Pres
f = np.ones((2, 2, 2, 2)) / 16

for P_mother in range(2):
    for P_father in range(2):
        for M_mother in range(2):
            for M_father in range(2):
                f[P_mother, P_father, M_mother, M_father] = f_P[P_mother] * f_P[P_father] * f_M[M_mother] * f_M[
                    M_father]

h_a = 0
h_m = 0
cf = 0.005
cr = 0.005
r = 0
mu = 0.5

n_eq = 100
n_inv = 100

L_delta = [0.02 * i for i in range(51)]
L_h_a = [-1 + 2 * 0.02 * i for i in range(51)]
L_rho = [0.02 * i for i in range(51)]

L_Result = np.zeros((51, 51))

for i_h_a, h_a in enumerate(L_h_a):
    for i_delta, delta in enumerate(L_delta):

        i_MM = -1

        print('delta', delta, flush=True)


        def s_aa(f, h):
            return - (1 - mu) * delta


        def s_ab(f, h):
            return 0


        def s_bb(f, h):
            return - (mu) * delta


        stop = False
        # i_MM = -1
        way = 1

        while stop == False and i_MM < len(L_rho) - 1:

            i_MM += 1

            print('i_rho', i_MM, flush=True)

            rho_MM = L_rho[i_MM]
            gen = generation(f, r, s_aa, s_ab, s_bb, rho_MM, rho_MM, h_a, h_m, cf, cr, 1)

            for t_inv in range(1, n_eq + 1):
                gen = gen.next_gen()

            f_ = gen.f

            if i_MM == 0:
                rho_mm = L_rho[1]

                gen_ = generation(f_, r, s_aa, s_ab, s_bb, rho_mm, rho_MM, h_a, h_m, cf, cr, 1)

                for t_inv in range(1, n_inv + 1):
                    gen_ = gen_.next_gen()

                pm = gen_.p_m

                if pm < 0.01:
                    way = -1

            if i_MM != len(L_rho) - 1 and i_MM != 0:

                rho_mm_1 = L_rho[i_MM - 1]
                rho_mm_2 = L_rho[i_MM + 1]

                gen1 = generation(f_, r, s_aa, s_ab, s_bb, rho_mm_1, rho_MM, h_a, h_m, cf, cr, 1)

                for t_inv in range(1, n_inv + 1):
                    gen1 = gen1.next_gen()

                pm1 = gen1.p_m

                gen2 = generation(f_, r, s_aa, s_ab, s_bb, rho_mm_2, rho_MM, h_a, h_m, cf, cr, 1)

                for t_inv in range(1, n_inv + 1):
                    gen2 = gen2.next_gen()

                pm2 = gen2.p_m

                if pm1 > pm2:
                    way = -1

            if way == -1:
                stop = True

        print(L_rho[i_MM], flush=True)
        L_Result[i_h_a, i_delta] = L_rho[i_MM]

np.save('results/ES_delta_h_a', L_Result)
