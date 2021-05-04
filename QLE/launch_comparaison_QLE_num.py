import matplotlib.pyplot as plt

from generation import *

f_P = [0.5, 0.5]
f_M = [0.01, 0.99]
f = np.ones((2, 2, 2, 2)) / 16

for P_mother in range(2):
    for P_father in range(2):
        for M_mother in range(2):
            for M_father in range(2):
                f[P_mother, P_father, M_mother, M_father] = f_P[P_mother] * f_P[P_father] * f_M[M_mother] * f_M[
                    M_father]


def s_aa(f, h):
    return 0


def s_ab(f, h):
    return 0.01


def s_bb(f, h):
    return 0


rho_mm = 0.01
rho_MM = 0
h_a = 1
h_m = -1
cf = 0
cr = 0
r = 0.1

gen = generation(f, r, s_aa, s_ab, s_bb, rho_mm, rho_MM, h_a, h_m, cf, cr, 1)
for i in range(100):
    gen = gen.next_gen()

f = gen.f

L_s_aa = [0.001 * i for i in range(1001)]
L_rho_mm = [0.001 * i for i in range(1001)]
L_cf = [0.001 * i for i in range(501)]
L_cr = [0.001 * i for i in range(501)]

## natural selection

L_s_aa_Delta_a_num = []
L_s_aa_Delta_a_QLE = []
L_s_aa_Delta_a_err = []

L_s_aa_Delta_m_num = []
L_s_aa_Delta_m_QLE = []
L_s_aa_Delta_m_err = []

for s_aa_ in L_s_aa:
    def s_aa(f, h):
        return s_aa_


    gen = generation(f, r, s_aa, s_ab, s_bb, rho_mm, rho_MM, h_a, h_m, cf, cr, 1)
    Delta_p_a_QLE, Delta_p_m_QLE = gen.QLE_Delta_a_and_m()
    p_a_ = gen.p_a
    p_m_ = gen.p_m
    gen = gen.next_gen()
    p_a = gen.p_a
    p_m = gen.p_m

    L_s_aa_Delta_a_num.append(p_a - p_a_)
    L_s_aa_Delta_a_QLE.append(Delta_p_a_QLE)
    L_s_aa_Delta_a_err.append(np.abs((Delta_p_a_QLE - (p_a - p_a_)) / (p_a - p_a_)))

    L_s_aa_Delta_m_num.append(p_m - p_m_)
    L_s_aa_Delta_m_QLE.append(Delta_p_m_QLE)
    L_s_aa_Delta_m_err.append(np.abs((Delta_p_m_QLE - (p_m - p_m_)) / (p_m - p_m_)))

plt.figure()
plt.plot(L_s_aa, L_s_aa_Delta_a_num, label=r'numerical $\Delta p_a$')
plt.plot(L_s_aa, L_s_aa_Delta_a_QLE, label=r'QLE $\Delta p_a$')
plt.legend()
plt.xlabel(r'$S_{aa}$')
plt.ylabel(r'$\Delta p_a$')

plt.figure()
plt.plot(L_s_aa, L_s_aa_Delta_a_err)
plt.xlabel(r'$S_{aa}$')
plt.ylabel(r'$\epsilon$')

plt.figure()
plt.plot(L_s_aa, L_s_aa_Delta_m_num, label=r'numerical $\Delta p_m$')
plt.plot(L_s_aa, L_s_aa_Delta_m_QLE, label=r'QLE $\Delta p_m$')
plt.legend()
plt.xlabel(r'$S_{aa}$')
plt.ylabel(r'$\Delta p_m$')

plt.figure()
plt.plot(L_s_aa, L_s_aa_Delta_m_err)
plt.xlabel(r'$S_{aa}$')
plt.ylabel(r'$\epsilon$')


def s_aa(f, h):
    return 0.01


## preference

L_rho_mm_Delta_a_num = []
L_rho_mm_Delta_a_QLE = []
L_rho_mm_Delta_a_err = []

L_rho_mm_Delta_m_num = []
L_rho_mm_Delta_m_QLE = []
L_rho_mm_Delta_m_err = []

for rho_mm in L_rho_mm:
    gen = generation(f, r, s_aa, s_ab, s_bb, rho_mm, rho_MM, h_a, h_m, cf, cr, 1)
    Delta_p_a_QLE, Delta_p_m_QLE = gen.QLE_Delta_a_and_m()
    p_a_ = gen.p_a
    p_m_ = gen.p_m
    gen = gen.next_gen()
    p_a = gen.p_a
    p_m = gen.p_m

    L_rho_mm_Delta_a_num.append(p_a - p_a_)
    L_rho_mm_Delta_a_QLE.append(Delta_p_a_QLE)
    L_rho_mm_Delta_a_err.append(np.abs((Delta_p_a_QLE - (p_a - p_a_)) / (p_a - p_a_)))

    L_rho_mm_Delta_m_num.append(p_m - p_m_)
    L_rho_mm_Delta_m_QLE.append(Delta_p_m_QLE)
    L_rho_mm_Delta_m_err.append(np.abs((Delta_p_m_QLE - (p_m - p_m_)) / (p_m - p_m_)))

plt.figure()
plt.plot(L_s_aa, L_rho_mm_Delta_a_num, label=r'numerical $\Delta p_a$')
plt.plot(L_s_aa, L_rho_mm_Delta_a_QLE, label=r'QLE $\Delta p_a$')
plt.legend()
plt.xlabel(r'$\rho_{mm}$')
plt.ylabel(r'$\Delta p_a$')

plt.figure()
plt.plot(L_s_aa, L_rho_mm_Delta_a_err, label=r'numerical $\Delta p_a$')
plt.xlabel(r'$\rho_{mm}$')
plt.ylabel(r'$\epsilon$')

plt.figure()
plt.plot(L_s_aa, L_rho_mm_Delta_m_num, label=r'numerical $\Delta p_m$')
plt.plot(L_s_aa, L_rho_mm_Delta_m_QLE, label=r'QLE $\Delta p_m$')
plt.legend()
plt.xlabel(r'$\rho_{mm}$')
plt.ylabel(r'$\Delta p_m$')

plt.figure()
plt.plot(L_s_aa, L_rho_mm_Delta_m_err, label=r'numerical $\Delta p_m$')
plt.xlabel(r'$\rho_{mm}$')
plt.ylabel(r'$\epsilon$')

rho_mm = 0.01

## fixed cost of choosiness

L_cf_Delta_m_num = []
L_cf_Delta_m_QLE = []
L_cf_Delta_m_err = []

for cf in L_cf:
    gen = generation(f, r, s_aa, s_ab, s_bb, rho_mm, rho_MM, h_a, h_m, cf, cr, 1)
    Delta_p_a_QLE, Delta_p_m_QLE = gen.QLE_Delta_a_and_m()
    p_a_ = gen.p_a
    p_m_ = gen.p_m
    gen = gen.next_gen()
    p_a = gen.p_a
    p_m = gen.p_m

    L_cf_Delta_m_num.append(p_m - p_m_)
    L_cf_Delta_m_QLE.append(Delta_p_m_QLE)
    L_cf_Delta_m_err.append(np.abs((Delta_p_m_QLE - (p_m - p_m_)) / (p_m - p_m_)))

plt.figure()
plt.plot(L_cf, L_cf_Delta_m_num, label=r'numerical $\Delta p_m$')
plt.plot(L_cf, L_cf_Delta_m_QLE, label=r'QLE $\Delta p_m$')
plt.legend()
plt.xlabel(r'$c_f$')
plt.ylabel(r'$\Delta p_m$')

plt.figure()
plt.plot(L_cf, L_cf_Delta_m_err, label=r'numerical $\Delta p_m$')
plt.xlabel(r'$c_f$')
plt.ylabel(r'$\epsilon$')

cf = 0

## relative cost of choosiness

L_cr_Delta_a_num = []
L_cr_Delta_a_QLE = []
L_cr_Delta_a_err = []

L_cr_Delta_m_num = []
L_cr_Delta_m_QLE = []
L_cr_Delta_m_err = []

for cr in L_cr:
    gen = generation(f, r, s_aa, s_ab, s_bb, rho_mm, rho_MM, h_a, h_m, cf, cr, 1)
    Delta_p_a_QLE, Delta_p_m_QLE = gen.QLE_Delta_a_and_m()
    p_a_ = gen.p_a
    p_m_ = gen.p_m
    gen = gen.next_gen()
    p_a = gen.p_a
    p_m = gen.p_m

    L_cr_Delta_a_num.append(p_a - p_a_)
    L_cr_Delta_a_QLE.append(Delta_p_a_QLE)
    L_cr_Delta_a_err.append(np.abs((Delta_p_a_QLE - (p_a - p_a_)) / (p_a - p_a_)))

    L_cr_Delta_m_num.append(p_m - p_m_)
    L_cr_Delta_m_QLE.append(Delta_p_m_QLE)
    L_cr_Delta_m_err.append(np.abs((Delta_p_m_QLE - (p_m - p_m_)) / (p_m - p_m_)))

plt.figure()
plt.plot(L_cf, L_cr_Delta_a_num, label=r'numerical $\Delta p_a$')
plt.plot(L_cf, L_cr_Delta_a_QLE, label=r'QLE $\Delta p_a$')
plt.legend()
plt.xlabel(r'$c_r$')
plt.ylabel(r'$\Delta p_a$')

plt.figure()
plt.plot(L_cf, L_cr_Delta_a_err, label=r'numerical $\Delta p_a$')
plt.xlabel(r'$c_r$')
plt.ylabel(r'$\epsilon$')
plt.show()

plt.figure()
plt.plot(L_cf, L_cr_Delta_m_num, label=r'numerical $\Delta p_m$')
plt.plot(L_cf, L_cr_Delta_m_QLE, label=r'QLE $\Delta p_m$')
plt.legend()
plt.xlabel(r'$c_r$')
plt.ylabel(r'$\Delta p_m$')

plt.figure()
plt.plot(L_cf, L_cr_Delta_m_err, label=r'numerical $\Delta p_m$')
plt.xlabel(r'$c_r$')
plt.ylabel(r'$\epsilon$')
plt.show()
