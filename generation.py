from copy import *

import numpy as np


class generation:

    def __init__(self, f, r, s_aa, s_ab, s_bb, rho_mm, rho_MM, h_a, h_m, cf, cr, g):
        self.f = deepcopy(f)  # np.array containing the frequency of every genotype
        self.r = r  # recombination rate between the locus C and the locus P
        self.s_aa = s_aa  # fitness advantage due to natural selection of individual of gentope aa at locus C. s_aa is a function of f
        self.s_ab = s_ab  # fitness advantage due to natural selection of individual of gentope ab at locus C. s_aa is a function of f
        self.s_bb = s_bb  # fitness advantage due to natural selection of individual of gentope bb at locus C. s_aa is a function of f
        self.rho_mm = rho_mm  # strenght of assortative mating of female with genotype mm at locus P
        self.rho_mM = (1 - h_m) / 2 * rho_MM + (
                    1 + h_m) / 2 * rho_mm  # strenght of assortative mating of female with genotype mM at locus P
        self.rho_MM = rho_MM  # strenght of assortative mating of female with genotype MM at locus P
        self.h_a = h_a  # dominance coefficiant at locus C
        self.h_M = h_m  # dominance coefficiant at locus P
        self.cf = cf  # fixed cost of choosiness
        self.cr = cr  # relative cost of choosiness
        self.g = g  # generation number

        self.f_aa_mm = f[0, 0, 0, 0]
        self.f_aa_mM = f[0, 0, 0, 1] + f[0, 0, 1, 0]
        self.f_aa_MM = f[0, 0, 1, 1]
        self.f_ab_mm = f[1, 0, 0, 0] + f[0, 1, 0, 0]
        self.f_ab_mM = f[1, 0, 1, 0] + f[0, 1, 0, 1]
        self.f_ab_Mm = f[1, 0, 0, 1] + f[0, 1, 1, 0]
        self.f_ab_MM = f[1, 0, 1, 1] + f[0, 1, 1, 1]
        self.f_bb_mm = f[1, 1, 0, 0]
        self.f_bb_mM = f[1, 1, 1, 0] + f[1, 1, 0, 1]
        self.f_bb_MM = f[1, 1, 1, 1]

        self.p_a = self.f_aa_mm + self.f_aa_mM + self.f_aa_MM + (
                self.f_ab_mm + self.f_ab_mM + self.f_ab_Mm + self.f_ab_MM) / 2
        self.p_b = 1 - self.p_a
        self.p_m = self.f_aa_mm + self.f_ab_mm + self.f_bb_mm + (
                self.f_aa_mM + self.f_ab_mM + self.f_ab_Mm + self.f_bb_mM) / 2
        self.p_M = 1 - self.p_m

        self.p_aa = self.f_aa_mm + self.f_aa_mM + self.f_aa_MM
        self.p_ab = self.f_ab_mm + self.f_ab_mM + self.f_ab_Mm + self.f_ab_MM
        self.p_bb = self.f_bb_mm + self.f_bb_mM + self.f_bb_MM

        self.p_mm = self.f_aa_mm + self.f_ab_mm + self.f_bb_mm
        self.p_mM = self.f_aa_mM + self.f_ab_mM + self.f_ab_Mm + self.f_bb_mM
        self.p_MM = self.f_aa_MM + self.f_ab_MM + self.f_bb_MM

        self.p_am = self.f_aa_mm + (self.f_aa_mM + self.f_ab_mM + self.f_ab_mm) / 2
        self.p_a_m = self.f_aa_mm + (self.f_aa_mM + self.f_ab_Mm + self.f_ab_mm) / 2
        self.p_am_a = self.f_aa_mm + self.f_aa_mM / 2

        self.D_a = self.p_aa - self.p_a ** 2
        self.D_am = self.p_am - self.p_a * self.p_m
        self.D_a_m = self.p_a_m - self.p_a * self.p_m
        self.D_am_a = self.p_am_a - self.p_a * (self.D_am + self.D_a_m) - self.p_m * self.D_a - self.p_a ** 2 * self.p_m

        Pref = np.zeros((2, 2, 2, 2, 2, 2))
        for P_mother_choosy in range(2):
            for P_father_choosy in range(2):
                if P_mother_choosy == 0 and P_father_choosy == 0:
                    rho = rho_mm
                elif P_mother_choosy == 1 and P_father_choosy == 1:
                    rho = rho_MM
                else:
                    rho = self.rho_mM
                for C_mother_choosy in range(2):
                    for C_father_choosy in range(2):
                        for C_mother_chosen in range(2):
                            for C_father_chosen in range(2):
                                if C_mother_choosy == C_father_choosy:  # the choosy individual is homozygote
                                    if C_mother_chosen == C_father_chosen:  # the chosen individual is homozygote
                                        if C_mother_choosy == C_mother_chosen:  # the choosy and the chosen individuals have the same genotype
                                            Pref[
                                                C_mother_choosy, C_father_choosy, P_mother_choosy, P_father_choosy, C_mother_chosen, C_father_chosen] = 1 - rho
                                        else:  # the choosy and the chosen individuals do not have the same genotype
                                            Pref[
                                                C_mother_choosy, C_father_choosy, P_mother_choosy, P_father_choosy, C_mother_chosen, C_father_chosen] = 1
                                    else:  # the chosen individual is heterozygote
                                        if C_mother_choosy == 0:
                                            Pref[
                                                C_mother_choosy, C_father_choosy, P_mother_choosy, P_father_choosy, C_mother_chosen, C_father_chosen] = 1 - (
                                                    1 + h_a) / 2 * rho
                                        else:
                                            Pref[
                                                C_mother_choosy, C_father_choosy, P_mother_choosy, P_father_choosy, C_mother_chosen, C_father_chosen] = 1 - (
                                                    1 - h_a) / 2 * rho
                                else:  # the choosy individual is heterozygote
                                    if C_mother_chosen == C_father_chosen:  # the chosen individual is homozygote
                                        if C_mother_chosen == 0:
                                            Pref[
                                                C_mother_choosy, C_father_choosy, P_mother_choosy, P_father_choosy, C_mother_chosen, C_father_chosen] = 1 - (
                                                    1 + h_a) / 2 * rho
                                        else:
                                            Pref[
                                                C_mother_choosy, C_father_choosy, P_mother_choosy, P_father_choosy, C_mother_chosen, C_father_chosen] = 1 - (
                                                    1 - h_a) / 2 * rho
                                    else:
                                        Pref[
                                            C_mother_choosy, C_father_choosy, P_mother_choosy, P_father_choosy, C_mother_chosen, C_father_chosen] = 1 - rho
        self.Pref = Pref

    def W(self):
        W = np.zeros((2, 2, 2, 2))

        for P_mother in range(2):
            for P_father in range(2):
                for C_mother in range(2):
                    for C_father in range(2):
                        if C_mother == C_father:
                            if C_mother == 0:
                                W[C_mother, C_father, P_mother, P_father] = (1 + self.s_aa(self.f, self.h_a))
                            else:
                                W[C_mother, C_father, P_mother, P_father] = (1 + self.s_bb(self.f, self.h_a))
                        else:
                            W[C_mother, C_father, P_mother, P_father] = (1 + self.s_ab(self.f, self.h_a))
        return W

    def T_F(self):
        T = np.zeros((2, 2, 2, 2))
        F = np.zeros((2, 2, 2, 2))
        for C_mother in range(2):
            for C_father in range(2):
                for P_mother in range(2):
                    for P_father in range(2):
                        if P_mother == 0 and P_father == 0:
                            rho = self.rho_mm
                        elif P_mother == 1 and P_father == 1:
                            rho = self.rho_MM
                        else:
                            rho = self.rho_mM
                        T_ = 0
                        for C_mother_chosen in range(2):
                            for C_father_chosen in range(2):
                                T_ += self.Pref[
                                          C_mother, C_father, P_mother, P_father, C_mother_chosen, C_father_chosen] * (
                                              self.f[C_mother_chosen, C_father_chosen, 0, 0] + self.f[
                                          C_mother_chosen, C_father_chosen, 0, 1] + self.f[
                                                  C_mother_chosen, C_father_chosen, 1, 0] + self.f[
                                                  C_mother_chosen, C_father_chosen, 1, 1])
                        T[C_mother, C_father, P_mother, P_father] = T_
                        F[C_mother, C_father, P_mother, P_father] = ((1 - self.cr) + self.cr * T_) * (
                                1 - self.cf * rho)
        return T, F

    def next_step_selection(self):
        W = self.W()
        Pref = self.Pref
        W_average = 0
        for C_mother in range(2):
            for C_father in range(2):
                for P_mother in range(2):
                    for P_father in range(2):
                        W_average += W[C_mother, C_father, P_mother, P_father] * self.f[
                            C_mother, C_father, P_mother, P_father]
        f_ = np.zeros((2, 2, 2, 2))
        for C_mother in range(2):
            for C_father in range(2):
                for P_mother in range(2):
                    for P_father in range(2):
                        f_[C_mother, C_father, P_mother, P_father] = W[
                                                                         C_mother, C_father, P_mother, P_father] / W_average * \
                                                                     self.f[C_mother, C_father, P_mother, P_father]
        return generation(f_ / np.sum(f_), self.r, self.s_aa, self.s_ab, self.s_bb, self.rho_mm,
                          self.rho_MM, self.h_a, self.h_M,
                          self.cf, self.cr, self.g + 1 / 2)

    def next_gen_mating(self):
        f_ = np.zeros((2, 2, 2, 2))
        T, F = self.T_F()
        Pref = self.Pref
        F_average = 0
        for C_mother in range(2):
            for C_father in range(2):
                for P_mother in range(2):
                    for P_father in range(2):
                        F_average += F[C_mother, C_father, P_mother, P_father] * self.f[
                            C_mother, C_father, P_mother, P_father]
        for C_mother_ in range(2):
            for C_father_ in range(2):
                for P_mother_ in range(2):
                    for P_father_ in range(2):
                        for C_mother_female in range(2):
                            for C_father_female in range(2):
                                for P_mother_female in range(2):
                                    for P_father_female in range(2):
                                        for C_mother_male in range(2):
                                            for C_father_male in range(2):
                                                for P_mother_male in range(2):
                                                    for P_father_male in range(2):
                                                        f_[C_mother_, C_father_, P_mother_, P_father_] += self.f[
                                                                                                              C_mother_female, C_father_female, P_mother_female, P_father_female] * \
                                                                                                          F[
                                                                                                              C_mother_female, C_father_female, P_mother_female, P_father_female] / F_average * \
                                                                                                          coef(
                                                                                                              C_mother_,
                                                                                                              C_father_,
                                                                                                              P_mother_,
                                                                                                              P_father_,
                                                                                                              C_mother_female,
                                                                                                              C_father_female,
                                                                                                              P_mother_female,
                                                                                                              P_father_female,
                                                                                                              C_mother_male,
                                                                                                              C_father_male,
                                                                                                              P_mother_male,
                                                                                                              P_father_male,
                                                                                                              self.r) * \
                                                                                                          Pref[
                                                                                                              C_mother_female, C_father_female, P_mother_female, P_father_female, C_mother_male, C_father_male] * \
                                                                                                          self.f[
                                                                                                              C_mother_male, C_father_male, P_mother_male, P_father_male] / \
                                                                                                          T[
                                                                                                              C_mother_female, C_father_female, P_mother_female, P_father_female]
        return generation(f_ / np.sum(f_), self.r, self.s_aa, self.s_ab, self.s_bb, self.rho_mm,
                          self.rho_MM, self.h_a, self.h_M,
                          self.cf, self.cr, self.g + 1 / 2)

    def next_gen(self):
        # print('before natural selection', 'paa',self.paa,'pab',self.pab,'pbb',self.pbb)
        gen = self.next_step_selection()
        # print('after natural selection', 'paa',gen.paa,'pab',gen.pab,'pbb',gen.pbb)
        return gen.next_gen_mating()

    def QLE_Delta_a_and_m(self):
        p_a = self.p_a
        p_b = self.p_b
        p_M = self.p_M
        p_m = self.p_m
        rho_ = p_M ** 2 * self.rho_MM + 2 * p_M * p_m * self.rho_mM + p_m ** 2 * self.rho_mm
        Delta_a = p_a * p_b * (p_a * (self.s_aa(self.f, self.h_a) - self.s_ab(self.f, self.h_a)) + p_b * (
                self.s_ab(self.f, self.h_a) - self.s_bb(self.f, self.h_a))) + rho_ * (1 + self.cr) * p_a * p_b * (
                          (p_b ** 4 - p_a ** 4) / 4 + self.h_a * (p_a ** 2 + p_b ** 2 - 4 * p_a * p_b) / 4)

        H_ns = self.s_aa(self.f, self.h_a) + self.s_bb(self.f, self.h_a) - 2 * self.s_ab(self.f, self.h_a)
        H1 = (1 - self.h_a) / 2 * rho_
        H2 = rho_
        H3 = (1 + self.h_a) / 2 * rho_
        H_ss = p_b ** 2 * H1 + 2 * p_a * p_b * H2 + p_a ** 2 * H3
        H_ss = H_ss / 2
        H_tot = H_ns + H_ss
        Delta_rho = p_m * (self.rho_mm - self.rho_mM) + p_M * (self.rho_mM - self.rho_MM)
        mate0 = p_a ** 4 + 4 * p_a ** 2 * p_b ** 2 + p_b ** 4
        matea = 4 * p_a ** 3 * p_b
        mateb = 4 * p_a * p_b ** 3
        mate1 = matea + mateb
        Delta_m = - Delta_rho / 2 * p_m * p_M * (
                    self.cf + self.cr * (mate0 + mate1 / 2 + self.h_a / 2 * (matea - mateb))) + self.D_am_a * H_tot + (
                              self.D_am + self.D_a_m) * Delta_a / (
                          p_a * p_b)
        return Delta_a, Delta_m


def coef_haplotype(C_offspring, P_offspring, C_mother_parent, C_father_parent, P_mother_parent, P_father_parent, r):
    if C_mother_parent == C_father_parent:  # the parent is homozygote ate locus P
        if C_offspring != C_mother_parent:
            return 0
        else:
            if P_mother_parent == P_father_parent:  # the parent is homozygote ate locus M
                if P_offspring != P_mother_parent:
                    return 0
                else:
                    return 1
            else:  # the parent is heterozygote ate locus M
                return 1 / 2
    else:  # the parent is heterozygote ate locus P
        if P_mother_parent == P_father_parent:  # the parent is homozygote ate locus M
            if P_offspring != P_mother_parent:
                return 0
            else:
                return 1 / 2
        else:  # the parent is heterozygote ate locus M
            if C_offspring == C_mother_parent and P_offspring == P_mother_parent:
                return (1 - r) / 2
            elif C_offspring == C_father_parent and P_offspring == P_father_parent:
                return (1 - r) / 2
            elif C_offspring == C_mother_parent and P_offspring == P_mother_parent:
                return r / 2
            else:
                return r / 2


def coef(C_mother, C_father, P_mother, P_father, C_mother_female, C_father_female, P_mother_female, P_father_female,
         C_mother_male, C_father_male, P_mother_male, P_father_male, r):
    return coef_haplotype(C_mother, P_mother, C_mother_female, C_father_female, P_mother_female, P_father_female,
                          r) * coef_haplotype(C_father, P_father, C_mother_male, C_father_male, P_mother_male,
                                              P_father_male, r)
