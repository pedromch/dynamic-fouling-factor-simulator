from numpy import pi, sum, exp, log as ln
from scipy.optimize import minimize, root


def calculate(Ds, Dte, Npt, Rtp, layout, Ltube, esp, ktube, Nc, fluid_inside_tubes, data_to_fit, tempo_final=None):
    Ltp = Rtp * Dte
    Ntt = get_total_number_of_tubes(Ds, Ltp, Npt, layout)
    At = Ntt * pi * Dte * Ltube
    Dti = Dte - 2 * esp
    Ntp = Ntt/Npt
    Deq = get_equivalent_diameter(Dte, Ltp, layout)
    Lbc = Ltube / (Nc + 1)
    Ac = get_shell_flow_area(Dte, Ds, Ltp, Lbc)
    
    def get_Rft(mt, rhot, Cpt, mut, kt, Tti, Tto, ms, rhos, Cps, mus, ks, Tsi, Tso):
        if fluid_inside_tubes == "frio":
            Tci, Tco = Tti, Tto
            Thi, Tho = Tsi, Tso
            mh = ms
            Cph = Cps
            operation_inside_tubes = "aquecimento" # check
        elif fluid_inside_tubes == "quente":
            Tci, Tco = Tsi, Tso
            Thi, Tho = Tti, Tto
            mh = mt
            Cph = Cpt
            operation_inside_tubes = "resfriamento" # check
        
        vt = mt / (rhot * Ntp * pi * Dti**2 / 4)
        Ret = get_reynolds(vt, rhot, Dti, mut)
        Prt = get_prandtl(Cpt, mut, kt)
        ht = get_nu_dittus_boelter(Ret, Prt, operation_inside_tubes) * kt / Dti
        
        vs = ms / (rhos * Ac)
        Res = get_reynolds(vs, rhos, Deq, mus)
        Prs = get_prandtl(Cps, mus, ks)
        hs = get_nu_kern(Res, Prs, mus, mus) * ks / Deq
        
        Uc = get_clean_global_heat_transfer_coefficient(ht, hs, Dti, Dte, ktube)
        
        Qd = mh * Cph * (Thi - Tho)
        deltaTLMd = get_logarithmic_mean_temperature_difference(Thi, Tho, Tci, Tco)
        F = get_mean_temperature_difference_correction_factor(Thi, Tho, Tci, Tco)
        Ud = get_dirty_global_heat_transfer_coefficient(Qd, At, deltaTLMd, F)
        
        return 1/Ud - 1/Uc
    
    def get_Rft_from_row(row):
        if fluid_inside_tubes == "frio":
            return get_Rft(row.mc, row.rhoc, row.Cpc, row.muc, row.kc, row.Tci, row.Tco, row.mh, row.rhoh, row.Cph, row.muh, row.kh, row.Thi, row.Tho)
        elif fluid_inside_tubes == "quente":
            return get_Rft(*row)
        
    data_to_fit["Rft"] = data_to_fit.apply(lambda row: get_Rft_from_row(row), axis=1)
    
    Rft_array = data_to_fit["Rft"].to_numpy()
    t_array = data_to_fit.index.to_numpy()
    
    def fobj(params):
        Rftinf, Rft0, S = params
        sum_to_minimize = sum((Rft_array - get_fitted_total_fouling_factor(t_array, Rftinf, Rft0, S))**2 / Rft_array)
        return sum_to_minimize
    
    guess = [Rft_array[-1], Rft_array[0], 1]
    result = minimize(fobj, guess)
    Rftinf, Rft0, S = result.x
    
    fobj_last_value = fobj([Rftinf, Rft0, S])
    
    data_to_fit["Rft_fitted"] = get_fitted_total_fouling_factor(data_to_fit.index.to_numpy(), Rftinf, Rft0, S)
    
    if (tempo_final is not None):
        data_to_extrapolate = data_to_fit.iloc[-1].to_numpy()
        if fluid_inside_tubes == "frio":
            ms, rhos, Cps, mus, ks, Tsi, Tso, mt, rhot, Cpt, mut, kt, Tti, Tto, _, _ = data_to_extrapolate
            Tci, Thi = Tti, Tsi
            Tco_guess, Tho_guess = Tto, Tso
            mh, Cph, mc, Cpc = ms, Cps, mt, Cpt
            operation_inside_tubes = "aquecimento" # check
        elif fluid_inside_tubes == "quente":
            mt, rhot, Cpt, mut, kt, Tti, Tto, ms, rhos, Cps, mus, ks, Tsi, Tso, _, _ = data_to_extrapolate
            Tci, Thi = Tsi, Tti
            Tco_guess, Tho_guess = Tso, Tto
            mh, Cph, mc, Cpc = mt, Cpt, ms, Cps
            operation_inside_tubes = "resfriamento" # check
            
        vt = mt / (rhot * Ntp * pi * Dti**2 / 4)
        Ret = get_reynolds(vt, rhot, Dti, mut)
        Prt = get_prandtl(Cpt, mut, kt)
        ht = get_nu_dittus_boelter(Ret, Prt, operation_inside_tubes) * kt / Dti
        
        vs = ms / (rhos * Ac)
        Res = get_reynolds(vs, rhos, Deq, mus)
        Prs = get_prandtl(Cps, mus, ks)
        hs = get_nu_kern(Res, Prs, mus, mus) * ks / Deq
        
        Uc = get_clean_global_heat_transfer_coefficient(ht, hs, Dti, Dte, ktube)
        
        def fobj(params, t):
            Q, Tho, Tco = params
            Rft = get_fitted_total_fouling_factor(t, Rftinf, Rft0, S)
            Ud = 1 / (1/Uc + Rft)
            deltaTLM = get_logarithmic_mean_temperature_difference(Thi, Tho, Tci, Tco)
            F = get_mean_temperature_difference_correction_factor(Thi, Tho, Tci, Tco)
            return [
                Q - Ud * At * deltaTLM * F,
                Q - mh * Cph * (Thi - Tho),
                Q - mc * Cpc * (Tco - Tci)
            ]
        
        Qguess = mh * Cph * (Thi - Tho_guess)
        guess = [Qguess, Tho_guess, Tco_guess]
        def get_Q_Tho_Tco(t):
            return root(fobj, guess, args=t).x
        
        t_list = list(range(0, tempo_final))
        Q_estimated, Tho_estimated, Tco_estimated = zip(*map(get_Q_Tho_Tco, t_list))
        Ud_estimated = [1 / (1/Uc + Rft) for Rft in [get_fitted_total_fouling_factor(t, Rftinf, Rft0, S) for t in t_list]]
        
        return Rftinf, Rft0, S, fobj_last_value, data_to_fit, t_list, Q_estimated, Tho_estimated, Tco_estimated, Ud_estimated
    
    return Rftinf, Rft0, S, fobj_last_value, data_to_fit


def get_nu_dittus_boelter(Re, Pr, operation):
    n_options = {
        "aquecimento": 0.4,
        "resfriamento": 0.3
    }
    
    n = n_options[operation]

    return 0.023 * Re**0.8 * Pr**n


def get_nu_kern(Re, Pr, mu, muw):
    return 0.36 * Re**0.55 * Pr**(1/3) * (mu/muw)**0.14


def get_equivalent_diameter(Dte, Ltp, layout):
    f_options = {
        "quadrado": 4,
        "triangular": 3.46
    }

    f = f_options[layout]

    return f * Ltp**2 / (pi * Dte) - Dte


def get_shell_flow_area(Dte, Ds, Ltp, Lbc):
    return Ds * (Ltp - Dte) * Lbc / Ltp


def get_clean_global_heat_transfer_coefficient(hi, he, Dti, Dte, kt):
    return 1 / (1/hi * Dte/Dti + Dte * ln(Dte/Dti) / (2 * kt) + 1/he)


def get_logarithmic_mean_temperature_difference(Thi, Tho, Tci, Tco):
    theta1 = Thi - Tco
    theta2 = Tho - Tci
    
    if theta1 == theta2:
        return theta1
    
    return (theta2-theta1)/(ln(theta2/theta1))


def get_mean_temperature_difference_correction_factor(Thi, Tho, Tci, Tco):
    R = (Thi-Tho)/(Tco-Tci)
    P = (Tco-Tci)/(Thi-Tci)
        
    if R != 1:
        return (R**2+1)**0.5 * ln((1-P)/(1-R*P)) / ((R-1) * ln((2 - P * (R+1-(R**2+1)**0.5))/(2 - P * (R+1+(R**2+1)**0.5))))

    return 2**0.5 * P / ((1-P) * ln((2 - P*(2-2**0.5))/(2 - P*(2+2**0.5))))


def get_reynolds(v, rho, D, mu):
    return D * v * rho / mu


def get_prandtl(Cp, mu, k):
    return Cp * mu / k


def get_dirty_global_heat_transfer_coefficient(Q, A, deltaTLM, correction_factor):
    return Q / (A * deltaTLM * correction_factor)


def get_fitted_total_fouling_factor(t, Rftinf, Rft0, S):
    return Rftinf - (Rftinf - Rft0) * exp(-S * t)


def get_total_number_of_tubes(Ds, Ltp, Npt, layout):
    Fc_options = {
        "quadrado": 1,
        "triangular": 0.866
    }
    
    Fs_options = {
        1: 0.93,
        2: 0.90
    }
    
    Fc = Fc_options[layout]
    Fs = Fs_options.get(Npt, Fs_options[2])
    
    return pi * Ds**2 * Fs / (4  * Ltp**2 * Fc)

if __name__ == "__main__":
    from calculations import calculate
    from utils import get_data_to_fit
    from plot_data import plot_data
    from inputs import *

    data = get_data_to_fit(filename = nome_arquivo_de_dados, dirpath = diretorio_arquivo_de_dados)

    results = calculate(
        Ds = diametro_casco,
        Dte = diametro_externo_tubo,
        Npt = numero_passes_tubos,
        Rtp = razao_passo_tubos,
        layout = arranjo_do_feixe,
        Ltube = comprimento_tubos,
        esp = espessura_tubos,
        ktube = condutividade_termica_tubos,
        Nc = numero_chicanas,
        fluid_inside_tubes = fluido_lado_dos_tubos,
        data_to_fit = data,
        tempo_final = tempo_final
    )

    plot_data(*results, dirpath = diretorio_de_saida, filename = nome_imagem_de_saida)