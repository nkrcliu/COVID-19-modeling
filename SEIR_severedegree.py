import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import datetime
import time

model_parameters = {
    'model': 'seir',
    'use_harvard_params': False,  # If True use the harvard parameters directly, if not calculate off the above
    'fix_r0': False,  # If True use the parameters that make R0 2.4, if not calculate off the above
    'days_to_model': 270,

    "population": 6949503,
    "hospitalization_rate": 0.073,
    "hospitalized_cases_requiring_icu_care": 0.3,  # ?????

    "duration_mild_infections": 6,
    "hospital_time_recovery": 11,  # Duration of hospitalization, In days

    'total_infected_period': 12,  # In days
    'icu_time_death': 7,  # Time from ICU admission to death, In days
    'case_fatality_rate': .0109341104294479,

    'beta': 0.5,
    'beta_hospitalized': 0.1,
    'beta_icu': 0.1,

    'presymptomatic_period': 5,
    'exposed_from_infected': True,

    'hospital_capacity_change_daily_rate': 1.05,
    'max_hospital_capacity_factor': 2.07,
    'initial_hospital_bed_utilization': .6,
    'exposed_infected_ratio': 1
}
pop_dict={
    "total":6949503,
    "infected":164,
    "recovered":0,
    "deaths":1
}

def deriv(y0, t, beta, alpha, gamma, rho, mu, N):
    dy = [0, 0, 0, 0, 0, 0]
    S = np.max([N - sum(y0), 0])

    dy[0] = np.min([(np.dot(beta[1:4], y0[1:4]) * S), S]) - (alpha * y0[0])  # Exposed
    dy[1] = (alpha * y0[0]) - (gamma[1] + rho[1]) * y0[1]  # Ia - Mildly ill
    dy[2] = (rho[1] * y0[1]) - (gamma[2] + rho[2]) * y0[2]  # Ib - Hospitalized
    dy[3] = (rho[2] * y0[2]) - ((gamma[3] + mu) * y0[3])  # Ic - ICU
    dy[4] = np.min([np.dot(gamma[1:4], y0[1:4]), sum(y0[1:4])])  # Recovered
    dy[5] = mu * y0[3]  # Deaths

    return dy

# for now just implement Harvard model, in the future use this to change
# key params due to interventions                             改变参数
def generate_epi_params(model_parameters):
    N = model_parameters["population"]

    fraction_critical = (
        model_parameters["hospitalization_rate"]
        * model_parameters["hospitalized_cases_requiring_icu_care"]
    )

    fraction_severe = model_parameters["hospitalization_rate"] - fraction_critical

    alpha = 1 / model_parameters["presymptomatic_period"]

    # assume hospitalized don't infect
    beta = [
        0,
        model_parameters["beta"] / N,
        model_parameters["beta_hospitalized"] / N,
        model_parameters["beta_icu"] / N,
    ]

    # have to calculate these in order and then put them into arrays
    gamma_0 = 0
    gamma_1 = (1 / model_parameters["duration_mild_infections"]) * (
        1 - model_parameters["hospitalization_rate"]
    )
    rho_0 = 0
    rho_1 = (1 / model_parameters["duration_mild_infections"]) - gamma_1

    rho_2 = (1 / model_parameters["hospital_time_recovery"]) * (
        (fraction_critical / (fraction_severe + fraction_critical))
    )

    gamma_2 = (1 / model_parameters["hospital_time_recovery"]) - rho_2

    mu = (1 / model_parameters["icu_time_death"]) * (
        model_parameters["case_fatality_rate"] / fraction_critical
    )

    gamma_3 = (1 / model_parameters["icu_time_death"]) - mu

    seir_params = {
        "beta": beta,
        "alpha": alpha,
        "gamma": [gamma_0, gamma_1, gamma_2, gamma_3],
        "rho": [rho_0, rho_1, rho_2],
        "mu": mu,
    }

    return seir_params

seir_params=generate_epi_params(model_parameters)
a = seir_params["alpha"]
b = seir_params["beta"]
p = seir_params["rho"]
g = seir_params["gamma"]
u = seir_params["mu"]

r0_dif_interventions = np.array([[3.7,3.7,3.7,3.7,3.7,3.7,3.7,3.7,3.7,3.7,3.7,3.7],[1.3,0.3,0.3,0.3,0.3,0.3,0.3,0.2,0.035,0.035,0.035,0.035],[1.1,1.1,1.1,1.1,1,1,1,1,0.8,0.8,0.8,0.8],[1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7]])
r0_dif_interventions


def inverse_calc_beta(seir_params, N):
    BETA = np.zeros(shape=(4, 12, 4))
    BETA[:,:,2] = seir_params["beta"][2]
    BETA[:,:,3] = seir_params["beta"][3]
    for i in range(0, 4):
        for j in range(0, 12):
            BETA[i][j][1] = (r0_dif_interventions[i][j] - (p[1] / (p[1] + g[1])) * (b[2] / (p[2] + g[2])) - (
                        p[2] / (p[2] + g[2])) * (b[3] / (u + g[3]))) / N * (p[1] + g[1])
    print(BETA[1,:,1])
    return BETA

def seir(
    pop_dict, model_parameters, BETA, alpha, gamma, rho, mu, intervention
):

    N = pop_dict["total"]
    # assume that the first time you see an infected population it is mildly so
    # after that, we'll have them broken out
    if "infected_b" in pop_dict:
        mild = pop_dict["infected_a"]
        hospitalized = pop_dict["infected_b"]
        icu = pop_dict["infected_c"]
    else:
        hospitalized = pop_dict["infected"] / 4
        mild = hospitalized / model_parameters["hospitalization_rate"]
        icu = hospitalized * model_parameters["hospitalized_cases_requiring_icu_care"]

    exposed = model_parameters["exposed_infected_ratio"] * mild

    susceptible = pop_dict["total"] - (
        pop_dict["infected"] + pop_dict["recovered"] + pop_dict["deaths"])

    y0 = [
        int(exposed),
        int(mild),
        int(hospitalized),
        int(icu),
        int(pop_dict.get("recovered", 0)),
        int(pop_dict.get("deaths", 0)),
    ]

    out = np.zeros((7, 6, 12))

    for i in range(0, 12):

        steps = 7

        t = np.arange(0, steps, 1)

        out[:, :, i] = odeint(deriv, y0, t, args=(BETA[intervention, i, :], alpha, gamma, rho, mu, N))

        y0 = out[-1, :, i]

    ret = out[:, :, 0]

    for i in range(1, 12):
        ret = np.vstack((ret, out[:, :, i]))

    ret = np.hstack((N - np.sum(ret, axis=1, keepdims=True), ret))
    print(ret)


    return ret

ret=seir(pop_dict, model_parameters, inverse_calc_beta(seir_params, model_parameters['population']), a, g, p, u, 2)
ret1=np.transpose(ret)
tvec=np.arange(0,84,1)
plt.figure(figsize=(13,5))
plt.plot(tvec,ret)
plt.xlabel("Time (days)")
plt.ylabel("Cases in  M")
plt.legend(("S","E","I1","I2","I3","R","D"))
plt.show()
