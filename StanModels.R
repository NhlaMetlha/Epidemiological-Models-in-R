setwd("C:/Users/Lebo/Desktop/thesis")

#load required packages
library(rstan)
library(ggplot2)
library(gridExtra)
library(StanHeaders)
library(rstudioapi)
library(tidyverse)
library(plotly)
library(dplyr)
library(lubridate)
library(bayesplot)
library(paramtest)
library(moments)
library(tidybayes)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#-----Importing the data-----

library(readxl)
Data = read_excel("DataRef1.xlsx", sheet = "Sheet1")
attach(Data)

#-----Data visualization-----

fig = plot_ly(Data, type = 'scatter', mode = 'lines') %>%
  add_trace(x = ~ date, y = ~ new_cases, name = 'Number of reported cases') %>%
  layout(showlegend = F)

fig = fig %>%
  layout(
    xaxis = list(zerolinewidth = 1,
                 title = "Date"),
    yaxis = list(zerolinewidth = 1,
                 title = "Confirmed Cases"),
    width = 900)

fig

Data %>% 
  ggplot(aes(x = week_cases) +
           geom_density(fill = "cornflowerblue", alpha = .5) +
           scale_x_continuous(labels = "week"))

#----DESCRIPTIVE STATISTICS----

summary(week_cases)
skewness(week_cases)
kurtosis(week_cases)
sd(week_cases)

#----FITTING SIR MODEL-----

# BW population estimated to be 2,346,179 (Statistics Botswana)
N = 2346179;

# times
n_days = length(week_cases) 
t = seq(1, n_days, by = 1)
t0 = 0 
t = t

# initial conditions
i0 = 1
s0 = N - i0
r0 = 0
y0 = c(S = s0, I = i0, R = r0)

# data for Stan
data_sir = list(n_days = n_days, y0 = y0, 
                t0 = t0, ts = t, N = N, 
                week_cases = week_cases)

# Next we compile the model, saved in the file `SIR.stan`

model = stan_model("SIR.stan")

# run MCMC
fit_sir = sampling(model,
                   data = data_sir,
                   iter = 2000,   # number of MCMC steps
                   warmup = 1000,
                   chains = 4,
                   seed = 0)

loo_SIR = loo(fit_sir)
loo_SIR

y_rep = as.matrix(fit_sir, par = "pred_cases")
dim(y_rep)

ppc_hist(week_cases, y_rep[1:8, ], binwidth = 1)
pred_cases_1 = as.matrix(fit_sir, pars = "pred_cases")
ppc_dens_overlay(week_cases, pred_cases_1[1:10, ])

### Checking the inference
check_hmc_diagnostics(fit_sir)

pars = c("lambda", "mu", "R0", "omega", "phi_inv")
print(fit_sir, pars = pars)


bayesplot_theme_set()
color_scheme_set("brewer-Spectral")
mcmc_trace(fit_sir, pars = c("lambda", "mu", "R0", "omega"))



#----FITTING SEIR MODEL----

# times
N = 2346179;
n_days = length(week_cases) 
t = seq(0, n_days, by = 1)
t0 = 0 
t = t[-1]


data_seir = list(
  cases = week_cases,
  n_days = length(week_cases),
  t0 = 0,
  N = N)
data_seir$ts = seq(1, data_seir$n_days, by = 1)


# Next we compile the model, saved in the file `StanSEIR.stan`

model_seir = stan_model("StanSEIR.stan")

# run MCMC

fit_seir = sampling(model_seir,
                   data = data_seir,
                   iter = 2000,   # number of MCMC steps
                   warmup = 1000,
                   chains = 4,
                   seed = 0)

loo_SEIR = loo(fit_seir)
loo_SEIR


loo_compare(loo_SIR, loo_SEIR)
library(loo)

help('pareto-k-diagnostic')

#checking the inference

check_hmc_diagnostics(fit_seir)

pars = c("lambda","sigma", "mu", "phi_inv", "R0", "omega")
print(fit_seir, pars = pars)

bayesplot_theme_set()
color_scheme_set("brewer-Spectral")
mcmc_dens_overlay(fit_seir, pars = c("lambda","sigma", "mu", "R0"))



library(ggdist)
library(magrittr)
library(ggplot2)
library(data.table)
library(deSolve)


extract(fit_sir, pars = "pred_cases")[[1]] %>%
  as.data.frame() %>%
  mutate(.draw = 1:n()) %>%
  tidyr::gather(key,value, -.draw) %>%
  mutate(step = readr::parse_number(stringr::str_extract(key,"\\d+"))) %>%
  group_by(step) %>%
  curve_interval(value, .width = c(.5, .8, .95)) %>%
  ggplot(aes(x = step, y = value)) +
  geom_hline(yintercept = 1, color = "gray75", linetype = "dashed") +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  scale_fill_brewer() +
  labs(
    title = "Simulated SIR Curve for Infections",
    y = "Cases"
  )+
  geom_point(data = tibble(week_cases = week_cases, t = 1:length(week_cases)),
             aes(t, week_cases), inherit.aes = FALSE, colour = "orange")+
  theme_minimal()


#----Change Point Detection----
library(onlineBcp)
library(bcp)

x = sqrt(week_cases + 3/5)

bcp = online_cp(x, theta = .7, alpha = 1.1, 
                beta = 1.2, th_cp = .3)


res = summary(bcp, norm.test = T)
bcp = summary(bcp)
print(res, digits = 4)

bcp.new = summary(bcp)
plot(bcp.new, ylab = "New Cases", xlab = "Weeks")
?plot
# change point plot
bcp_bot = bcp(week_cases, return.mcmc = TRUE)
plot(bcp_bot)

bcp_sum = as.data.frame(summary(bcp_bot))


# Let's filter the data frame and identify the year:
bcp_sum$id = 1:length(week_cases)
(sel = bcp_sum[which(bcp_bot$posterior.prob > 0.9), ])
# Get the year:
time(new_cases)[sel$id]


#----SIR WITH CONTROL MEASURES----

week_switch = "4" #week of introduction of control measures
tswitch = Data %>% filter(week < week_switch) %>% nrow() + 1 # convert time to number

data_forcing = list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, 
                    cases = week_cases, tswitch = tswitch)

model_forcing = stan_model("SIRchange_point.stan")

fit_forcing = sampling(model_forcing, 
                       data_forcing, 
                       iter = 2000,
                       seed = 4)

check_hmc_diagnostics(fit_forcing)

pars = c("lambda","sigma", "mu", "phi_inv", "R0",
         "recovery_time", "eta", "nu", "xi_raw")
print(fit_forcing, pars = pars)

#----SEIR WITH CONTROL MEASURES----

week_switch = "10" #week of introduction of control measures
tswitch = Data %>% filter(week < week_switch) %>% nrow() + 1 # convert time to number

data_forcing = list(n_days = n_days, t0 = t0, ts = t, N = N, 
                    cases = week_cases, tswitch = tswitch)

model_forcing = stan_model("SEIRChange_point.stan")

fit_forcing = sampling(model_forcing, 
                        data_forcing, 
                        iter = 2000,
                        seed = 4)

#checking the inference

check_hmc_diagnostics(fit_forcing)

pars = c("lambda","sigma", "mu", "phi_inv", "R0", "incubation_time", 
         "recovery_time", "eta", "nu", "xi_raw")
print(fit_forcing, pars = pars)



week_switch = "4" #week of introduction of control measures
tswitch = Data %>% filter(week < week_switch) %>% nrow() + 1 # convert time to number

data_forcing = list(n_days = n_days, t0 = t0, ts = t, N = N, 
                    cases = week_cases, tswitch = tswitch)

model_forcing = stan_model("SEIRChange_point.stan")

fit_forcing = sampling(model_forcing, 
                       data_forcing, 
                       iter = 2000,
                       seed = 4)

check_hmc_diagnostics(fit_forcing)

pars = c("lambda","sigma", "mu", "phi_inv", "R0", "incubation_time", 
         "recovery_time", "eta", "nu", "xi_raw")
print(fit_forcing, pars = pars)


week_switch_2 = "75" #week of introduction of control measures
tswitch = Data %>% filter(week < week_switch_2) %>% nrow() + 1 # convert time to number

data_forcing = list(n_days = n_days, t0 = t0, ts = t, N = N, 
                    cases = week_cases, tswitch = tswitch)

model_forcing = stan_model("SEIRChange_point.stan")

fit_forcing = sampling(model_forcing, 
                       data_forcing, 
                       iter = 2000,
                       seed = 4)

check_hmc_diagnostics(fit_forcing)
pars = c("lambda","sigma", "mu", "phi_inv", "R0", "incubation_time", 
         "recovery_time", "eta", "nu", "xi_raw", "reporting_D")
print(fit_forcing, pars = pars)



print(N)
print(week_cases)
qplot(week_cases)

