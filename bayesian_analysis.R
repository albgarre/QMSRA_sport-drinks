
library(tidyverse)
library(readxl)
library(brms)
library(MASS)
library(tidybayes)
library(ggridges)
library(extraDistr)
library(RColorBrewer)

## Load data

d <- excel_sheets("./data/Datos MB- 6 años (Global).xlsx") %>%
    map(.,
        ~ read_excel("./data/Datos MB- 6 años (Global).xlsx",
                     sheet = ., skip = 1)
    ) %>%
    imap_dfr(., ~ mutate(.x, year = as.numeric(.y))) %>%
    set_names(c("season", "area", "N", "logN", "N_2", "logN_2", "year")) %>%
    mutate(year = as.factor(year))

d %>%
    ggplot() +
    geom_density_ridges(aes(x = N, y = season, fill = area), alpha = .5) +
    facet_wrap("year") 

## Multilevel model for the mould concentration in the air

# get_prior(
#     family = negbinomial,
#     formula = N ~ 0 + year + (1|area/season),
#     data = d
# )

year_model <- brm(data = d, family = negbinomial,
                   N ~ 0 + year + (1|area/season),
                   iter = 10000, warmup = 1000, cores = 2, chains = 2, 
                  seed = 1211
)

year_model
# plot(year_model)

## Supp. Table 1

year_model$fit

## Extract draws 

draws_eps2 <- as_draws_df(year_model) %>%
    dplyr::select(starts_with("r_area:season"), ".draw") %>%
    pivot_longer(-.draw, names_to = "area_season", values_to = "eps_2") %>%
    mutate(area_season = gsub("r_area:season\\[", "", area_season)) %>%
    separate(area_season, into = c("area_season"), sep = ",") %>%
    separate(area_season, into = c("area", "season"), sep = "_")

draws_eps1 <- as_draws_df(year_model) %>% 
    dplyr::select(starts_with("r_area["), ".draw") %>%
    pivot_longer(-.draw, names_to = "area", values_to = "eps_1") %>%
    mutate(area = gsub("r_area\\[", "", area)) %>%
    separate(area, into = c("area"), sep = ",")

draws_pop <- as_draws_df(year_model) %>% 
    dplyr::select(".draw", starts_with("b_year"), "shape") %>%
    pivot_longer(-c(.draw, shape), names_to = "year", values_to = "b") %>%
    mutate(year = gsub("b_year", "", year))


all_draws <- full_join(draws_pop, draws_eps1,
          by = ".draw") %>%
    full_join(draws_eps2, by = c(".draw", "area"))

## Posterior of the parameter estimates

as_draws_df(year_model) %>%
    dplyr::select(starts_with("r_area:season"), starts_with("r_area["), starts_with("b_year"), "shape") %>%
    pivot_longer(everything()) %>%
    ggplot() +
    stat_gradientinterval(aes(y = name, x = value), .width = .66)

## Marginal medians - Table 2

all_draws %>%
    mutate(logmu = b + eps_1 + eps_2,
           m = exp(logmu)) %>%
    mutate(year = 2014 + as.numeric(year)) %>%
    group_by(year) %>%
    summarize(median(m),
              quantile(m, .1),
              quantile(m, .9),
              sd(m))

all_draws %>%
    mutate(logmu = b + eps_1 + eps_2,
           m = exp(logmu)) %>%
    mutate(year = 2014 + as.numeric(year)) %>%
    group_by(area) %>%
    summarize(median(m),
              quantile(m, .1),
              quantile(m, .9),
              sd(m))

all_draws %>%
    mutate(logmu = b + eps_1 + eps_2,
           m = exp(logmu)) %>%
    mutate(year = 2014 + as.numeric(year)) %>%
    filter(area == "Tapado") %>%
    group_by(season) %>%
    summarize(median(m),
              quantile(m, .1),
              quantile(m, .9),
              sd(m))

## Posterior of mu of gamma - Figure 2

p <- all_draws %>%
    mutate(logmu = b + eps_1 + eps_2,
           m = exp(logmu)) %>%
    mutate(year = 2014 + as.numeric(year)) %>%
    full_join(.,
              tibble(season = c("Invierno", "Verano", "Otoño", "Primavera"),
                     season_plot = c("Winter", "Summer", "Autumn", "Spring"))
              ) %>%
    mutate(area = ifelse(area == "Enjuagado.de.botellas", "Rinsing",
                         ifelse(area == "Llenado", "Filling", "Capping")
                         )) %>%
    mutate(season_plot = factor(season_plot, levels = c("Spring", "Summer", "Autumn", "Winter"))) %>%
    ggplot() +
    geom_density_ridges(aes(x = m, y = season_plot, fill = area),
                        alpha = .5) +
    facet_wrap("year") +
    xlab("Fungal load (CFU/100 L)") + ylab("") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          legend.title = element_blank())

ggsave(plot = p,
       "Figure_2.png",
       width = 9, height = 6, units = "in"
       )

## Posterior gamma density

all_draws %>%
    mutate(logmu = b + eps_1 + eps_2,
           lambda = exp(logmu)) %>%
    sample_n(1e4) %>%
    mutate(i = row_number()) %>%
    split(.$i) %>%
    map(.,
        ~ tibble(
            year = .$year,
            area = .$area,
            season = .$season,
            x = seq(.1, 10, lenght = 100),
            y = dgamma(x, shape = .$shape,  scale = .$lambda/.$shape)
            )
        ) %>%
    imap_dfr(., ~ mutate(.x, iter = .y)) %>%
    group_by(year, area, season, x) %>%
    summarize(
        m_y = median(y, na.rm = TRUE),
        q10 = quantile(y, .1, na.rm = TRUE),
        q90 = quantile(y, .9, na.rm = TRUE)
    ) %>%
    ggplot(aes(x = x, fill = season, colour = season)) +
    geom_ribbon(aes(ymin = q10, ymax = q90), alpha = .1) +
    facet_grid(area ~ year) +
    scale_x_log10()

## Get the spoilage data

d_spoil <- read_excel("./data/spore_spoilage.xlsx") %>%
    mutate(temperature = as.factor(temperature),
           # N0 = as.factor(N0),
           ensayo = as.factor(ensayo),
           total = 60) # %>%
    # filter(N0 == 1)

d_spoil %>%
    mutate(percentage = ifelse(N0 == 10, percentage^1, percentage)) %>%
    ggplot(aes(x = hour, y = percentage, colour = factor(N0))) +
    geom_point() +
    geom_line(aes(linetype = factor(temperature)))

## Spoilage model

spoil_model <- brm(positives | trials(total) ~ 1 + hour + hour:temperature,  
                   data = d_spoil %>% filter(N0 == 1), 
                   family = binomial(link = "logit"),
                   warmup = 500, 
                   iter = 4000, 
                   chains = 2,  
                   cores = 2,
                   seed = 2124)
spoil_model
spoil_model$fit
plot(spoil_model)

## Supp. table 2

spoil_model$fit

## Figure 3

p <- as_draws_df(spoil_model) %>% 
    mutate(hour = list(seq(0, 200, length = 100))) %>%
    unnest(hour) %>%
    mutate(
        intercept = b_Intercept,
        slope_20 = b_hour,
        slope_30 = b_hour + `b_hour:temperature30`,
        mu_20 = intercept + slope_20*hour,
        mu_30 = intercept + slope_30*hour,
        pred_20 = exp(mu_20)/(1 + exp(mu_20)),
        pred_30 = exp(mu_30)/(1 + exp(mu_30))
    ) %>%
    dplyr::select(hour, starts_with("pred")) %>%
    pivot_longer(-hour, values_to = "pred") %>%
    group_by(hour, name) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low = quantile(pred, prob = 0.05/2, na.rm = TRUE),
              pred_high = quantile(pred, prob = 1 - 0.05/2, na.rm = TRUE)) %>%
    separate(name, into = c("p", "temperature")) %>%
    mutate(temperature = paste0(temperature, "ºC")) %>%
    mutate(temperature = as.factor(temperature)) %>%
    ggplot(aes(x = hour, y = pred_m, colour = temperature,
               fill = temperature)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_low, ymax = pred_high), alpha = .5) +
    # facet_wrap("temperature") +
    geom_point(aes(x = hour, y = percentage),
               data = d_spoil %>% 
                   filter(N0 == 1) %>%
                   mutate(temperature = paste0(as.character(temperature), "ºC"))
                   ) +
    xlab("Storage time (h)") + 
    ylab("Probability of visual spoilage from 1 spore") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          legend.title = element_blank())

ggsave(plot = p,
       "Figure_3.png",
       width = 9, height = 6, units = "in"
)

## Values

aa <- as_draws_df(spoil_model) %>% 
    mutate(hour = list(seq(0, 200, length = 100))) %>%
    unnest(hour) %>%
    mutate(
        intercept = b_Intercept,
        slope_20 = b_hour,
        slope_30 = b_hour + `b_hour:temperature30`,
        mu_20 = intercept + slope_20*hour,
        mu_30 = intercept + slope_30*hour,
        pred_20 = exp(mu_20)/(1 + exp(mu_20)),
        pred_30 = exp(mu_30)/(1 + exp(mu_30))
    ) %>%
    dplyr::select(hour, starts_with("pred")) %>%
    pivot_longer(-hour, values_to = "pred") %>%
    group_by(hour, name) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low = quantile(pred, prob = 0.05/2, na.rm = TRUE),
              pred_high = quantile(pred, prob = 1 - 0.05/2, na.rm = TRUE)) %>%
    separate(name, into = c("p", "temperature"))
    
aa %>% 
    split(.$temperature) %>%
    map(.,
        ~ approx(x = .$hour, y = .$pred_m, xout = 100)
        )

## Get the draws

spoil_draws <- as_draws_df(spoil_model)  %>%
    dplyr::select(.draw, starts_with("b_")) %>%
    pivot_longer(-c(.draw, b_Intercept), values_to = "b_T",
                 names_to = "temperature") %>%
    mutate(temperature = ifelse(temperature == "b_hour", 20, 30)) %>%
    mutate(temperature = as.factor(temperature))
    
## Spoilage Monte Carlo

set.seed(1244)

sims_bottle <- all_draws %>%
    filter(year == 1, season == "Verano", area == "Tapado") %>%
    sample_n(5000) %>%
    mutate(logmu = b + eps_1 + eps_2,
           mu = exp(logmu)*.1  # convert from /100 L to /m3
           ) %>%
    mutate(
        Cair = list(rgamma(50, shape = shape, scale = mu/shape)),
    ) %>%
    unnest(Cair) %>%
    mutate(.,
           d_spore = rtriang(nrow(.), 2, 30, 2.6),  # um
           vs_cm = (d_spore/18.02)^2,  # cm/s
           vs = vs_cm/100,  #m/s
           d_bottle = 38*1e-3,  # m
           Area = pi*(d_bottle/2)^2,  # m2
           t_exp = rtriang(nrow(.), 3, 10, 4)  # s
    ) %>%
    mutate(lambda_bottle = Cair*vs*Area*t_exp) 
        
# View(sims_bottle)

sims_bottle %>%
    summarize(mean = mean(lambda_bottle), 
              median = median(lambda_bottle),
              q90 = quantile(lambda_bottle, .9),
              q10 = quantile(lambda_bottle, .1))

## Figure 4

p1 <- ggplot(sims_bottle) + 
    geom_density(aes(lambda_bottle), fill = "black", alpha = .5) +
    geom_vline(aes(xintercept = x), linetype = 2,
               data = 
                   sims_bottle %>% 
                   summarize(x = mean(lambda_bottle))
    ) +
    geom_vline(aes(xintercept = x), linetype = 3,
               data = 
                   sims_bottle %>% 
                   summarize(x = median(lambda_bottle))
    ) +
    xlab("Expected mould contamination per bottle (CFU/bottle)") +
    ylab("Probability density") +
    theme_bw(base_size = 14)

p2 <- ggplot(sims_bottle) + 
    geom_density(aes(lambda_bottle), fill = "black", alpha = .5) +
    geom_vline(aes(xintercept = x), linetype = 2,
               data = 
                   sims_bottle %>% 
                   summarize(x = mean(lambda_bottle))
    ) +
    geom_vline(aes(xintercept = x), linetype = 3,
               data = 
                   sims_bottle %>% 
                   summarize(x = median(lambda_bottle))
    ) +
    scale_x_log10() +
    xlab("Expected mould contamination per bottle (CFU/bottle)") +
    ylab("Probability density") +
    theme_bw(base_size = 14)

cowplot::plot_grid(p1, p2, labels = "AUTO")

ggsave(
       "Figure_4.png",
       width = 12, height = 4, units = "in"
)

##############

aa <- all_draws %>%
    filter(year == 1, area == "Tapado") %>%
    sample_n(5000) %>%
    mutate(logmu = b + eps_1 + eps_2,
           mu = exp(logmu)*.1  # convert from /100 L to /m3
    ) %>%
    mutate(
        Cair = list(rgamma(50, shape = shape, scale = mu/shape)),
    ) %>%
    unnest(Cair) %>%
    mutate(.,
           d_spore = rtriang(nrow(.), 2, 30, 2.6),  # um
           vs_cm = (d_spore/18.02)^2,  # cm/s
           vs = vs_cm/100,  #m/s
           d_bottle = 38*1e-3,  # m
           Area = pi*(d_bottle/2)^2,  # m2
           t_exp = rtriang(nrow(.), 3, 10, 4)  # s
    ) %>%
    mutate(lambda_bottle = Cair*vs*Area*t_exp) %>%
    mutate(
        Cbottle = rpois(nrow(.), lambda = lambda_bottle)
    ) 

aa %>%
    filter(Cbottle > 0) %>%
    # pull(Cbottle)
    ggplot() + 
    geom_histogram(aes(Cbottle))

aa %>% filter(Cbottle > 0)
aa$Cbottle %>% mean()

## Spoilage Monte Carlo - seasons

set.seed(1244)

sims_bottle <- all_draws %>%
    filter(area == "Tapado") %>%
    sample_n(50000) %>%
    mutate(logmu = b + eps_1 + eps_2,
           mu = exp(logmu)*.1  # convert from /100 L to /m3
    ) %>%
    mutate(
        Cair = list(rgamma(50, shape = shape, scale = mu/shape)),
    ) %>%
    unnest(Cair) %>%
    mutate(.,
           d_spore = rtriang(nrow(.), 2, 30, 2.6),  # um
           vs_cm = (d_spore/18.02)^2,  # cm/s
           vs = vs_cm/100,  #m/s
           d_bottle = 38*1e-3,  # m
           Area = pi*(d_bottle/2)^2,  # m2
           t_exp = rtriang(nrow(.), 3, 10, 4)  # s
    ) %>%
    mutate(lambda_bottle = Cair*vs*Area*t_exp) 

## Sup. Table 4

sims_bottle %>%
    group_by(year, season) %>%
    summarize(
              median = median(lambda_bottle),
              q90 = quantile(lambda_bottle, .9),
              q10 = quantile(lambda_bottle, .1)) %>%
    write_excel_csv("supp_table4.csv")

ggplot(sims_bottle) + 
    geom_density(aes(lambda_bottle, fill = season), alpha = .5) +
    scale_x_log10() +
    xlab("Expected mould contamination per bottle (CFU/bottle)") +
    ylab("Probability density")

## Medians

all_draws %>%
    filter(year == 1, area == "Tapado") %>%
    sample_n(5000) %>%
    mutate(logmu = b + eps_1 + eps_2,
           mu = exp(logmu)*.1  # convert from /100 L to /m3
    ) %>%
    mutate(
        Cair = list(rgamma(50, shape = shape, scale = mu/shape)),
    ) %>%
    unnest(Cair) %>%
    mutate(.,
           d_spore = rtriang(nrow(.), 2, 30, 2.6),  # um
           vs_cm = (d_spore/18.02)^2,  # cm/s
           vs = vs_cm/100,  #m/s
           d_bottle = 38*1e-3,  # m
           Area = pi*(d_bottle/2)^2,  # m2
           t_exp = rtriang(nrow(.), 3, 10, 4)  # s
    ) %>%
    mutate(lambda_bottle = Cair*vs*Area*t_exp) %>%
    dplyr::select(Cair, vs, Area, t_exp, lambda_bottle) %>%
    pivot_longer(everything()) %>%
    group_by(name) %>%
    summarize(median(value))

## Orders of magnitude

(10/18.02)^2/100
pi*(38*1e-3/2)^2

## Figure 5

my_cols <- brewer.pal(5, "RdYlGn")

tibble(
    O_air = -3:1,
    O_t_m6 = -6 + 6 - O_air,
    O_t_m12 = -12 + 6 - O_air,
    O_t_m9 = -9 + 6 - O_air
) %>%
    pivot_longer(-O_air) %>%
    ggplot(aes(x = O_air, y = value, colour = name)) +
    geom_point(size = 4) +
    geom_line(linetype = 1, size = 1) +
    geom_label(aes(x = 1, y = 4, label = "NON-COMPLIANT"),
               inherit.aes = FALSE,
               size = 8,
               # fill = "darkorange3",
               fill = my_cols[1], hjust = 1,
               colour = "white") +
    geom_label(aes(x = -3, y = -6, label = "COMPLIANT"),
               inherit.aes = FALSE,
               size = 8, hjust = 0,
               # fill = "chartreuse4",
               fill = my_cols[5],
               colour = "white") +
    geom_point(aes(x = -1, y = 0),
               size = 5, shape = 9, colour = "black",
               inherit.aes = FALSE) +
    geom_text(aes(x = 1, y = -1.5, label = "ALOS: -6"),
              hjust = 1,
              colour = "tan1"
    ) +
    geom_text(aes(x = 1, y = -4.5, label = "ALOS: -9"),
              hjust = 1,
              colour = "tan2"
    ) +
    geom_text(aes(x = 1, y = -7.5, label = "ALOS: -12"),
              hjust = 1,
              colour = "tan3"
    ) +
    theme_bw(base_size = 14) +
    scale_y_continuous(breaks = seq(-6, 4, by = 2),
                       name = "Order of magnitude of the exposure time (s)") +
    scale_x_continuous(breaks = seq(-3, 1, by = 1),
                       name = bquote(Order~of~magnitude~of~the~air~mould~concentration~(CFU/m^3))
                       ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = -1, linetype = 3) +
    scale_colour_manual(values = c("tan1", "tan2", "tan3")) +
    theme(legend.position = "none",
          axis.ticks = element_line(size = 0),
          panel.grid.minor = element_line(size = 0))

ggsave(
    "Figure_5.png",
    width = 9, height = 6, units = "in"
)




