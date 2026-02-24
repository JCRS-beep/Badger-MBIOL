### Density tests:

dens_posterior<- c(75.55, 79.89, 113.06, 114.24, 133.25, 143.50, 132.65, 143.44, 
                   154.41, 163.43, 167.22, 156.11, 151.17, 115.78, 149.87, 169.58,
                   142.61, 129.66, 112.25, 91.65, 89.41, 80.12, 75.82)
dens_posterior_mean<- mean(dens_posterior)
dens_posterior_sd<- sd(dens_posterior)

beta_reported<- (-0.239)
beta_rawdens<- beta_reported*(1-dens_posterior_mean)/dens_posterior_sd

# So now, log10(recruitment) = beta_rawdens * population_size

