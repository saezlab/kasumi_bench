# Set R repository
options(repos = list(CRAN = 'https://cloud.r-project.org/'))

# Install and update packages not present in the Docker image by default
install.packages("remotes")
remotes::install_github("jtanevski/mistyR@devel")
install.packages("extraDistr")
install.packages("ggrepel")
install.packages("cowplot")
install.packages("spdep")
install.packages("tidyverse")
install.packages("knn.covertree")
install.packages("ClusterR")