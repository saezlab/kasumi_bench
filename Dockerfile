FROM tanevski/mistyr:latest

# Requirements for different libraries
RUN apt update ; apt install -y libglpk-dev libgmp-dev libxml2-dev libudunits2-dev libgdal-dev libgeos-dev
# Install required R Modules
COPY requirements.R /tmp/
RUN Rscript /tmp/requirements.R
RUN Rscript -e "install.packages('tidyverse')"