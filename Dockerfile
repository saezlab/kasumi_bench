FROM tanevski/mistyr:latest

# Requirements for different libraries
RUN apt update ; apt install -y libglpk-dev libgmp-dev libxml2-dev libudunits2-dev libgdal-dev libgeos-dev
# Install required R Modules
COPY requirements_interpretation.R /tmp/
RUN Rscript /tmp/requirements_interpretation.R