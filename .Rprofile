# Modified by rsconnect package 0.8.29 on 2023-03-16 13:28:04:
# renv initialization disabled in published application
# source("renv/activate.R")

VIRTUALENV_NAME = 'example_env_name'

if (Sys.info()[['user']] == 'shiny'){
  
  # Running on shinyapps.io
  Sys.setenv(PYTHON_PATH = 'python3')
  Sys.setenv(VIRTUALENV_NAME = VIRTUALENV_NAME) # Installs into default shiny virtualenvs dir
  Sys.setenv(RETICULATE_PYTHON = paste0('/home/shiny/.virtualenvs/', VIRTUALENV_NAME, '/bin/python'))
  
} else if (Sys.info()[['user']] == 'rstudio-connect'){
  
  # Running on remote server
  Sys.setenv(PYTHON_PATH = '/opt/python/3.7.7/bin/python3')
  Sys.setenv(VIRTUALENV_NAME = paste0(VIRTUALENV_NAME, '/')) # include '/' => installs into rstudio-connect/apps/
  Sys.setenv(RETICULATE_PYTHON = paste0(VIRTUALENV_NAME, '/bin/python'))
  
} else {
  
  # Running locally
  options(shiny.port = 7450)
  Sys.setenv(PYTHON_PATH = 'python')
  Sys.setenv(VIRTUALENV_NAME = VIRTUALENV_NAME) # exclude '/' => installs into ~/.virtualenvs/
  # RETICULATE_PYTHON is not required locally, RStudio infers it based on the ~/.virtualenvs path
}
