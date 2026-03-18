#!/bin/bash

CONFIG_FILE="$HOME/.bashrc"
ISCE_HOME="$CONDA_PREFIX/lib/python3.11/site-packages/isce"
PROJ_PATH="$CONDA_PREFIX/share/proj"

# Set ISCE command line variables
if ! grep -q "export ISCE_HOME=" "$CONFIG_FILE"; then
  echo "export ISCE_HOME=$ISCE_HOME" >> "$CONFIG_FILE"
  echo "export PATH=$ISCE_HOME/bin:$ISCE_HOME/applications:$PATH" >> "$CONFIG_FILE"
  echo "export PYTHONPATH=$ISCE_HOME:$PYTHONPATH" >> "$CONFIG_FILE"

  echo "Variable ISCE_HOME added to $CONFIG_FILE"

else
  echo "Variable ISCE_HOME already exists in $CONFIG_FILE"
fi

# Set pyproj data and lib database variables
if ! grep -q "export PROJ_DATA=" "$CONFIG_FILE"; then
  echo "export PROJ_DATA=$PROJ_PATH" >> "$CONFIG_FILE"
  
  echo "Variable PROJ_DATA added to $CONFIG_FILE"

else
  echo "Variable PROJ_DATA already exists in $CONFIG_FILE"
fi

if ! grep -q "export PROJ_LIB=" "$CONFIG_FILE"; then
  echo "export PROJ_LIB=$PROJ_PATH" >> "$CONFIG_FILE"
  
  echo "Variable PROJ_LIB added to $CONFIG_FILE"

else
  echo "Variable PROJ_LIB already exists in $CONFIG_FILE"
fi

# Create ISCE2 kernel with proper environment variables
python -m ipykernel install --user --name isce2 --display-name "ISCE2" \
    --env ISCE_HOME "$ISCE_HOME" \
    --env PATH "$ISCE_HOME/bin:$ISCE_HOME/applications:$PATH" \
    --env PYTHONPATH "$ISCE_HOME:$PYTHONPATH" \
    --env PROJ_DATA  "$PROJ_HOME" \
    --env PROJ_LIB   "$PROJ_HOME"