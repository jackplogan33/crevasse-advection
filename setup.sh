#!/bin/bash

ACTIVATE_SCRIPT="$CONDA_PREFIX/etc/conda/activate.d/isce2-activate.sh"
DEACTIVATE_SCRIPT="$CONDA_PREFIX/etc/conda/deactivate.d/isce2-deactivate.sh"
PROJ_HOME="$CONDA_PREFIX/share/proj"

# --- Patch activation script ---
if [[ ! -f "$ACTIVATE_SCRIPT" ]]; then
  echo "ERROR: Activation script not found at $ACTIVATE_SCRIPT"
  exit 1
fi

if grep -q "_CONDA_SET_PATH" "$ACTIVATE_SCRIPT"; then
  echo "PATH backup already present in $ACTIVATE_SCRIPT, skipping."
else
  # Insert PATH backup block before the ISCE_HOME export line
  sed -i 's|^export ISCE_HOME=|if [[ -n "$PATH" ]]; then\n    export _CONDA_SET_PATH=$PATH\nfi\n\nexport ISCE_HOME=|' "$ACTIVATE_SCRIPT"
  echo "PATH backup added to $ACTIVATE_SCRIPT"
fi

if grep -q "ISCE_HOME/bin" "$ACTIVATE_SCRIPT"; then
  echo "PATH export already present in $ACTIVATE_SCRIPT, skipping."
else
  echo '' >> "$ACTIVATE_SCRIPT"
  echo 'export PATH="$ISCE_HOME/bin:$ISCE_HOME/applications:$PATH"' >> "$ACTIVATE_SCRIPT"
  echo "PATH export added to $ACTIVATE_SCRIPT"
fi

# --- Patch deactivation script ---
if [[ ! -f "$DEACTIVATE_SCRIPT" ]]; then
  echo "ERROR: Deactivation script not found at $DEACTIVATE_SCRIPT"
  exit 1
fi

if grep -q "_CONDA_SET_PATH" "$DEACTIVATE_SCRIPT"; then
  echo "PATH restore already present in $DEACTIVATE_SCRIPT, skipping."
else
  cat >> "$DEACTIVATE_SCRIPT" << 'EOF'

if [[ -n "$_CONDA_SET_PATH" ]]; then
    export PATH=$_CONDA_SET_PATH
    unset _CONDA_SET_PATH
fi
EOF
  echo "PATH restore added to $DEACTIVATE_SCRIPT"
fi

echo ""

# Create ISCE2 kernel with proper environment variables
python -m ipykernel install --user --name isce2 --display-name "ISCE2" \
    --env ISCE_HOME "$ISCE_HOME" \
    --env PATH "$ISCE_HOME/bin:$ISCE_HOME/applications:$PATH" \
    --env PYTHONPATH "$ISCE_HOME:$PYTHONPATH" \
    --env PROJ_DATA  "$PROJ_HOME" \
    --env PROJ_LIB   "$PROJ_HOME"

echo "Done. Reactivate your environment to apply changes:"
echo "  conda activate $(basename $CONDA_PREFIX)"
