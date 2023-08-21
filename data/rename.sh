#!/bin/bash

# Function to rename specific files
rename_files() {
  # Check if the file exists and is a regular file
  if [[ -f "$1" ]]; then
    # Check if the file name matches the pattern "brownianLDeltahgammaD"
    if [[ $1 == *brownianLDeltahgammaD* ]]; then
      # Extract the prefix and suffix of the file name
      prefix="${1%brownianLDeltahgammaD*}"
      suffix="${1#*brownianLDeltahgammaD}"

      # Construct the new file name
      new_name="${prefix}brownianL100Delta1h0.5gamma0D1500${suffix}"

      # Rename the file
      mv "$1" "$new_name"
      echo "Renamed: $1 -> $new_name"
    fi
  fi
}

# Call the function for each file
rename_files "ee_brownianL100Delta1h0.5gamma0D1500"
rename_files "S2symbrownianLDeltahgammaD.json"
rename_files "PMbrownianLDeltahgammaD.json"
rename_files "MPbrownianLDeltahgammaD.json"

