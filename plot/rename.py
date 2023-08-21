import os

# Specify the directory path
directory = '../data/'

# Get the list of files in the directory
files = os.listdir(directory)

# Loop through each file
for filename in files:
    # Check if the file is a JSON file
    if filename.endswith('.json') or filename.startswith('ee'):
        # Extract the gamma value from the filename
        gamma_str = filename.split('gamma')[1].split('D')[0]

        if len(gamma_str) > 10:
            gamma = float(gamma_str)

            for num_decimals in range(3):
                gamma_approx = round(gamma, num_decimals)
                if abs(gamma - gamma_approx) < 1e-13:
                    # If the difference is within tolerance, keep the corresponding decimal places
                    modified_gamma = format(gamma_approx, f'.{num_decimals}f')
                    break
            else:
                # If no approximation is within tolerance, use the original gamma with one decimal place
                modified_gamma = format(gamma, '.5f') if gamma % 1 else str(int(gamma))


	    # Construct the new filename with modified gamma value
            new_filename = filename.replace('gamma' + gamma_str, 'gamma' + modified_gamma)

            # Rename the file with the new filename
            os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))
