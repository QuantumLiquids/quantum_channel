import os
import matplotlib.pyplot as plt
import struct

# Set the values of L, Delta, h, gamma, and D
L = 100
Delta = 1
h = 0.5
gamma = 0
D = 1001

# Construct the filename based on the variable values
filename = f"../data/ee_brownianL{L}Delta{Delta}h{h}gamma{gamma}D{D}"

# Load data from the binary file
with open(filename, "rb") as file:
    binary_data = file.read()

# Calculate the number of double values in the binary data
num_doubles = len(binary_data) // 8

# Truncate the binary data to ensure it has a length divisible by 8
binary_data = binary_data[:num_doubles * 8]

# Convert binary data to a list of double values
data = struct.unpack(f"{num_doubles}d", binary_data)

# Generate X-axis data as a sequence of integers from 1
x = list(range(1, len(data) + 1))

# Plot the entanglement entropy
plt.plot(x, data, 'o')
plt.xlabel(r'Site', fontsize=18, fontname='Times New Roman')
plt.ylabel(r'Entanglement Entropy', fontsize=18, fontname='Times New Roman')
plt.tick_params(labelsize=18)
plt.xticks(fontsize=18, fontname='Times New Roman')
plt.yticks(fontsize=18, fontname='Times New Roman')
plt.grid(True)
plt.show()

