#!/usr/bin/env python3
import math

# Beam size
beam_size = float(input("Input beam size(mm) : "))

# Calculate phi value
tan_phi = (beam_size/10) / 40 # Length of beam pipe is 40 cm
phi = math.degrees(math.atan(tan_phi))

# Calculate theta value
theta_min = 90 - phi
theta_max = 90 + phi

# Calculate eta value
eta_min = -math.log(math.tan(math.radians(theta_min) / 2))
eta_max = -math.log(math.tan(math.radians(theta_max) / 2))

# Output
print(f"phi: {phi} degree")
print(f"theta_min：{theta_min} degree")
print(f"theta_max：{theta_max} degree")
print(f"eta_min：{eta_min}")
print(f"eta_max：{eta_max}")
