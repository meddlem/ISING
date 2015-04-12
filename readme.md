# 2D Ising model simulation
Monte Carlo simulation of the 2D Ising model. Two different cluster algorithms may be used: Swendsen-Wang and Wolff. This program relies on gnuplot for visualization of the system.  

##Use flags:  

1. Use -S or --SW for Swendsen Wang algorithm
2. Use -W or --Wolff for Wolff algorithm (default)
3. Use -c to calclate the correlation function
4. Use -a or --Auto to scan temperature range around critical point

Example: $./main --Wolff -a
