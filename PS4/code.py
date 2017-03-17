# Import libraries
import numpy as np #Load numpy library for matrices operations
import matplotlib.pyplot as plt # library for graphics
from matplotlib.backends.backend_pdf import PdfPages # save graphics as pdf
import control #Load control for solving Riccati Equation
np.random.seed(1234)

# System components that are not going to be modified
N = 100 # horizon (time periods)
n = 2

A = np.array([[2, 0], [1, 0]])
B = np.array([[0, 2], [1, 1]])
C = np.array([[0, 3]])

Q = np.matmul(np.transpose(C), C)
#np.linalg.eigvals(Q)

# Behaviour of the system
# ii --------------------------------------------------
# Fix R and x[0], D_1 = normal - D_2 = "much larger"
R = np.diag(np.repeat(1, [n], axis = 0))
x_1 = np.empty([N + 1, n])
x_1[0] = [1, 2]
D_1 = np.diag(np.repeat(1, [n], axis = 0))
w_1 = np.random.multivariate_normal(mean = np.repeat(0, [n], axis = 0), cov = D_1, size = N)
x_2 = np.empty([N + 1, n])
x_2[0] = x_1[0]
D_2 = np.diag(np.repeat(5, [n], axis = 0))
w_2 = np.random.multivariate_normal(mean = np.repeat(0, [n], axis = 0), cov = D_2, size = N)
#np.linalg.eigvals(K)
K = np.empty([N + 1, n, n])
K[N] = Q
for i in range(100,0,-1):
    K[i-1] = np.matmul(np.matmul(np.transpose(A), K[i] - np.matmul(np.matmul(K[i], B), np.linalg.inv(np.matmul(np.matmul(np.transpose(B), K[i]), B) + R), np.matmul(np.transpose(B), K[i]))), A) + Q
    
L = np.empty([N, n, n])
for i in range(0, N):
    L[i] = - np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(B), K[i+1]), B) + R), np.matmul(np.matmul(np.transpose(B), K[i+1]), A))

for i in range(0, N):
    # Values normal
    x_1[i+1] = np.matmul(A + np.matmul(B, L[i]), x_1[i]) + w_1[i]
    # Values much larger
    x_2[i+1] = np.matmul(A + np.matmul(B, L[i]), x_2[i]) + w_2[i]

# plot
plot_fig2 = plt.figure(2)
plt.subplot(211)
plt.title('Behavior of the system for two covariance matrices (D)')
plt.ylabel('first element')
plt.plot(x_1.T[0])
plt.plot(x_2.T[0])
plt.subplot(212)
plt.ylabel('second element')
plt.plot(x_1.T[1], label = 'D "normal"')
plt.plot(x_2.T[1], label = 'D "much larger"')
plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4), ncol=2)

# iii --------------------------------------------------
# Fix R and D, x_1[0] = normal - x_2[0] = "much larger"
R = np.diag(np.repeat(1, [n], axis = 0))
x_1 = np.empty([N + 1, n])
x_1[0] = [1, 2]
D = np.diag(np.repeat(1, [n], axis = 0))
w = np.random.multivariate_normal(mean = np.repeat(0, [n], axis = 0), cov = D, size = N)
x_2 = np.empty([N + 1, n])
x_2[0] = [35, 41]
K = np.empty([N + 1, n, n])
K[N] = Q
for i in range(100,0,-1):
    K[i-1] = np.matmul(np.matmul(np.transpose(A), K[i] - np.matmul(np.matmul(K[i], B), np.linalg.inv(np.matmul(np.matmul(np.transpose(B), K[i]), B) + R), np.matmul(np.transpose(B), K[i]))), A) + Q
    
L = np.empty([N, n, n])
for i in range(0, N):
    L[i] = - np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(B), K[i+1]), B) + R), np.matmul(np.matmul(np.transpose(B), K[i+1]), A))

for i in range(0, N):
    # Values normal
    x_1[i+1] = np.matmul(A + np.matmul(B, L[i]), x_1[i]) + w[i]
    # Values much larger
    x_2[i+1] = np.matmul(A + np.matmul(B, L[i]), x_2[i]) + w[i]

# plot
plot_fig3 = plt.figure(3)
plt.subplot(211)
plt.title('Behavior of the system for two initial states (x0)')
plt.ylabel('first element')
plt.plot(x_1.T[0])
plt.plot(x_2.T[0])
plt.subplot(212)
plt.ylabel('second element')
plt.plot(x_1.T[1], label = 'x0 "normal"')
plt.plot(x_2.T[1], label = 'x0 "much larger"')
plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4), ncol=2)

# iv --------------------------------------------------
# Fix x[0] and D, R_1[0] = normal - R_2[0] = "much larger"
R_1 = np.diag(np.repeat(1, [n], axis = 0))
R_2 = np.diag(np.repeat(100, [n], axis = 0))
x_1 = np.empty([N + 1, n])
x_1[0] = [1, 2]
D = np.diag(np.repeat(1, [n], axis = 0))
w = np.random.multivariate_normal(mean = np.repeat(0, [n], axis = 0), cov = D, size = N)
x_2 = np.empty([N + 1, n])
x_2[0] = x_1[0]

K_1 = np.empty([N + 1, n, n])
K_1[N] = Q
for i in range(100,0,-1):
    K_1[i-1] = np.matmul(np.matmul(np.transpose(A), K_1[i] - np.matmul(np.matmul(K_1[i], B), np.linalg.inv(np.matmul(np.matmul(np.transpose(B), K_1[i]), B) + R_1), np.matmul(np.transpose(B), K_1[i]))), A) + Q
    
L_1 = np.empty([N, n, n])
for i in range(0, N):
    L_1[i] = - np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(B), K_1[i+1]), B) + R_1), np.matmul(np.matmul(np.transpose(B), K_1[i+1]), A))

K_2 = np.empty([N + 1, n, n])
K_2[N] = Q
for i in range(100,0,-1):
    K_2[i-1] = np.matmul(np.matmul(np.transpose(A), K_2[i] - np.matmul(np.matmul(K_2[i], B), np.linalg.inv(np.matmul(np.matmul(np.transpose(B), K_2[i]), B) + R_2), np.matmul(np.transpose(B), K_2[i]))), A) + Q
    
L_2 = np.empty([N, n, n])
for i in range(0, N):
    L_2[i] = - np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(B), K_2[i+1]), B) + R_2), np.matmul(np.matmul(np.transpose(B), K_2[i+1]), A))

for i in range(0, N):
    # Values normal
    x_1[i+1] = np.matmul(A + np.matmul(B, L_1[i]), x_1[i]) + w[i]
    # Values much larger
    x_2[i+1] = np.matmul(A + np.matmul(B, L_2[i]), x_2[i]) + w[i]

# plot
plot_fig4 = plt.figure(4)
plt.subplot(211)
plt.title('Behavior of the system for two input-cost matrices (R)')
plt.ylabel('first element')
plt.plot(x_1.T[0])
plt.plot(x_2.T[0])
plt.subplot(212)
plt.ylabel('second element')
plt.plot(x_1.T[1], label = 'R "normal"')
plt.plot(x_2.T[1], label = 'R "much larger"')
plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4), ncol=2)

# v --------------------------------------------------
# Fix x[0] and D, R_1[0] = normal - R_2[0] = "much larger"
R = np.diag(np.repeat(1, [n], axis = 0))
x_1 = np.empty([N + 1, n])
x_1[0] = [1, 2]
D = np.diag(np.repeat(1, [n], axis = 0))
w = np.random.multivariate_normal(mean = np.repeat(0, [n], axis = 0), cov = D, size = N)
x_2 = np.empty([N + 1, n])
x_2[0] = x_1[0]

K_1 = np.empty([N + 1, n, n])
K_1[N] = Q
for i in range(100,0,-1):
    K_1[i-1] = np.matmul(np.matmul(np.transpose(A), K_1[i] - np.matmul(np.matmul(K_1[i], B), np.linalg.inv(np.matmul(np.matmul(np.transpose(B), K_1[i]), B) + R), np.matmul(np.transpose(B), K_1[i]))), A) + Q
    
L_1 = np.empty([N, n, n])
for i in range(0, N):
    L_1[i] = - np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(B), K_1[i+1]), B) + R), np.matmul(np.matmul(np.transpose(B), K_1[i+1]), A))

K_2, G, E = control.dare(A = A, B = B, Q = Q, R = R)

L_2 = - np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(B), K_2), B) + R), np.matmul(np.matmul(np.transpose(B), K_2), A))

for i in range(0, N):
    # Values normal
    x_1[i+1] = np.matmul(A + np.matmul(B, L_1[i]), x_1[i]) + w[i]
    # Values much larger
    x_2[i+1] = np.matmul(A + np.matmul(B, L_2), x_2[i]) + w[i]

# plot
plot_fig5 = plt.figure(5)
plt.subplot(211)
plt.title('Behavior of the system optimal control vs steady-state control')
plt.ylabel('first element')
plt.plot(x_1.T[0])
plt.plot(x_2.T[0])
plt.subplot(212)
plt.ylabel('second element')
plt.plot(x_1.T[1], label = 'Optimal control')
plt.plot(x_2.T[1], label = 'Steady-state control')
plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4), ncol=2)

pp = PdfPages('/home/chpmoreno/Dropbox/Documents/BGSE/Second_Term/SMO/Problemsets/PS4/figures.pdf')
pp.savefig(plot_fig2)
pp.savefig(plot_fig3)
pp.savefig(plot_fig4)
pp.savefig(plot_fig5)
pp.close()
