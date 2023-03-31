import numpy as np
import matplotlib.pyplot as plt

def read(filename_v, filename_t):
    # Read v and t from file
    with open(filename_v) as f:
        v = np.array([float(x) for x in f.readlines()])

    with open(filename_t) as f:
        t = np.array([float(x) for x in f.readlines()])

    return v, t

def plot(v, t):
    # Plot v vs t
    plt.figure(figsize=(10, 5))
    plt.plot(t, v)
    plt.xlabel('t (ms)')
    plt.ylabel('V (mV)')
    plt.show()


# Min model
v, t = read('V.txt', 't.txt')
plot(v, t)
