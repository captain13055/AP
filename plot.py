import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', family = 'Times New Roman')
plt.rc('text', usetex = True)


def BL():
    df = np.genfromtxt('BLF.csv', dtype = str, delimiter = ',')
    dd = np.genfromtxt('BLD.csv', dtype = str, delimiter = ',')

    rf = np.asarray(map(lambda x: 100 * float(x), df[:, 0]))
    rd = np.asarray(map(lambda x: 100 * float(x), dd[:, 0]))
    kf = np.asarray(map(lambda x: float(x), df[:, 1]))
    kd = np.asarray(map(lambda x: float(x), dd[:, 1]))


    plt.subplot(121)
    plt.plot(rf, kf, 'o',
             color = '#00adfe')
    plt.ylim(4, 10)
    plt.xlabel('$r$ [cm]')
    plt.ylabel('$local k + 1$ [arbs]')
    plt.title('F')
    plt.subplot(122)
    plt.plot(rd, kd, 'o',
             color = '#ffae13')
    plt.ylim(6, 15)
    plt.xlabel('$r$ [cm]')
    plt.ylabel('$local k + 1$ [arbs]')
    plt.title('D')
    plt.savefig('fig/FD_BL.pdf')
    plt.show()
def tune():
    f = open('tuneShiftBeta.csv', 'r')
    f.readline()
    d = np.genfromtxt(f, dtype = str, delimiter = ',')
    f.close()

    E = np.asarray(map(lambda x: float(x), d[:, 0]))
    Qh = np.asarray(map(lambda x: float(x), d[:, 1]))
    Qv = np.asarray(map(lambda x: float(x), d[:, 2]))
    r = np.asarray(map(lambda x: float(x) * 100, d[:, 3]))
    

    plt.plot(Qh, Qv, 'x',
             color = '#ff2b82')
    plt.xlabel('$\\nu_h$')
    plt.ylabel('$\\nu_v$')
    plt.title('Tune Diagram')
    plt.savefig('fig/tuneDiagram.pdf')
    plt.show()
    plt.subplot(121)
    plt.plot(E, Qh, 'o',
             color = '#00adfe')
    plt.xlabel('Energy [MeV]')
    plt.ylabel('$\\nu_h$')
    plt.subplot(122)
    plt.plot(E, Qv, 'o',
             color = '#ffae13')
    plt.xlabel('Energy [MeV]')
    plt.ylabel('$\\nu_v$')
    plt.savefig('fig/E_Q.pdf')
    plt.show()
    plt.plot(E, r, '-',
             color = '#0089ff')
    plt.xlabel('Energy')
    plt.ylabel('$r$ at the center of F [cm]')
    plt.savefig('fig/E_r.pdf')
    plt.show()
    
