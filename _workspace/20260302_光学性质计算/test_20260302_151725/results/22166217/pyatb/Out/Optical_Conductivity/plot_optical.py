import numpy as np
import matplotlib.pyplot as plt

direction = {
    'xx' : 1,
    'xy' : 2,
    'xz' : 3,
    'yx' : 4,
    'yy' : 5,
    'yz' : 6,
    'zx' : 7,
    'zy' : 8,
    'zz' : 9
}

oc_real_part = np.loadtxt('/home/input_lbg-2570537-22166217/pyatb/Out/Optical_Conductivity/optical_conductivity_real_part.dat')
oc_imag_part = np.loadtxt('/home/input_lbg-2570537-22166217/pyatb/Out/Optical_Conductivity/optical_conductivity_imag_part.dat')
df_real_part = np.loadtxt('/home/input_lbg-2570537-22166217/pyatb/Out/Optical_Conductivity/dielectric_function_real_part.dat')
df_imag_part = np.loadtxt('/home/input_lbg-2570537-22166217/pyatb/Out/Optical_Conductivity/dielectric_function_imag_part.dat')

x = oc_real_part[:, 0]

for key, value in direction.items():
    figure = plt.figure()
    plt.title('Optical conductivity')
    plt.xlim(x[0], x[-1])
    plt.xlabel('$\omega (eV)$')
    plt.ylabel('$\sigma_{%s} (S/m)$'%(key))
    real = oc_real_part[:, value]
    imag = oc_imag_part[:, value]
    plt.plot(x, real, label='real part', color='b', linewidth=1, linestyle='-')
    plt.plot(x, imag, label='imag part', color='r', linewidth=1, linestyle='-')

    plt.legend()
    plt.savefig('/home/input_lbg-2570537-22166217/pyatb/Out/Optical_Conductivity/' + 'oc-%s.pdf'%(key))
    plt.close('all')

for key, value in direction.items():
    figure = plt.figure()
    plt.title('dielectric function')
    plt.xlim(x[0], x[-1])
    plt.xlabel('$\omega (eV)$')
    plt.ylabel('$\epsilon_{%s}$'%(key))
    real = df_real_part[:, value]
    imag = df_imag_part[:, value]
    plt.plot(x, real, label='real part', color='b', linewidth=1, linestyle='-')
    plt.plot(x, imag, label='imag part', color='r', linewidth=1, linestyle='-')

    plt.legend()
    plt.savefig('/home/input_lbg-2570537-22166217/pyatb/Out/Optical_Conductivity/' + 'df-%s.pdf'%(key))
    plt.close('all')

