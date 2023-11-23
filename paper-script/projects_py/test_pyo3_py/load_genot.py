

import numpy as np
import genetics_pyo3


def load_genot_py():
    # None
    use_samples = np.array([], dtype=bool)
    use_snvs = np.array([True, False, True], dtype=bool)
    genot = genetics_pyo3.load_genot("./test/data/toy3/genot",
                                     # "./test/data/toy1/genot.phe", "phe",
                                     use_samples,
                                     use_snvs,
                                     True)
                                     #False)
    print('genot', genot)


def axpy():
    a = 0.1
    x = np.array([1, 2, 3], dtype=np.float64)
    y = np.array([4, 5, 6], dtype=np.float64)
    z = genetics_pyo3.axpy(a, x, y)

    print('z', z)


def main():
    # axpy()
    load_genot_py()


if __name__ == '__main__':
    main()
