import numpy as np


def refine_spectrum(spectrum: np.ndarray) -> np.ndarray:
    """
    Refines the spectrum by rounding the eigenvalues to one decimal place and adding 0.25 to the eigenvalues.
    Tries to circumvent the problem of the eigenvalues being slightly off from the actual eigenvalues.

    Args:
        spectrum (np.ndarray): barcode spectrum

    Returns:
        np.ndarray: refined barcode spectrum
    """

    for i, eps in enumerate(spectrum):

        deceps = np.abs(eps - int(eps))

        dev = deceps - 0.5

        if eps < 0.0 and dev < 0.0 and np.abs(dev) < 1e-2:

            spectrum[i] = int(eps) - 0.75

        elif eps > 0.0 and dev < 0.0 and np.abs(dev) < 1e-2:

            spectrum[i] = int(eps) + 0.75

    return spectrum
