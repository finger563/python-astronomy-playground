from astroquery.sdss import SDSS
from astropy import coordinates as coords
import numpy as np
import matplotlib.pyplot as plt

speed_of_light = 299792.458  # Speed of light in km/s

def get_galaxy_data(ra, dec, radius=0.1, max_galaxies=100):
    """Query SDSS for galaxies in a given region and return their redshifts and magnitudes.

    Parameters:
    - ra: Right Ascension of the center of the search region (in degrees).
    - dec: Declination of the center of the search region (in degrees).
    - radius: Radius of the search region (in degrees).
    - max_galaxies: Maximum number of galaxies to return.

    Returns:
    - redshifts: Array of redshifts for the galaxies.
    - magnitudes: Array of r-band magnitudes for the galaxies.
    - speeds: Array of relative speeds for the galaxies (in km/s).

    Note: The relative speed is calculated from the redshift using the formula:
    speed = redshift * speed_of_light

    The r-band magnitude is a measure of the brightness of an astronomical
    object as observed through a filter that allows light with wavelengths in
    the "r-band" range to pass through. The r-band is part of the optical
    spectrum and is centered around a wavelength of approximately 620
    nanometers. Magnitudes are a logarithmic scale, and the r-band magnitude
    specifically refers to the brightness of an object as measured in this
    wavelength range.
    """
    sky_coord = coords.SkyCoord(ra, dec, unit="deg")
    query = SDSS.query_region(sky_coord, radius=coords.Angle(radius, unit='deg'), spectro=True, photoobj_fields=['ra', 'dec', 'z', 'petroMag_r'])
    # query = SDSS.query_region(sky_coord, radius=radius, spectro=True, photoobj_fields=['ra', 'dec', 'z', 'petroMag_r'])

    if query is None:
        return np.array([]), np.array([])

    redshifts = query['z'][:max_galaxies]
    magnitudes = query['petroMag_r'][:max_galaxies]
    speeds = redshifts * speed_of_light
    return redshifts, magnitudes, speeds

def plot_redshift_magnitude(redshifts, magnitudes):
    """
    Plot a scatter plot of redshift vs. magnitude for a set of galaxies.

    Parameters:
    - redshifts: Array of redshifts for the galaxies.
    - magnitudes: Array of r-band magnitudes for the galaxies.
    """
    plt.figure(figsize=(8, 6))
    plt.scatter(redshifts, magnitudes, alpha=0.5)
    plt.xlabel('Redshift')
    plt.ylabel('r-band Magnitude')
    plt.title('Redshift vs. Magnitude')
    plt.gca().invert_yaxis()  # Magnitudes are brighter for lower values
    plt.show()

def plot_redshift_magnitude_speed(redshifts, magnitudes, speeds):
    """
    Plot a scatter plot of redshift vs. magnitude and redshift vs. speed for a set of galaxies.

    Parameters:
    - redshifts: Array of redshifts for the galaxies.
    - magnitudes: Array of r-band magnitudes for the galaxies.
    - speeds: Array of relative speeds for the galaxies (in km/s).
    """
    fig, axs = plt.subplots(2, 1, figsize=(8, 12))

    axs[0].scatter(redshifts, magnitudes, alpha=0.5)
    axs[0].set_xlabel('Redshift')
    axs[0].set_ylabel('r-band Magnitude')
    axs[0].set_title('Redshift vs. Magnitude')
    axs[0].invert_yaxis()  # Magnitudes are brighter for lower values

    axs[1].scatter(redshifts, speeds, alpha=0.5, color='r')
    axs[1].set_xlabel('Redshift')
    axs[1].set_ylabel('Relative Speed (km/s)')
    axs[1].set_title('Redshift vs. Relative Speed')

    plt.tight_layout()
    plt.show()

def main():
    ra, dec = 180, 0  # Example coordinates (in degrees)
    radius = 1  # Search radius (in degrees)
    redshifts, magnitudes, speeds = get_galaxy_data(ra, dec, radius)
    if len(redshifts) > 0:
        plot_redshift_magnitude_speed(redshifts, magnitudes, speeds)
    else:
        print("No galaxies found in the specified region.")

if __name__ == '__main__':
    main()
