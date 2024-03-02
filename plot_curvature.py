import sys
import numpy as np
import matplotlib.pyplot as plt
from astroquery.sdss import SDSS
from astropy import units as u
from astropy.coordinates import SkyCoord

def plot_sky_patch(ra_min, ra_max, dec_min, dec_max, scale_distance=False):
    """
    Plot a specified patch of the sky using redshift data to indicate hypothetical spacetime curvature
    and optionally scale point size by redshift as a proxy for distance.

    Parameters:
    - ra_min, ra_max: Minimum and maximum right ascension in degrees.
    - dec_min, dec_max: Minimum and maximum declination in degrees.
    - scale_distance: If True, scale point sizes by redshift.
    """
    # Query SDSS for galaxies within the specified RA and Dec range
    sky_coord = SkyCoord(ra=np.mean([ra_min, ra_max])*u.degree, dec=np.mean([dec_min, dec_max])*u.degree, frame='icrs')
    radius = np.max([ra_max - ra_min, dec_max - dec_min]) * u.degree
    query = SDSS.query_region(sky_coord, radius=radius, spectro=True, photoobj_fields=['ra', 'dec', 'z'])

    # If no galaxies are found, return
    if query is None:
        print("No galaxies found in the specified region.")
        return

    # ensure the redshifts are all positive
    # print the min, average, max, and median redshifts
    print("Min redshift:", np.min(query['z']))
    print("Average redshift:", np.mean(query['z']))
    print("Max redshift:", np.max(query['z']))
    print("Median redshift:", np.median(query['z']))
    # print the number of galaxies before removing any
    print("Number of galaxies before removing any:", len(query))
    # remove any data with redshifts less than 0.01
    query = query[query['z'] > 0.01]
    # print the number of galaxies remaining
    print("Number of galaxies:", len(query))

    # Convert RA and Dec to Mollweide projection
    ra_rad = np.radians(query['ra'] - 180)  # Shift RA range to [-180, 180] for Mollweide projection
    dec_rad = np.radians(query['dec'])
    x, y = ra_rad, dec_rad

    # Compute point sizes based on redshift, if scale_distance is True
    point_sizes = (query['z'] / np.max(query['z']) * 100) ** 2 if scale_distance else 10

    # Plot the Mollweide projection
    plt.figure(figsize=(10, 5))
    ax = plt.subplot(111, projection='mollweide')
    sc = ax.scatter(x, y, c=query['z'], cmap='viridis', s=point_sizes, alpha=0.7, edgecolors='none')
    plt.colorbar(sc, orientation='horizontal', pad=0.07, label='Redshift')
    plt.title('Map of Hypothetical Spacetime Curvature')

    # Set grid and labels
    ax.grid(True)

    # Define the tick positions and labels for RA and Dec
    ra_ticks = np.linspace(-np.pi, np.pi, 8)  # Evenly spaced around the circle
    ra_tick_labels = ['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h']
    ax.set_xticks(ra_ticks)
    ax.set_xticklabels(ra_tick_labels)

    dec_ticks = np.linspace(-np.pi/2, np.pi/2, 5)  # Evenly spaced from -90 to +90 degrees
    dec_tick_labels = ['−90°', '−45°', '0°', '45°', '90°']
    ax.set_yticks(dec_ticks)
    ax.set_yticklabels(dec_tick_labels)

    plt.show()

def main():
    # use the arguments to specify the patch of the sky to plot, and if no
    # arguments are provided, plot the entire sky
    if len(sys.argv) == 5:
        ra_min, ra_max, dec_min, dec_max = map(float, sys.argv[1:])
        plot_sky_patch(ra_min, ra_max, dec_min, dec_max, scale_distance=False)
    else:
        plot_sky_patch(0, 360, -90, 90, scale_distance=False)

if __name__ == '__main__':
    main()
