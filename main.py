
from dbfread import DBF
import csv
import tables
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd
from geopandas import GeoDataFrame
import geodatasets
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import shapely as shp
from shapely.geometry import Polygon
import numpy as np
import matplotlib.cm as cm

# Open the DBF file
table = DBF('pireps_202409010000_202410262359.dbf', encoding='utf-8')
# Open the hdf file
hdf_file = 'MOD11C3.A2022001.061.2022032100901.hdf'
csv_file = 'temperature.csv'

"""
# Write to CSV
with open('turbLatLonOnly.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    # Write header
    writer.writerow(['VALID', 'TURB', 'LAT', 'LON', 'color'])
    # Write records
    print("CSV made...")
    i = 0
    for record in table:
        obj = type('', (), {})
        obj.time = record['VALID']
        obj.turbulence = record['TURB']
        obj.latitude = record['LAT']
        obj.longitude = record['LON']
        # print("\nRecord", i)
        # print("Time", obj.time)
        # print("Turbulence", obj.turbulence)
        # print("Latitude", obj.latitude)
        # print("Longitude", obj.longitude)
        if i == 0:
            print([obj.time, obj.turbulence, obj.latitude, obj.longitude])
        i += 1
        if obj.turbulence is not "":
            if "LGT" in obj.turbulence:
                color = "g"
            elif "MOD" in obj.turbulence:
                color = "y"
            elif "NEG" in obj.turbulence:
                color = "b"
            elif "SEV" in obj.turbulence:
                color = "r"
            else:
                color = "k"
            writer.writerow([obj.time, obj.turbulence, obj.latitude, obj.longitude, color])
    # reader = csv.DictReader(csvfile)
    # for row in reader:
"""
df = pd.read_csv("turbLatLonOnly.csv", delimiter=',', skiprows=0, low_memory=False)
#print(df['color'])
geometry = [Point(xy) for xy in zip(df['LON'], df['LAT'])]
gdf = GeoDataFrame(df, geometry=geometry)   

#this is a simple map that goes with geopandas
# deprecated: world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

world = gpd.read_file(geodatasets.data.naturalearth.land['url'])
gdf.plot(ax=world.plot(figsize=(10, 6)), marker='o', color=df['color'], markersize=1);

plt.show()

# make values for a typical coordinate grid with a specified step
step = 1
lons_full = np.arange(-179.125, 179.125 + step, step)
lats_full = np.arange(89.125, -89.125 - step, -step) # latitudes range from 90 to -90 to match array row indices growing from top to bottom

# get values for a special grid, representing central points of cells in the normal grid
halfstep = step/2
lons = lons_full[:-1] + halfstep
lats = lats_full[:-1] - halfstep

# make a 2D array to keep square polygons corresponding to those grid cells
squares = [[None for col in range(len(lons))] for row in range(len(lats))]
squares_shape = np.asarray(squares).shape

# make shapely polygons and store them in that array
for row, col in np.ndindex(squares_shape):
    # longitudes correspond to columns, latitudes to rows
    lon = lons[col]
    lat = lats[row]
    # define vertices of a square polygon
    bottom_left = (lon - halfstep, lat - halfstep)
    bottom_right = (lon + halfstep, lat - halfstep)
    top_right = (lon + halfstep, lat + halfstep)
    top_left = (lon - halfstep, lat + halfstep)
    coords = (bottom_left, bottom_right, top_right, top_left, bottom_left)
    # make a polygon
    squares[row][col] = shp.Polygon(coords)

routes_file = 'temperature_change_2022.csv'
routes_gdf = gpd.read_file(routes_file)

df = pd.read_csv("temperature_change_2022.csv", delimiter=',', skiprows=0, low_memory=False)
#print(df['color'])
geometry = [Point(xy) for xy in zip(df['LON'], df['LAT'])]
gdf = GeoDataFrame(df, geometry=geometry)   
                
# get coastlines and countries
countries_file = gpd.datasets.get_path('naturalearth_lowres')
countries_gdf = gpd.read_file(countries_file)

# show everything on one plot
fig, ax = plt.subplots(figsize=(15, 8))
plt.title('Distribuzione del traffico marittimo globale')
# show coastlines and countries
countries_gdf.plot(ax=ax, facecolor='none', edgecolor='grey')

# make bounds for 10 sections of a colour bar
cmap = cm.Reds
trips_min = -1
trips_max = +5
bounds = np.linspace(trips_min, trips_max + 1, 11)
norm = mcolors.BoundaryNorm(bounds, ncolors=cmap.N)
# show squares, whose colour represents the trip count
im = ax.imshow(gdf, cmap=cmap, alpha=0.5, extent=(min(lons_full), max(lons_full), min(lats_full), max(lats_full)), norm=norm)
# add a color bar
cb = fig.colorbar(im)
cb.set_label('Numero di trip_count unici per cella')

plt.show()