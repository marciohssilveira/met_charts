import datetime as dt

from siphon.catalog import TDSCatalog


class GetGFSData:

    def __init__(self, variables):
        self.variables = variables
        self.URL = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml'
        self.dataset = 'Latest Collection for GFS Quarter Degree Forecast'

    def get_data(self):
        """
        :param coordinates: tuple like (lon, lat)
        :param variables: chosen list of variables based on the variables list for the dataset
        :param n_hours: number of hours for the prediction
        :return: a subset of the netCDF4 dataset based on the given coordinates and variables
        """
        ###########################################
        # First we construct a TDSCatalog instance using the url
        gfs_catalog = TDSCatalog(self.URL)

        # We see this catalog contains three datasets.
        # print(gfs_catalog.datasets)

        gfs_subset = gfs_catalog.datasets[self.dataset].subset()

        ###########################################
        # Define sub_point to proceed with the query
        query = gfs_subset.query()

        ###########################################
        # Then we construct a query asking for data corresponding to desired latitude and longitude and
        # for the time interval. We also ask for NetCDF version 4 data and choose the variables.
        # This request will return all vertical levels for a single point and for the time interval.
        # Note the string representation of the query is a properly encoded query string.
        # lonlat_box(west, east, south, north)
        query.lonlat_box(north=0, south=-60, east=-30, west=-80)
        now = dt.datetime.utcnow()
        n_hours = 30
        query.time_range(now, now + dt.timedelta(hours=n_hours))
        query.accept('netcdf4')

        ###########################################
        # We'll be pulling out the variables we want to use in the future,
        # as well as the values of pressure levels.
        # To get the name of the correct variable we look at the 'variables' attribute on.
        # The last of the variables listed in `coordinates` is the pressure dimension.
        # print(gfs_subset.variables)

        query.variables(*self.variables)

        ###########################################
        # We now request data from the server using this query.
        print(f'Downloading data...')
        raw_data = gfs_subset.get_data(query)
        return raw_data
