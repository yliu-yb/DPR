import numpy
import numpy as np
import datetime
import h5py
import draw

# DPR decode and draw
# GPM_2ADPR GPM DPR Precipitation Profile L2A 1.5 hours 5 km V07

def region_ind(lon, lat, extents, form='mask'):
    '''
    返回落入给定经纬度方框范围内的索引.

    Parameters
    ----------
    lon : ndarray
        经度数组. 若form='mask'则要求形状与lat一致.

    lat : ndarray
        纬度数组. 若form='mask'则要求形状与lon一致.

    extents : 4-tuple of float
        经纬度方框的范围[lonmin, lonmax, latmin, latmax].

    form : {'mask', 'ix'}
        索引的形式.
        'mask': 使用与lon和lat同形状的布尔数组进行索引.
        'ix': 使用下标构成的开放网格进行索引.

    Returns
    -------
    ind : ndarray or 2-tuple of ndarray
        索引的布尔数组或下标数组构成的元组.

    Examples
    --------
    mask = region_ind(lon2d, lat2d, extents)
    data1d = data2d[mask]

    ixgrid = region_ind(lon1d, lat1d, extents, form='ix')
    data2d_subset = data2d[ixgrid]
    data3d_subset = data3d[:, ixgrid[0], ixgrid[1]]
    '''
    lonmin, lonmax, latmin, latmax = extents
    mask_lon = (lon >= lonmin) & (lon <= lonmax)
    mask_lat = (lat >= latmin) & (lat <= latmax)
    if form == 'mask':
        if lon.shape != lat.shape:
            raise ValueError('lon和lat的形状不匹配.')
        ind = mask_lon & mask_lat
    elif form == 'ix':
        ind = np.ix_(mask_lon, mask_lat)
    else:
        raise ValueError('form不支持')

    return ind

# Get DPR data
class DPR():
    def __init__(self, file):
        # 打开HDF文件
        self.hdf_file = h5py.File(file)

    def close(self):
        '''关闭文件.'''
        self.hdf_file.close()

    def __exit__(self, *args):
        self.close()

    # CR
    def get_Cross_CR_Ref(self, box='false'):
        if box == 'false':
            pass
        else:
            # lon lat
            lon = self.hdf_file['FS']['Longitude'][:]
            lat = self.hdf_file['FS']['Latitude'][:]
            # height = self.hdf_file['FS']['PRE']['height']

            # datetime
            year = self.hdf_file['FS']['ScanTime']['Year'][:]
            month = self.hdf_file['FS']['ScanTime']['Month'][:]
            day = self.hdf_file['FS']['ScanTime']['DayOfMonth'][:]
            hour = self.hdf_file['FS']['ScanTime']['Hour'][:]
            minute = self.hdf_file['FS']['ScanTime']['Minute'][:]
            dt = []
            for y, m, d, h, min in zip(year, month, day, hour, minute):
                dt.append([datetime.datetime(year=y, month=m, day=d, hour=h, minute=min)] * 49)
            # ref data
            ref = self.hdf_file['FS']['SLV']['zFactorFinal'][:].astype(np.float64)
            ref_0 = ref[:, :, :, 0]
            ref_1 = ref[:, :, :, 1]

            ref_0 = ref_0[:, :, ::-1]
            ref_1 = ref_1[:, :, ::-1]
            ref_cr_0 = np.max(ref_0, axis=2)
            ref_cr_1 = np.max(ref_1, axis=2)

            ref_cr_0[ref_cr_0 < 0] = np.nan
            ref_cr_1[ref_cr_1 < 0] = np.nan

            # 根据范围裁剪数据
            mask = region_ind(lon, lat, box)
            dt = numpy.array(dt)
            lon = lon[mask]
            lat = lat[mask]
            ref_cr_0 = ref_cr_0[mask]
            ref_cr_1 = ref_cr_1[mask]
            dt = dt[mask]

            return ref_cr_0, ref_cr_1, lon, lat, dt

    # 截面 CAPPI
    def get_Cross_Ref_lon_lat_height_time_At(self, h_idx, box = 'false'):
        if h_idx >= 176:
            raise ValueError('0 <= h_idx < 176')
        if box == 'false':
            pass
        else:
            # lon lat
            lon = self.hdf_file['FS']['Longitude'][:]
            lat = self.hdf_file['FS']['Latitude'][:]
            # height = self.hdf_file['FS']['PRE']['height']

            # datetime
            year = self.hdf_file['FS']['ScanTime']['Year'][:]
            month = self.hdf_file['FS']['ScanTime']['Month'][:]
            day = self.hdf_file['FS']['ScanTime']['DayOfMonth'][:]
            hour = self.hdf_file['FS']['ScanTime']['Hour'][:]
            minute = self.hdf_file['FS']['ScanTime']['Minute'][:]
            dt = []
            for y, m, d, h, min in zip(year, month, day, hour, minute):
                dt.append([datetime.datetime(year=y, month=m,day=d,hour=h,minute=min)] * 49)
            # ref data
            ref = self.hdf_file['FS']['SLV']['zFactorFinal'][:].astype(np.float64)
            ref_0 = ref[:,:,:, 0]
            ref_1 = ref[:,:,:, 1]

            ref_0 = ref_0[:,:,::-1]
            ref_1 = ref_1[:,:,::-1]

            ref_0 = ref_0[:,:,h_idx]
            ref_1 = ref_1[:,:,h_idx]

            ref_0[ref_0 < 0] = np.nan
            ref_1[ref_1 < 0] = np.nan

            # ref_0 = np.ma.masked_where(ref_0 < 0, ref_0)
            # ref_1 = np.ma.masked_where(ref_1 < 0, ref_1)

            # 根据范围裁剪数据
            mask = region_ind(lon, lat, box)
            dt = numpy.array(dt)
            lon = lon[mask]
            lat = lat[mask]
            ref_0 = ref_0[mask]
            ref_1 = ref_1[mask]
            dt = dt[mask]
            height = (h_idx + 1) * 125
            # print(numpy.min(dt), numpy.max(dt))

            return ref_0, ref_1, lon, lat, height, dt

    # 廓线
    def get_Profile_Ref_At(self, ray_idx, box = 'false'):
        if ray_idx >= 49:
            raise ValueError('0 <= h_idx < 49')
        if box == 'false':
            pass
        else:
            # lon lat
            lon = self.hdf_file['FS']['Longitude'][:, ray_idx]
            lat = self.hdf_file['FS']['Latitude'][:, ray_idx]
            # datetime
            year = self.hdf_file['FS']['ScanTime']['Year'][:]
            month = self.hdf_file['FS']['ScanTime']['Month'][:]
            day = self.hdf_file['FS']['ScanTime']['DayOfMonth'][:]
            hour = self.hdf_file['FS']['ScanTime']['Hour'][:]
            minute = self.hdf_file['FS']['ScanTime']['Minute'][:]
            dt = []
            for y, m, d, h, min in zip(year, month, day, hour, minute):
                dt.append([datetime.datetime(year=y, month=m,day=d,hour=h,minute=min)])
            # ref data
            ref = self.hdf_file['FS']['SLV']['zFactorFinal'][:, ray_idx, :, :].astype(np.float64)

            ref_0 = ref[:, :, 0]
            ref_1 = ref[:, :, 1]

            ref_0 = ref_0[:, ::-1]
            ref_1 = ref_1[:, ::-1]

            ref_0[ref_0 < 0] = np.nan
            ref_1[ref_1 < 0] = np.nan

            # ref_0 = np.ma.masked_where(ref_0 < 0, ref_0)
            # ref_1 = np.ma.masked_where(ref_1 < 0, ref_1)
            #

            height = [h * 0.125 for h in range(176)]

            # other var
            # temp
            airT = self.hdf_file['FS']['VER']['airTemperature'][:, ray_idx, :].astype(np.float64)
            airT = airT[:, ::-1]
            airT = numpy.array(airT) - 273.15
            # phase
            phase = self.hdf_file['FS']['DSD']['phase'][:, ray_idx, :].astype(np.float64)
            phase = phase[:, ::-1]
            # 找到数组中不为255的元素的索引
            phase[phase == 255] = np.nan
            # height of zero deg
            heightZeroDeg = self.hdf_file['FS']['VER']['heightZeroDeg'][:, ray_idx]
            heightBB = self.hdf_file['FS']['CSF']['heightBB'][:, ray_idx]
            precipWaterIntegrated_0 = self.hdf_file['FS']['SLV']['precipWaterIntegrated'][:, ray_idx, 0]
            precipWaterIntegrated_1 = self.hdf_file['FS']['SLV']['precipWaterIntegrated'][:, ray_idx, 1]

            # surface rainrate
            rain_rate_surface = self.hdf_file['FS']['SLV']['precipRateNearSurface'][:, ray_idx]
            # average of precipitation rate for 2 to 4km height.
            precipRateAve24 = self.hdf_file['FS']['SLV']['precipRateAve24'][:, ray_idx]

            # DSD
            # Parameters of DSD functions, Nw and Dm. Nw in 1/m3 mm
            # ref = self.hdf_file['FS']['SLV']['zFactorFinal'][:, ray_idx, :, :].astype(np.float64)

            Nw = self.hdf_file['FS']['SLV']['paramDSD'][:, ray_idx, :, 0]
            Dm = self.hdf_file['FS']['SLV']['paramDSD'][:, ray_idx, :, 1]
            Nw = Nw[:, ::-1]
            Dm = Dm[:, ::-1]
            Nw[Nw < 0] = np.nan
            Dm[Dm < 0] = np.nan
            heightBB[heightBB <= 0] = np.nan
            # 根据范围裁剪数据
            mask = region_ind(lon, lat, box)
            dt = numpy.array(dt)
            lon = lon[mask]
            lat = lat[mask]
            ref_0 = ref_0[mask]
            ref_1 = ref_1[mask]
            airT = airT[mask]
            Nw = Nw[mask]
            Dm = Dm[mask]
            phase = phase[mask]
            # phase[phase < 0] = np.nan

            heightZeroDeg = heightZeroDeg[mask]
            heightBB = heightBB[mask]
            rain_rate_surface = rain_rate_surface[mask]
            precipRateAve24 = precipRateAve24[mask]
            precipWaterIntegrated_0 = precipWaterIntegrated_0[mask]
            precipWaterIntegrated_1 = precipWaterIntegrated_1[mask]
            dt = dt[mask]
            # print(numpy.min(dt), numpy.max(dt))
            return ref_0, ref_1, lon, lat, height, dt, heightZeroDeg, heightBB, \
                   rain_rate_surface, precipRateAve24, precipWaterIntegrated_0, precipWaterIntegrated_1, \
                   airT, Nw, Dm, phase
        pass

    def get_Profile_Lon_at(self, ray_idx):
        pass

    def get_Profile_Lat_at(self, ray_idx):
        pass

    def get_Profile_height(self):
        pass

if __name__ == "__main__":
    file_name = '2A.GPM.DPR.V9-20211125.20230908-S003033-E020302.054129.V07B.HDF5'
    DPR = DPR("./" + file_name)
    box1 = [100, 130, 0, 30]
    ref_0_cr, ref_1_cr, lon_cr, lat_cr, dt_cr = DPR.get_Cross_CR_Ref(box1)

    # ref_0_cappi, ref_1_cappi, lon_cappi, lat_cappi, height_cappi, dt_cappi = DPR.get_Cross_Ref_lon_lat_height_time_At(39, box)

    box2 = [100, 130, 18, 24]

    ref_0_profile, ref_1_profile, lon_profile, lat_profile, height_profile, dt_profile, heightZeroDeg, heightBB, \
                   rain_rate_surface, precipRateAve24, precipWaterIntegrated_0, precipWaterIntegrated_1, \
                   airT, Nw, Dm, phase = DPR.get_Profile_Ref_At(10, box2)

    if len(dt_profile) == 0:
        print("*"*5, "-1", "*"*5)
        exit()

    draw.draw(ref_0_cr, lon_cr, lat_cr, dt_cr,
              ref_0_profile, lon_profile, lat_profile, height_profile, dt_profile,
              heightZeroDeg, heightBB,
              rain_rate_surface, precipRateAve24, precipWaterIntegrated_0,
              airT, Nw, Dm, phase, file_name, 'KuPR', box1
              )

    # draw.draw(ref_0_cappi, lon_cappi, lat_cappi, height_cappi, dt_cappi,
    #           ref_0_profile, lon_profile, lat_profile, height_profile, dt_profile,
    #           heightZeroDeg, heightBB,
    #           rain_rate_surface, precipRateAve24, precipWaterIntegrated_0,
    #           airT, Nw, Dm, phase, file_name, 'KuPR', box
    #             )

    # draw.draw(ref_1_cappi, lon_cappi, lat_cappi, height_cappi, dt_cappi,
    #           ref_1_profile, lon_profile, lat_profile, height_profile, dt_profile,
    #           heightZeroDeg, heightBB,
    #           rain_rate_surface, precipRateAve24, precipWaterIntegrated_0,
    #           airT, Nw, Dm, phase, file_name, 'KaPR'
    #           )
