import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib as mpl
import cartopy.crs as ccrs
import numpy
from matplotlib import cm as mat_cm
import frykit.plot as fplt
import numpy as np
from matplotlib import colors
import os
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
import cartopy.io.shapereader as shpreader

mpl.rcParams['font.family'] = 'Arial'

# ref
idx_REF = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
ref_bar_ticks = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

# ref color
rgb_value = 255 / float(len(idx_REF) + 1)
ref_color = []
for i in range(1, len(idx_REF) + 1):
    ref_color.append(mat_cm.jet(int(i * rgb_value)))
ref_cmap = (colors.ListedColormap(ref_color)
            .with_extremes(over='#8B008B'))
ref_norm = mpl.colors.BoundaryNorm(idx_REF, ref_cmap.N)

phase_colors = [
    'blue', 'gold', 'green'
]
phase_ticklabels = [
    'solid', 'mixed', 'liquid']
phase_cmap, phase_norm, phase_ticks = fplt.make_qualitative_cmap(phase_colors)

# 定义自定义colormap的颜色映射
# phase_cmap = ListedColormap(phase_colors)


cm = 1 / 2.54  # centimeters in inches
crs = ccrs.PlateCarree()


def draw(ref_cappi, lon_cappi, lat_cappi, height_cappi, dt_cappi,
              ref_profile, lon_profile, lat_profile, height_profile, dt_profile,
              heightZeroDeg, heightBB,
              rain_rate_surface,
              airT, Nw, Dm,
              file_name):
    figsize = [8 * cm, 8 * cm]
    fig = plt.figure(figsize=figsize)
    # 绘制地图.
    extents = [100, 130, 0, 30]
    ax1 = fig.add_axes([0.1, 0.4, 0.9, 0.9], projection=crs)
    ax1.coastlines(resolution='10m', lw=0.5)
    fplt.set_extent_and_ticks(
        ax1, extents=extents,
        yticks=np.arange(0, 31, 5),
        xticks=np.arange(100, 131, 5),
        nx=1, ny=1
    )
    ax1.tick_params(labelsize='medium')

    # shapefile
    filepath = 'D:/code/web/data_download/V2.0/Shp/Re_XA/修复后雄安县界.shp'
    reader = shpreader.Reader(filepath)
    geoms = reader.geometries()
    print(geoms)
    exit()
    ax1.add_geometries(geoms, crs, lw=0.5, fc='none')
    reader.close()

    # CAPPI
    ax1.scatter(lon_cappi, lat_cappi, color = 'grey')
    ax1.scatter(lon_cappi, lat_cappi, c=ref_cappi, norm = ref_norm, cmap = ref_cmap, marker = 's',s = .1, transform=crs)

    ax1.plot(lon_profile, lat_profile, '-', color='black', lw='1.5', transform=ccrs.Geodetic())
    date_str = numpy.min(dt_profile).strftime('%Y-%m-%d')
    ax1.set_title('(a) ' + date_str + ' (UTC) \n' + ' '*6 +numpy.min(dt_cappi).strftime('%H:%M:%S~') +
                  numpy.max(dt_cappi).strftime('%H:%M:%S')
                  , loc='left', fontsize='large')

    ax1.set_title('CAPPI: ' + str(height_cappi * 0.001) + 'km', loc='right', fontsize='large')

    ax2 = fplt.add_side_axes(ax1, loc='right', pad=0.2, depth=0.9)
    # 构造截面图所需的x轴刻度.
    x, xticks, xticklabels = fplt.get_slice_xticks(
        lon_profile, lat_profile, ntick=5, decimals=1
    )
    X, Y = np.meshgrid(x, height_profile, indexing='ij')

    pc = ax2.pcolormesh(X, Y, ref_profile, norm = ref_norm, cmap = ref_cmap)

    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticklabels, fontsize='medium')
    ax2.set_ylim(0, 20)
    ax2.set_ylabel('Height (km)', fontsize='medium')
    ax2.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax2.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    # reflectivity
    ax2.set_title('(b) ' + 'Reflectivity', loc='left', fontsize='large')

    ax2.plot(heightZeroDeg * 0.001, ls='--', c='r', lw=1, label='$0\degree $C height')
    # ax2.plot(heightBB * 0.001, ls='--', c='black', lw=1, label='BB height')
    ax2.scatter(x, heightBB * 0.001, marker = 'x', c='black', label='brightband')

    ax2.legend(loc='upper right', fontsize='large')

    # 设置colorbar.
    cax12 = fplt.add_side_axes([ax1, ax2], loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax12, ticks=ref_bar_ticks, orientation='horizontal', label="dBZ")
    cbar.set_label(label="dBZ", fontsize='large')
    #
    # Dm
    ax3 = fplt.add_side_axes(ax2, loc='right', pad=0.2, depth=0.9)
    pc = ax3.pcolormesh(X, Y, Dm, cmap = 'jet')
    ax3.set_xticks(xticks)
    ax3.set_xticklabels(xticklabels, fontsize='medium')
    ax3.set_ylim(0, 20)
    ax3.set_ylabel('Height (km)', fontsize='medium')
    ax3.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax3.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    ax3.set_title('(c) ' + 'Parameters of DSD functions, Dm', loc='left', fontsize='large')
    cax3 = fplt.add_side_axes(ax3, loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax3, orientation='horizontal', label="mm")
    cbar.set_label(label="mm", fontsize='large')

    # Nw
    Nw = 10 ** (numpy.array(Nw) * 0.1)
    ax4 = fplt.add_side_axes(ax3, loc='right', pad=0.2, depth=0.9)
    pc = ax4.pcolormesh(X, Y, Nw, cmap = 'jet')
    ax4.set_xticks(xticks)
    ax4.set_xticklabels(xticklabels, fontsize='medium')
    ax4.set_ylim(0, 20)
    ax4.set_ylabel('Height (km)', fontsize='medium')
    ax4.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax4.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    cax4 = fplt.add_side_axes(ax4, loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax4, orientation='horizontal', label="m$^{-3}$ mm")
    ax4.set_title('(d) ' + 'Parameters of DSD functions, Nw', loc='left', fontsize='large')
    cbar.set_label(label="m$^{-3}$ mm", fontsize='large')

    # temp
    ax5 = fplt.add_side_axes(ax4, loc='right', pad=0.2, depth=0.9)
    pc = ax5.pcolormesh(X, Y, airT, cmap = 'jet')
    ax5.set_xticks(xticks)
    ax5.set_xticklabels(xticklabels, fontsize='medium')
    ax5.set_ylim(0, 20)
    ax5.set_ylabel('Height (km)', fontsize='medium')
    ax5.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax5.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    cax5 = fplt.add_side_axes(ax5, loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax5, orientation='horizontal', label="$\degree$C")
    ax5.set_title('(e) ' + 'Air Temperature', loc='left', fontsize='large')
    cbar.set_label(label="$\degree$C", fontsize='large')

    # 绘制九段线
    with open(r'E:\doctor\china-geospatial-data-GB2312\ten-dash-line.gmt') as src:
        context = ''.join([line for line in src if not line.startswith('#')])
        blocks = [cnt for cnt in context.split('>') if len(cnt) > 0]
        borders = [np.fromstring(block, dtype=float, sep=' ') for block in blocks]
    for line in borders:
        ax1.plot(line[0::2], line[1::2], '-', color='black', lw='1', transform=ccrs.Geodetic())

    plt.savefig(file_name + ".png", dpi=900, bbox_inches='tight')

# new fig
def draw(ref_cappi, lon_cappi, lat_cappi, height_cappi, dt_cappi,
              ref_profile, lon_profile, lat_profile, height_profile, dt_profile,
              heightZeroDeg, heightBB,
              rain_rate_surface, precipRateAve24, precipWaterIntegrated_0,
              airT, Nw, Dm, phase,
              file_name,
              type,
              box):
    figsize = [8 * cm, 8 * cm]
    fig = plt.figure(figsize=figsize)
    # 绘制地图.
    extents = box
    ax1 = fig.add_axes([0.1, 0.4, 0.9, 0.9], projection=crs)
    ax1.coastlines(resolution='10m', lw=0.5)
    fplt.set_extent_and_ticks(
        ax1, extents=extents,
        yticks=np.arange(box[2], box[3] + 1, 5),
        xticks=np.arange(box[0], box[1] + 1, 5),
        nx=1, ny=1
    )
    ax1.tick_params(labelsize='medium')

    # shapefile
    filepath = 'D:/code/web/data_download/V2.0/Shp/Re_XA/修复后雄安县界.shp'
    reader = shpreader.Reader(filepath)
    geoms = reader.geometries()
    # print(geoms)
    # exit()
    ax1.add_geometries(geoms, crs, lw=0.5, fc='none')
    reader.close()

    # CAPPI
    ax1.scatter(lon_cappi, lat_cappi, color = 'grey')
    ax1.scatter(lon_cappi, lat_cappi, c=ref_cappi, norm = ref_norm, cmap = ref_cmap, marker = 's',s = .1, transform=crs)

    ax1.plot(lon_profile, lat_profile, '-', color='black', lw='1.5', transform=ccrs.Geodetic())
    date_str = numpy.min(dt_profile).strftime('%Y-%m-%d')
    ax1.set_title('(a) ' + date_str + ' (UTC) \n' + ' '*6 +numpy.min(dt_cappi).strftime('%H:%M:%S~') +
                  numpy.max(dt_cappi).strftime('%H:%M:%S')
                  , loc='left', fontsize='large')

    ax1.set_title('CAPPI: ' + str(height_cappi * 0.001) + 'km', loc='right', fontsize='large')

    ax2 = fplt.add_side_axes(ax1, loc='right', pad=0.2, depth=0.9)
    # 构造截面图所需的x轴刻度.
    x, xticks, xticklabels = fplt.get_slice_xticks(
        lon_profile, lat_profile, ntick=5, decimals=1
    )
    X, Y = np.meshgrid(x, height_profile, indexing='ij')

    pc = ax2.pcolormesh(X, Y, ref_profile, norm = ref_norm, cmap = ref_cmap)

    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticklabels, fontsize='medium')
    ax2.set_ylim(0, 20)
    ax2.set_ylabel('Height (km)', fontsize='medium')
    ax2.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax2.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    # reflectivity
    ax2.set_title('(b) ' + 'Reflectivity', loc='left', fontsize='large')

    # 0层度高度
    ax2.plot(heightZeroDeg * 0.001, ls='--', c='r', lw=1, label='$0\degree $C height')
    # 亮带
    ax2.plot(heightBB * 0.001, ls='--', c='black', lw=1, label='brightband')
    # ax2.scatter(x, heightBB * 0.001, marker = 'x', c='black', label='brightband')

    ax2.legend(loc='upper right', fontsize='large')

    # 设置colorbar.
    cax12 = fplt.add_side_axes([ax1, ax2], loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax12, ticks=ref_bar_ticks, orientation='horizontal', label="dBZ")
    cbar.set_label(label="dBZ", fontsize='large')
    #
    # Dm
    ax3 = fplt.add_side_axes(ax2, loc='right', pad=0.2, depth=0.9)
    pc = ax3.pcolormesh(X, Y, Dm, cmap = 'jet')
    ax3.set_xticks(xticks)
    ax3.set_xticklabels(xticklabels, fontsize='medium')
    ax3.set_ylim(0, 20)
    ax3.set_ylabel('Height (km)', fontsize='medium')
    ax3.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax3.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    ax3.set_title('(c) ' + 'Parameters of DSD functions, Dm', loc='left', fontsize='large')
    cax3 = fplt.add_side_axes(ax3, loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax3, orientation='horizontal', label="mm")
    cbar.set_label(label="mm", fontsize='large')

    # Nw
    Nw = 10 ** (numpy.array(Nw) * 0.1)
    ax4 = fplt.add_side_axes(ax3, loc='right', pad=0.2, depth=0.9)
    pc = ax4.pcolormesh(X, Y, Nw, cmap = 'jet', norm=mcolors.LogNorm())
    ax4.set_xticks(xticks)
    ax4.set_xticklabels(xticklabels, fontsize='medium')
    ax4.set_ylim(0, 20)
    ax4.set_ylabel('Height (km)', fontsize='medium')
    ax4.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax4.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    cax4 = fplt.add_side_axes(ax4, loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax4, orientation='horizontal', label="m$^{-3}$ mm", format='%g')
    ax4.set_title('(d) ' + 'Parameters of DSD functions, Nw', loc='left', fontsize='large')
    cbar.set_label(label="m$^{-3}$ mm", fontsize='large')

    # temp
    ax5 = fplt.add_side_axes(ax4, loc='right', pad=0.2, depth=0.9)
    pc = ax5.pcolormesh(X, Y, airT, cmap = 'jet')
    ax5.set_xticks(xticks)
    ax5.set_xticklabels(xticklabels, fontsize='medium')
    ax5.set_ylim(0, 20)
    ax5.set_ylabel('Height (km)', fontsize='medium')
    ax5.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax5.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    cax5 = fplt.add_side_axes(ax5, loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax5, orientation='horizontal', label="$\degree$C")
    ax5.set_title('(e) ' + 'Air Temperature', loc='left', fontsize='large')
    cbar.set_label(label="$\degree$C", fontsize='large')

    # precipitation line
    rain_rate_surface = numpy.ma.masked_where(rain_rate_surface < 0, rain_rate_surface)
    precipRateAve24 = numpy.ma.masked_where(precipRateAve24 < 0, precipRateAve24)
    precipWaterIntegrated_0 = numpy.ma.masked_where(precipWaterIntegrated_0 < 0, precipWaterIntegrated_0)

    ax6 = fplt.add_side_axes(cax12, loc='bottom', pad=0.2, depth=0.9)
    # rain_rate_surface, precipRateAve24, precipWaterIntegrated_0,
    ax6.plot(x, rain_rate_surface, marker = '+', color = '#E2AE79', label = 'precipRateNearSurface')
    ax6.plot(x, precipRateAve24, marker = '+', color = '#D0DCAA', label = 'precipRateAve(2~4 km)')
    ax6.set_xticks(xticks)
    ax6.set_xticklabels(xticklabels, fontsize='medium')
    ax6.set_ylabel('mm hr$^{-1}$', color='black')

    ax66 = ax6.twinx()
    ax66.plot(x, precipWaterIntegrated_0, marker = '+', color = '#9B3A4D', label = 'precipWaterIntegrated')
    ax66.set_ylabel('g m$^{-2}$', color='black')
    ax6.set_title('(f) ' + 'Precipitation', loc='left', fontsize='large')

    # 合并左边y轴和右边y轴的图例
    lines, labels = ax6.get_legend_handles_labels()
    lines2, labels2 = ax66.get_legend_handles_labels()
    ax6.legend(lines + lines2, labels + labels2, loc='upper right')

    # phase
    # ax7 = fplt.add_side_axes(ax6, loc='right', pad=0.35, depth=0.9)
    # pc = ax7.pcolormesh(X, Y, phase, cmap=phase_cmap, norm=phase_norm)
    # ax7.set_xticks(xticks)
    # ax7.set_xticklabels(xticklabels, fontsize='medium')
    # ax7.set_ylim(0, 20)
    # ax7.set_ylabel('Height (km)', fontsize='medium')
    # ax7.yaxis.set_major_locator(mticker.MultipleLocator(5))
    # ax7.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    # cax7 = fplt.add_side_axes(ax7, loc='bottom', pad=0.15, depth=0.03)
    # cbar = fig.colorbar(pc, cax=cax7, orientation='horizontal', label="")
    #
    # # 创建colorbar，并设置标签
    # cbar.set_ticks(phase_ticks)
    # cbar.set_ticklabels(phase_ticklabels)
    # ax7.set_title('(g) ' + 'Phase', loc='left', fontsize='large')


    # 绘制九段线
    with open(r'E:\doctor\china-geospatial-data-GB2312\ten-dash-line.gmt') as src:
        context = ''.join([line for line in src if not line.startswith('#')])
        blocks = [cnt for cnt in context.split('>') if len(cnt) > 0]
        borders = [np.fromstring(block, dtype=float, sep=' ') for block in blocks]
    for line in borders:
        ax1.plot(line[0::2], line[1::2], '-', color='black', lw='1', transform=ccrs.Geodetic())

    save_folder = './fig'
    if not os.path.exists(save_folder):
        # 如果目录不存在，创建目录
        os.makedirs(save_folder)

    plt.savefig(save_folder + '/' + file_name + "_"+type+".png", dpi=900, bbox_inches='tight')

# cr
def draw(ref_cr, lon_cr, lat_cr, dt_cr,
              ref_profile, lon_profile, lat_profile, height_profile, dt_profile,
              heightZeroDeg, heightBB,
              rain_rate_surface, precipRateAve24, precipWaterIntegrated_0,
              airT, Nw, Dm, phase,
              file_name,
              type,
              box1):
    figsize = [8 * cm, 8 * cm]
    fig = plt.figure(figsize=figsize)
    # 绘制地图.
    extents = box1
    ax1 = fig.add_axes([0.1, 0.4, 0.9, 0.9], projection=crs)
    ax1.coastlines(resolution='10m', lw=0.5)
    fplt.set_extent_and_ticks(
        ax1, extents=extents,
        yticks=np.arange(box1[2], box1[3] + 1, 5),
        xticks=np.arange(box1[0], box1[1] + 1, 5),
        nx=1, ny=1
    )
    ax1.tick_params(labelsize='medium')

    # CAPPI
    ax1.scatter(lon_cr, lat_cr, color = 'grey')
    ax1.scatter(lon_cr, lat_cr, c=ref_cr, norm = ref_norm, cmap = ref_cmap, marker = 's',s = .1, transform=crs)

    ax1.plot(lon_profile, lat_profile, '-', color='black', lw='1.5', transform=ccrs.Geodetic())
    date_str = numpy.min(dt_profile).strftime('%Y-%m-%d')
    ax1.set_title('(a) ' + date_str + ' (UTC) \n' + ' '*6 +numpy.min(dt_cr).strftime('%H:%M:%S~') +
                  numpy.max(dt_cr).strftime('%H:%M:%S')
                  , loc='left', fontsize='large')

    ax1.set_title('CR', loc='right', fontsize='large')

    ax2 = fplt.add_side_axes(ax1, loc='right', pad=0.2, depth=0.9)
    # 构造截面图所需的x轴刻度.
    x, xticks, xticklabels = fplt.get_slice_xticks(
        lon_profile, lat_profile, ntick=5, decimals=1
    )
    X, Y = np.meshgrid(x, height_profile, indexing='ij')

    pc = ax2.pcolormesh(X, Y, ref_profile, norm = ref_norm, cmap = ref_cmap)

    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticklabels, fontsize='medium')
    ax2.set_ylim(0, 20)
    ax2.set_ylabel('Height (km)', fontsize='medium')
    ax2.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax2.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    # reflectivity
    ax2.set_title('(b) ' + 'Reflectivity', loc='left', fontsize='large')

    # 0层度高度
    ax2.plot(heightZeroDeg * 0.001, ls='--', c='r', lw=1, label='$0\degree $C height')
    # 亮带
    ax2.plot(heightBB * 0.001, ls='--', c='black', lw=1, label='brightband')
    # ax2.scatter(x, heightBB * 0.001, marker = 'x', c='black', label='brightband')

    ax2.legend(loc='upper right', fontsize='large')

    # 设置colorbar.
    cax12 = fplt.add_side_axes([ax1, ax2], loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax12, ticks=ref_bar_ticks, orientation='horizontal', label="dBZ")
    cbar.set_label(label="dBZ", fontsize='large')
    #
    # Dm
    ax3 = fplt.add_side_axes(ax2, loc='right', pad=0.2, depth=0.9)
    pc = ax3.pcolormesh(X, Y, Dm, cmap = 'jet')
    ax3.set_xticks(xticks)
    ax3.set_xticklabels(xticklabels, fontsize='medium')
    ax3.set_ylim(0, 20)
    ax3.set_ylabel('Height (km)', fontsize='medium')
    ax3.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax3.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    ax3.set_title('(c) ' + 'Parameters of DSD functions, Dm', loc='left', fontsize='large')
    cax3 = fplt.add_side_axes(ax3, loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax3, orientation='horizontal', label="mm")
    cbar.set_label(label="mm", fontsize='large')

    # Nw
    Nw = 10 ** (numpy.array(Nw) * 0.1)

    normNw = mcolors.LogNorm(vmin=1E2, vmax=1E5)

    ax4 = fplt.add_side_axes(ax3, loc='right', pad=0.2, depth=0.9)
    pc = ax4.pcolormesh(X, Y, Nw, cmap = 'jet', norm=normNw)
    ax4.set_xticks(xticks)
    ax4.set_xticklabels(xticklabels, fontsize='medium')
    ax4.set_ylim(0, 20)
    ax4.set_ylabel('Height (km)', fontsize='medium')
    ax4.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax4.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    cax4 = fplt.add_side_axes(ax4, loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax4, orientation='horizontal', label="m$^{-3}$ mm")
    ax4.set_title('(d) ' + 'Parameters of DSD functions, Nw', loc='left', fontsize='large')
    cbar.set_label(label="m$^{-3}$ mm", fontsize='large')

    # temp
    ax5 = fplt.add_side_axes(ax4, loc='right', pad=0.2, depth=0.9)
    pc = ax5.pcolormesh(X, Y, airT, cmap = 'jet')
    ax5.set_xticks(xticks)
    ax5.set_xticklabels(xticklabels, fontsize='medium')
    ax5.set_ylim(0, 20)
    ax5.set_ylabel('Height (km)', fontsize='medium')
    ax5.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax5.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    cax5 = fplt.add_side_axes(ax5, loc='bottom', pad=0.15, depth=0.03)
    cbar = fig.colorbar(pc, cax=cax5, orientation='horizontal', label="$\degree$C")
    ax5.set_title('(e) ' + 'Air Temperature', loc='left', fontsize='large')
    cbar.set_label(label="$\degree$C", fontsize='large')

    # precipitation line
    rain_rate_surface = numpy.ma.masked_where(rain_rate_surface < 0, rain_rate_surface)
    precipRateAve24 = numpy.ma.masked_where(precipRateAve24 < 0, precipRateAve24)
    precipWaterIntegrated_0 = numpy.ma.masked_where(precipWaterIntegrated_0 < 0, precipWaterIntegrated_0)

    ax6 = fplt.add_side_axes(cax12, loc='bottom', pad=0.2, depth=0.9)
    # rain_rate_surface, precipRateAve24, precipWaterIntegrated_0,
    ax6.plot(x, rain_rate_surface, marker = '+', color = '#E2AE79', label = 'precipRateNearSurface')
    ax6.plot(x, precipRateAve24, marker = '+', color = '#D0DCAA', label = 'precipRateAve(2~4 km)')
    ax6.set_xticks(xticks)
    ax6.set_xticklabels(xticklabels, fontsize='medium')
    ax6.set_ylabel('mm hr$^{-1}$', color='black')

    ax66 = ax6.twinx()
    ax66.plot(x, precipWaterIntegrated_0, marker = '+', color = '#9B3A4D', label = 'precipWaterIntegrated')
    ax66.set_ylabel('g m$^{-2}$', color='black')
    ax6.set_title('(f) ' + 'Precipitation', loc='left', fontsize='large')

    # 合并左边y轴和右边y轴的图例
    lines, labels = ax6.get_legend_handles_labels()
    lines2, labels2 = ax66.get_legend_handles_labels()
    ax6.legend(lines + lines2, labels + labels2, loc='upper right')


    # 绘制九段线
    with open('./ten-dash-line.gmt') as src:
        context = ''.join([line for line in src if not line.startswith('#')])
        blocks = [cnt for cnt in context.split('>') if len(cnt) > 0]
        borders = [np.fromstring(block, dtype=float, sep=' ') for block in blocks]
    for line in borders:
        ax1.plot(line[0::2], line[1::2], '-', color='black', lw='1', transform=ccrs.Geodetic())

    save_folder = './fig'
    if not os.path.exists(save_folder):
        # 如果目录不存在，创建目录
        os.makedirs(save_folder)

    plt.savefig(save_folder + '/' + file_name + "_"+type+".png", dpi=300, bbox_inches='tight')
    plt.close()


