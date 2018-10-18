#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


# create a new 'Legacy VTK Reader'
for x in range(1,101):
	FileList = []
	for n in range(1,7):
		t = (n-1)*2000
		FileNow = '/home/bjacob/Desktop/SSA_CORRECTED/examples/ssa_tsdpd/cylinder_annihilation/temporal_analysis/65x1/%d/dump%d.vtk' % (x,t)
		FileList.append(FileNow)
	
        dump = LegacyVTKReader(FileNames = FileList)

        # get animation scene
        animationScene1 = GetAnimationScene()

        # update animation scene based on data timesteps
        animationScene1.UpdateAnimationUsingDataTimeSteps()

        # get active view
        renderView1 = GetActiveViewOrCreate('RenderView')

        # create a new 'Point Volume Interpolator'
        pointVolumeInterpolator1 = PointVolumeInterpolator(Input=dump,Source='Bounded Volume')
        pointVolumeInterpolator1.Kernel = 'VoronoiKernel'
        pointVolumeInterpolator1.Locator = 'Static Point Locator'

        # init the 'Bounded Volume' selected for 'Source'
        pointVolumeInterpolator1.Source.Origin = [0.0, 0.0, 0.0]
        pointVolumeInterpolator1.Source.Scale = [1.0, 0.0, 0.0]

	# create a new 'Plot Over Line'
        plotOverLine1 = PlotOverLine(Input=pointVolumeInterpolator1,Source='High Resolution Line Source')

        # init the 'High Resolution Line Source' selected for 'Source'
        plotOverLine1.Source.Point1 = [0.0, 0.0, 0.0]
        plotOverLine1.Source.Point2 = [1.0, 0.0, 0.0]

        # Properties modified on plotOverLine1.Source
        plotOverLine1.Source.Resolution = 64

        # Properties modified on plotOverLine1
        plotOverLine1.Tolerance = 2.22044604925031e-16

        # Properties modified on plotOverLine1.Source
        plotOverLine1.Source.Resolution = 64

        # save data
        SaveData('/home/bjacob/Desktop/SSA_CORRECTED/examples/ssa_tsdpd/cylinder_annihilation/temporal_analysis/65x1/%d/cylinder_annihilation_results.csv' % (x), proxy=plotOverLine1, Precision=16,WriteAllTimeSteps=1)

        # destroy plotOverLine1
        Delete(plotOverLine1)
        del plotOverLine1

        # destroy pointVolumeInterpolator1
        Delete(pointVolumeInterpolator1)
        del pointVolumeInterpolator1

        # destroy dump
        Delete(dump)
        del dump
