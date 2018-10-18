#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


# create a new 'Legacy VTK Reader'
for x in range(1,101):
	FileList = []
	for n in range(1,101):
		t = (n-1)*100
		FileNow = '/home/bjacob/Desktop/SSA_CORRECTED/examples/ssa_tsdpd/diffusion1d/17x1/%d/dump%d.vtk' % (x,t)
		FileList.append(FileNow)
	
 	dump = LegacyVTKReader(FileNames = FileList)

	# get animation scene
	animationScene1 = GetAnimationScene()

	# update animation scene based on data timesteps
	animationScene1.UpdateAnimationUsingDataTimeSteps()

	# get active view
	renderView1 = GetActiveViewOrCreate('RenderView')
	# uncomment following to set a specific view size
	# renderView1.ViewSize = [1715, 1128]

	# show data in view
	dumpDisplay = Show(dump, renderView1)
	# trace defaults for the display properties.
	dumpDisplay.ColorArrayName = [None, '']
	dumpDisplay.OSPRayScaleArray = 'C_[0]'
	dumpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
	dumpDisplay.GlyphType = 'Arrow'

	# reset view to fit data
	renderView1.ResetCamera()

	#changing interaction mode based on data extents
	renderView1.InteractionMode = '2D'
	renderView1.CameraPosition = [0.5, 10000.0, 10000.0]
	renderView1.CameraFocalPoint = [0.5, 0.0, 0.0]
	renderView1.CameraViewUp = [1.0, 1.0, 0.0]

	# save data
	SaveData('/home/bjacob/Desktop/SSA_CORRECTED/examples/ssa_tsdpd/diffusion1d/17x1/%d/diffusion1d_results.csv' % (x), proxy=dump, Precision=16,
		WriteAllTimeSteps=1)

	#### saving camera placements for all active views

	# current camera placement for renderView1
	renderView1.InteractionMode = '2D'
	renderView1.CameraPosition = [0.5, 10000.0, 10000.0]
	renderView1.CameraFocalPoint = [0.5, 0.0, 0.0]
	renderView1.CameraViewUp = [0.7071067811865475, 0.7071067811865475, 0.0]
	renderView1.CameraParallelScale = 0.5

	#### uncomment the following to render all views
	# RenderAllViews()
	# alternatively, if you want to write images, you can use SaveScreenshot(...).

