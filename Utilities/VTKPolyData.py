#!/usr/bin/env python

"""
Description: Render a list of VTK data  and (or) an nifti image, then view or save PNG or save WebGL.
The code is modified from ITKExamples.

When an output PNG or WebGL file is not specified, an interactive windows is
displayed.  To get the camera position from the interactive window, press the
"c" key.

e: exit the application
q: quit the application

s: surface representation
w: wireframe representation

Examples:
VTKPolyData.py --vtk file1.vtk file2.vtk --image im.nii.gz
VTKPolyData.py --vtk file1.vtk file2.vtk --image im.nii.gz --sliceX 50 --sliceY 40 --sliceZ 60
VTKPolyData.py --vtk file1.vtk file2.vtk --image im.nii.gz --png out.png

Author(s): Jian Cheng (jian.cheng.1983@gmail.com)
"""

import argparse
import sys

import vtk



def arg_values(value, typefunc, numberOfValues):
    '''set arguments based using comma. If numberOfValues<0, it supports arbitrary number of inputs.'''
    values = value.split(',')
    if numberOfValues>0 and len(values) != numberOfValues:
        raise argparse.ArgumentError
    return map(typefunc, values)

def arg_bool(parser, boolarg, defaultbool, helpdoc):
    '''set bool argment'''
    parser.add_argument('--' + boolarg, dest=boolarg, action='store_true', help='with ' + boolarg +'. '+ helpdoc)
    parser.add_argument('--no-' + boolarg, dest=boolarg, action='store_false', help='without ' + boolarg +'. '+ helpdoc)
    exec ''.join(['parser.set_defaults(', boolarg, '=', 'True' if defaultbool else 'False', ')'])

def ConvertDataSetToSurface(algorithmOutputPort):
    dataSetSurfaceFilter = vtk.vtkDataSetSurfaceFilter()
    dataSetSurfaceFilter.SetInputConnection(algorithmOutputPort)
    dataSetSurfaceFilter.UseStripsOn()
    dataSetSurfaceFilter.Update()
    polyData = vtk.vtkPolyData()
    polyData.ShallowCopy(dataSetSurfaceFilter.GetOutput())
    return polyData

def ReadLegacyVTK(file_name):
    '''Support POLYDATA and some other vtk data (*.vtk).'''
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(file_name)
    reader.Update()
    if None != reader.GetPolyDataOutput():
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(reader.GetPolyDataOutput())
        return polyData, 'PolyData'
    if None != reader.GetUnstructuredGridOutput():
        return ConvertDataSetToSurface(reader.GetOutputPort()), 'UnstructuredGrid'
    if None != reader.GetStructuredPointsOutput():
        return ConvertDataSetToSurface(reader.GetOutputPort()), 'StructuredPoints'
    if None != reader.GetStructuredGridOutput():
        return ConvertDataSetToSurface(reader.GetOutputPort()), 'StructuredGrid'
    if None != reader.GetRectilinearGridOutput():
        return ConvertDataSetToSurface(reader.GetOutputPort()), 'RectilinearGrid'
    else:
        raise Exception("Unsupported type!\n")

def PrintImage(image):
    '''print image'''
    dim = image.GetDimensions()
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                if image.GetScalarComponentAsFloat(i, j, k, 0)>0:
                    print i, j, k, image.GetScalarComponentAsFloat(i, j, k, 0)

def FlipVTKImageData(image, image_flip):
    '''flip image for visualization.'''

    if image_flip[0]==1 and image_flip[1]==1 and image_flip[2]==1:
        return image

    image2 = vtk.vtkImageData()
    image2.DeepCopy(image)

    for i in range(len(image_flip)):
        if image_flip[i]==-1:
            flipFilter = vtk.vtkImageFlip()
            flipFilter.SetInputData(image2)
            flipFilter.SetFilteredAxis(i)
            flipFilter.Update()
            image2 = flipFilter.GetOutput()

    return image2


def main():

    #  work for arguments with minus sign
    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

    parser = argparse.ArgumentParser(description=__doc__,  formatter_class=argparse.RawTextHelpFormatter)

    # vtk data
    parser.add_argument('--vtk', help='.vtk PolyData input', nargs='*')
    arg_bool(parser, 'frame', False, 'Visualize frame.')
    arg_bool(parser, 'surface', True, 'Visualize surface.')
    arg_bool(parser, 'normal', True, 'Use vtkPolyDataNormals for polydata visualization.')

    # nifti data
    parser.add_argument('--image', help='nifti image file', nargs=1)
    parser.add_argument('--image-flip', help='flip x,y,z axis in image. Default: (1,1,1)', default=(1,1,1),
                        type=(lambda value: arg_values(value, int, 3)), metavar=('flipZ,flipY,flipZ') )
    arg_bool(parser, 'interpolate', False, 'Image interpolation.')
    parser.add_argument('--sliceX', help='x-axis slices of the image', type=(lambda value: arg_values(value, int, -1)), metavar=('x1,x2,...'))
    parser.add_argument('--sliceY', help='y-axis slices of the image', type=(lambda value: arg_values(value, int, -1)), metavar=('y1,y2,...'))
    parser.add_argument('--sliceZ', help='z-axis slices of the image', type=(lambda value: arg_values(value, int, -1)), metavar=('z1,z2,...'))
    parser.add_argument('--valuerange', help='lowest and highest contrast value for the image visualization. If not set, use the minimal and maximal values in the image.',
                        type=(lambda value: arg_values(value, float, 2)), metavar=('lowest value, highest value') )

    # camera
    parser.add_argument('--angle', help='azimuth and elevation for camera',
                        type=(lambda value: arg_values(value, float, 2)), metavar=('azimuth,elevation'))
    parser.add_argument('--position', help='Camera Position',
                        type=(lambda value: arg_values(value, float, 3)), metavar=('x,y,z'))
    parser.add_argument('--focal-point', help='Camera FocalPoint',
                        type=(lambda value: arg_values(value, float, 3)), metavar=('x,y,z'))
    parser.add_argument('--view-up', help='Camera ViewUp',
                        type=(lambda value: arg_values(value, float, 3)), metavar=('x,y,z'))
    parser.add_argument('--bgcolor', help='back ground color',
                        type=(lambda value: arg_values(value, float, 3)), metavar=('r,g,b'), default=(0,0,0))
    parser.add_argument('--size', help='Window size in pixels. Default: (600,600)',default=(600,600),
                        type=(lambda value: arg_values(value, int, 2)), metavar=('width,height') )
    parser.add_argument('--clipping-range', help='Window size in pixels',
                        type=(lambda value: arg_values(value, float, 2)), metavar=('near,far'))

    parser.add_argument('--png', help='Output PNG file',
                        metavar='file.png')
    parser.add_argument('--zoom', help='camera zoom factor',
                        type=float, metavar=('zoomfactor'), default=1.0)
    parser.add_argument('--webgl', help='File prefix for WebGL output',
                        metavar='webglFilePrefix')

    args = parser.parse_args()
    # print args

    if not args.vtk and not args.image:
        print "need inputs for --vtk and (or) --image"
        raise argparse.ArgumentError

    render_window = vtk.vtkRenderWindow()
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(args.bgcolor[0],args.bgcolor[1],args.bgcolor[2])
    render_window.AddRenderer(renderer)
    render_window.SetSize(args.size)

    if args.vtk:
        for inputFile in args.vtk:
            polyData, dataType = ReadLegacyVTK(inputFile)

            if args.frame:
                frame_mapper = vtk.vtkDataSetMapper()
                frame_mapper.SetInputData(polyData)
                frame_actor = vtk.vtkLODActor()
                frame_actor.SetMapper(frame_mapper)
                prop = frame_actor.GetProperty()
                prop.SetRepresentationToWireframe()
                prop.SetColor(0.0, 0.0, 1.0)
                renderer.AddActor(frame_actor)


            if args.surface:
                surface_mapper = vtk.vtkDataSetMapper()

                if args.normal and polyData.GetPointData().GetNormals() is None:
                    polyDataNormals = vtk.vtkPolyDataNormals()
                    try:
                        polyDataNormals.SetInputData(polyData)
                    except:
                        polyDataNormals.SetInput(polyData)
                    # polyDataNormals.SetFeatureAngle(90.0)
                    surface_mapper.SetInputConnection(
                        polyDataNormals.GetOutputPort())
                else:
                    try:
                        surface_mapper.SetInputData(polyData)
                    except:
                        surface_mapper.SetInput(polyData)

                surface_actor = vtk.vtkLODActor()
                surface_actor.SetMapper(surface_mapper)
                prop = surface_actor.GetProperty()
                prop.SetRepresentationToSurface()
                renderer.AddActor(surface_actor)


    if args.image:

        imagefile = args.image[0]
        reader = vtk.vtkNIFTIImageReader()
        reader.SetFileName(imagefile)
        reader.Update()

        # set origin based on itk::NiftiImageIO
        im = reader.GetOutput()
        niftiheader = reader.GetNIFTIHeader()
        im.SetOrigin(-niftiheader.GetQOffsetX(), -niftiheader.GetQOffsetY(), niftiheader.GetQOffsetZ())

        #  for 2D image, set sliceX or sliceY or slizeZ automatically
        x1, x2, y1, y2, z1, z2 = im.GetExtent()
        if x1==x2:
            args.sliceX = [x1]
        if y1==y2:
            args.sliceY = [y1]
        if z1==z2:
            args.sliceZ = [z1]

        # for 3D image, slice needs to be set manually
        if not args.sliceX and not args.sliceY and not args.sliceZ:
            raise Exception("need to choose one or more slice for the 3D image!\n")

        # flip image if needed
        image = FlipVTKImageData(im, args.image_flip)

        lut = vtk.vtkLookupTable()
        valueRange = args.valuerange if args.valuerange else image.GetScalarRange()
        lut.SetTableRange(valueRange[0], valueRange[1])
        lut.SetSaturationRange(0, 0)
        lut.SetHueRange(0, 0)
        lut.SetValueRange(0, 1)
        lut.SetRampToLinear()
        lut.Build()

        planeColors = vtk.vtkImageMapToColors()
        planeColors.SetInputData(image)
        planeColors.SetLookupTable(lut)
        planeColors.Update()

        assem = vtk.vtkAssembly()

        if args.sliceX:
            for x in args.sliceX:
                act = vtk.vtkImageActor()
                act.SetInputData(planeColors.GetOutput())
                act.SetDisplayExtent(x, x, y1, y2, z1, z2)
                act.InterpolateOn() if args.interpolate else act.InterpolateOff()
                assem.AddPart(act)

        if args.sliceY:
            for y in args.sliceY:
                act = vtk.vtkImageActor()
                act.SetInputData(planeColors.GetOutput())
                act.SetDisplayExtent(x1, x2, y, y, z1, z2)
                act.InterpolateOn() if args.interpolate else act.InterpolateOff()
                assem.AddPart(act)

        if args.sliceZ:
            for z in args.sliceZ:
                act = vtk.vtkImageActor()
                act.SetInputData(planeColors.GetOutput())
                act.SetDisplayExtent(x1, x2, y1, y2, z, z)
                act.InterpolateOn() if args.interpolate else act.InterpolateOff()
                assem.AddPart(act)

        renderer.AddViewProp(assem)


    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)
    interactorStyle = render_window_interactor.GetInteractorStyle()
    interactorStyle.SetCurrentStyleToTrackballCamera()


    render_window.Render()

    camera = renderer.GetActiveCamera()


    def print_camera_position(interactor, event):
        def cmd_line_friendly(xyz):
            return '[{0:+8.4e},{1:+8.4e},{2:+8.4e}]'.format(*xyz)

        def cmd_line_friendly2(clip_range):
            return '[{0:+8.4e},{1:+8.4e}]'.format(*clip_range)
        if interactor.GetKeySym() == 'c':
            print('\nPosition:    ' + cmd_line_friendly(camera.GetPosition()))
            print('FocalPoint:   ' + cmd_line_friendly(camera.GetFocalPoint()))
            print('ViewUp:       ' + cmd_line_friendly(camera.GetViewUp()))
            print('ClippingRange:' + cmd_line_friendly2(camera.GetClippingRange()))
            sys.stdout.write('\n--position ')
            sys.stdout.write(cmd_line_friendly(camera.GetPosition()))
            sys.stdout.write(' --focal-point ')
            sys.stdout.write(cmd_line_friendly(camera.GetFocalPoint()))
            sys.stdout.write(' --view-up ')
            sys.stdout.write(cmd_line_friendly(camera.GetViewUp()))
            sys.stdout.write(' --clipping-range ')
            sys.stdout.write(cmd_line_friendly2(camera.GetClippingRange()))
            sys.stdout.write('\n')
            sys.stdout.flush()

    render_window_interactor.AddObserver('KeyPressEvent', print_camera_position)

    if args.position:
        camera.SetPosition(args.position)
    if args.focal_point:
        camera.SetFocalPoint(args.focal_point)
    if args.view_up:
        camera.SetViewUp(args.view_up)
    if args.clipping_range:
        camera.SetClippingRange(args.clipping_range)
    if args.angle:
        camera.Roll(args.angle[0])
        camera.Elevation(args.angle[1])


    camera.Zoom(args.zoom)

    # re-render after setting camera
    render_window.Render()

    if args.png:
        window_to_image = vtk.vtkWindowToImageFilter()
        window_to_image.SetInput(render_window)
        window_to_image.SetMagnification(2)
        png_writer = vtk.vtkPNGWriter()
        png_writer.SetInputConnection(window_to_image.GetOutputPort())
        png_writer.SetFileName(args.png)
        png_writer.Write()

    if args.webgl:
        webgl_exporter = vtk.vtkPVWebGLExporter()
        webgl_exporter.SetFileName(args.webgl)
        webgl_exporter.SetRenderWindow(render_window)
        webgl_exporter.Update()

    if not args.png and not args.webgl:
        render_window_interactor.Start()


if __name__ == '__main__':
    main()
