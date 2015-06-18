#!/usr/bin/env python

"""
Description: Render a list of VTK Polydata and view or save PNG or save WebGL. Modified from ITKExamples.

When an output PNG or WebGL file is not specified, an interactive windows is
displayed.  To get the camera position from the interactive window, press the
"c" key.

e: exit the application
q: quit the application

s: surface representation
w: wireframe representation
"""

import argparse
import sys

import vtk


def three_floats(value):
    values = value[1:-1].split(',')
    if len(values) != 3:
        raise argparse.ArgumentError
    return map(float, values)


def two_floats(value):
    values = value[1:-1].split(',')
    if len(values) != 2:
        raise argparse.ArgumentError
    return map(float, values)

def arg_bool(parser, boolarg, defaultbool):
    '''set bool argment'''
    parser.add_argument('--' + boolarg, dest=boolarg, action='store_true', help='with ' + boolarg)
    parser.add_argument('--no-' + boolarg, dest=boolarg, action='store_false', help='without ' + boolarg)
    exec ''.join(['parser.set_defaults(', boolarg, '=', 'True' if defaultbool else 'False', ')'])

def main():

    parser = argparse.ArgumentParser(description=__doc__,  formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('inputFiles', help='.vtk PolyData input', nargs='*')
    arg_bool(parser, 'frame', False)
    arg_bool(parser, 'surface', True)
    arg_bool(parser, 'normal', True)
    parser.add_argument('--position', help='Camera Position',
                        type=three_floats, metavar=('x,y,z'))
    parser.add_argument('--focal-point', help='Camera FocalPoint',
                        type=three_floats, metavar=('x,y,z'))
    parser.add_argument('--view-up', help='Camera ViewUp',
                        type=three_floats, metavar=('x,y,z'))
    parser.add_argument('--bgcolor', help='back ground color',
                        type=float, metavar=('r', 'g', 'b'), nargs=3)
    parser.add_argument('--size', help='Window size in pixels',
                        type=int, metavar=('width', 'height'), nargs=2)
    parser.add_argument('--clipping-range', help='Window size in pixels',
                        type=two_floats, metavar=('near,far'))
    parser.add_argument('--png', help='Output PNG file',
                        metavar='file.png')
    parser.add_argument('--webgl', help='File prefix for WebGL output',
                        metavar='webglFilePrefix')

    args = parser.parse_args()

    if len(args.inputFiles) == 0:
        print "need inputs"
        raise argparse.ArgumentError

    render_window = vtk.vtkRenderWindow()
    renderer = vtk.vtkRenderer()
    if args.bgcolor:
        renderer.SetBackground(args.bgcolor[0],args.bgcolor[1],args.bgcolor[2])
    else:
        renderer.SetBackground(0,0,0)
    render_window.AddRenderer(renderer)
    render_window.SetSize(600, 600)

    for inputFile in args.inputFiles:
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(inputFile)

        if args.size:
            render_window.SetSize(args.size)

        reader.Update()

        if args.frame:
            frame_mapper = vtk.vtkDataSetMapper()
            frame_mapper.SetInputConnection(reader.GetOutputPort())
            frame_actor = vtk.vtkActor()
            frame_actor.SetMapper(frame_mapper)
            prop = frame_actor.GetProperty()
            prop.SetRepresentationToWireframe()
            prop.SetColor(0.0, 0.0, 1.0)
            renderer.AddActor(frame_actor)


        if args.surface:
            surface_mapper = vtk.vtkDataSetMapper()
            polyData = reader.GetOutput()
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

            surface_actor = vtk.vtkActor()
            surface_actor.SetMapper(surface_mapper)
            prop = surface_actor.GetProperty()
            prop.SetRepresentationToSurface()
            renderer.AddActor(surface_actor)

    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

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

    if args.png:
        window_to_image = vtk.vtkWindowToImageFilter()
        window_to_image.SetInput(render_window)
        window_to_image.SetMagnification(2);
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
